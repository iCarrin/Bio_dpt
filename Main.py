from Output import *
from Primer_functions import *
from Multiplex import *
import re # run 'pip install regex' if not already installed
import time # to handle rate limiting
import requests
from typing import List # for type hinting

# Ensembl REST API base URL
ENSEMBL_REST = "https://rest.ensembl.org"

def Fetch_SNP_Data(rsids: List[str], flank_length: int = 800, use_all = False) -> list[dict]:
    """
    Retrieves SNP data from the Ensembl API, including flanking sequences and alleles.

    Args:
        rsids: List of SNP identifiers (e.g., ["rs1799971"]).
        flank_length: Number of base pairs to include on either side of the SNP.

    Returns:
        DataFrame containing SNP ID, allele, modified sequence, and SNP position.
    """
    snp_data = []
    headers = {"Content-Type": "application/json"}  # Required by Ensembl for JSON responses

    for rsid in rsids:
        try:
            # Step 1: Get SNP mapping (chromosome + position info)
            var_resp = requests.get(f"{ENSEMBL_REST}/variation/homo_sapiens/{rsid}?", headers=headers)
            print(f"getting info on {rsid}")
            var_resp.raise_for_status()
            var_data = var_resp.json()

            mappings = var_data.get("mappings", [])
            if not mappings:
                print(f"Warning: No mappings found for {rsid}.")
                continue

            # We'll just take the first mapping (usually sufficient for common SNPs)
            mapping = mappings[0]
            chrom = mapping["seq_region_name"]  # e.g., "11"
            pos = int(mapping["start"])   

            # Extract allele string like "A/G" or "C/T"
            allele_str = mapping.get("allele_string", "")
            # ancestral = mapping.get("ancestral_allele")

            #Isaiah change: I dropped any ancestral Alleles that we knew were normal
            #In the future we should pull the list, split it, and have the user pick which ones they want.
            raw_alleles = allele_str.split("/") if allele_str else []
       

            def ask_user(use_all = False):
                if use_all:
                    alleles_wanted = "ALL"
                else:
                    for i, allele in enumerate(raw_alleles):
                        print(f"{i+1}) {allele}")
                    alleles_wanted = input("type the corresponding numbers for the alleles you want separated by spaces, or just \"All\" for all of them")

                if alleles_wanted.strip().upper() == "ALL" or alleles_wanted.strip().upper() == "":
                    print("using all alleles")
                    return raw_alleles
                elif re.fullmatch(r'[0-9\s]+', alleles_wanted):
                    indices = [int(x) for x in alleles_wanted.strip().split(" ")]

                    if max(indices) > len(raw_alleles) or min(indices) < 1:
                        print("you asked for an allele that wasn't in the list")
                        return ask_user()
                    else:
                        return [raw_alleles[i-1] for i in indices]
                else:
                    print("you typed more than just numbers and spaces. Try again")
                    return ask_user()
                
            wanted_alleles = ask_user(use_all)

            # Ensure there are at least two alleles to work with
            if len(wanted_alleles) < 1:
                print(f"Warning: Less than 2 alleles for {rsid}. Skipping.")
                continue

            # Step 2: Fetch the flanking DNA sequence around the SNP
            seq_start = max(1, pos - flank_length)  # 1-based for Ensemble 
            seq_end = pos + flank_length #might run off the end of the chromosome if very unlucky
            seq_url = f"{ENSEMBL_REST}/sequence/region/human/{chrom}:{seq_start}..{seq_end}:1?"
            # print (f'chrom: {chrom}\nseq_start: {seq_start}\nseq_end: {seq_end}')
            seq_resp = requests.get(seq_url, headers={"Content-Type": "text/plain"})
            
            seq_resp.raise_for_status()
          
            template_seq = seq_resp.text.strip()

            # Position of the SNP relative to the start of the fetched sequence
            rel_pos = flank_length if seq_start > 1 else pos #should just be flanking length
            



            # Step 3: Replace the SNP base with each possible allele to simulate variation
            for allele in wanted_alleles:
                # Validate that allele contains valid DNA characters only
                if not re.fullmatch("[ACGTNacgtn]+", allele):
                    print(f"Skipping non-standard allele '{allele}' for {rsid}")
                    continue

                # Insert allele at the SNP site
                modified_seq = template_seq[:rel_pos] + allele.upper() + template_seq[rel_pos + 1:]

                # Append to results
                snp_data.append({
                    "snpID": rsid,
                    "allele": allele.upper(),
                    "sequence": modified_seq,
                    "position": rel_pos
                })

            # Sleep to respect Ensembl's rate limit (max 15 req/sec)
            time.sleep(0.34)

        except Exception as e:
            print(f"Error processing {rsid}: {e}")

    # If nothing was successfully retrieved, return an empty DataFrame
    if not snp_data:
        print("No valid SNP data could be retrieved.")
        return []
  

    return snp_data
# ran 5 test of 63 snps. 1 snp call takes 12.7 seconds

from datetime import datetime 
 
def Main():
    """
        Main function to generate, filter, rank, pair, and export primers.
        TODO: Add comprehensive error handling and logging.
        - Log progress and errors to a file.
        - Add input validation for rsIDs and output format.
        TODO: Test with real SNP data.
        - Validate output with biological experts.
        - Benchmark performance for large SNP sets.
        """
    start = datetime.now()

    snp_df = Fetch_SNP_Data([
                            "rs1799971", "rs12184297", "rs116801199", "rs12565286", "rs2977670", "rs28454925", "rs116587930", "rs116720794", "rs4951859",\
                            "rs148120343", "rs529266287", "rs79643588", "rs17396518", "rs983166", "rs28842593", "rs7014597", "rs599839", "rs28428499",\
                            "rs756427959", "rs1045833", "rs28503599", "rs531646671", "rs541940975", "rs62635297", "rs6682375", "rs6682385", "rs11807796",\
                            "rs11580262", "rs201026389", "rs11489794", "rs113141985", "rs62636498", "rs148220436", "rs141130360", "rs199745162", "rs544798311",\
                            "rs199740902", "rs201535981", "rs372841554", "rs201057270", "rs369546232", "rs373030326", "rs746313125", "rs768288163", "rs13302934",\
                            "rs773821537", "rs371996741", "rs4099234", "rs62028215", "rs62637791", "rs62637793", "rs200677948", "rs796517829", "rs62637813",\
                            "rs2691277", "rs3091274", "rs10399749", "rs3107975", "rs28396308", "rs13343114", "rs6658003", "rs201418760", "rs199609147",\
                            "rs2691334", "rs2691335", "rs75024357", "rs2531261", "rs76735897", "rs77573425", "rs3844233", "rs62639101", "rs62639102"], 1500, True)
    #                         "rs111440589", "rs574185548", "rs567315672", "rs771510052", "rs749441654", "rs12402780", "rs13328683", "rs201488854" , "rs3964475",\
    #                         "rs200715572", "rs200113991", "rs201609557", "rs796405068", "rs376621245", "rs78627600", "rs201323924", "rs368526587", "rs12409693",\
    #                         "rs62642104", "rs62642105", "rs375749629", "rs374391784", "rs200298631", "rs373871696", "rs62642109", "rs796779951", "rs3115860",\
    #                         "rs62642112", "rs4287120", "rs200280933", "rs201078696", "rs74747225", "rs79073610", "rs200678306", "rs796221635", "rs3131965",\
    #                         "rs368487483", "rs372912307", "rs182468771", "rs1851943", "rs377389630", "rs370064167", "rs79114531", "rs368732900", "rs4951929",\
    #                         "rs199894737", "rs376726476", "rs377161483", "rs372768449", "rs371080566", "rs367730352", "rs369820305", "rs375850923", "rs2286139",\
    #                         "rs371677125", "rs370365546", "rs375833141", "rs373582709", "rs375144581", "rs376555721", "rs147252685", "rs370691115", "rs3115847",\
    #                         "rs371864317", "rs199583388", "rs370188889", "rs10047119", "rs10047174", "rs10047175", "rs10047121", "rs200482415", "rs10047230",\
    #                         "rs112821789", "rs112455420", "rs200066978", "rs8179466", "rs200864529", "rs8179455", "rs8179403", "rs201234721", "rs6603780",\
    #                         "rs6422503", "rs62635272", "rs184796003", "rs768683158", "rs112186068", "rs6603779", "rs6670732", "rs8179401", "rs199966892",\
    #                         "rs184451216", "rs201702841", "rs200000439", "rs202156508", "rs143434963", "rs147990930", "rs79974410", "rs11490246", "rs4970334",\
    #                         "rs201921502", "rs200295906", "rs62635275", "rs140116466", "rs796211217", "rs796317578", "rs112390583", "rs199660279", "rs202020431"\
    #                         "rs200930606", "rs768053647", "rs62635277", "rs62639107", "rs376115996", "rs62639108", "rs202029170", "rs375277195", "rs62639110",\
    #                         "rs372766634", "rs11489922", "rs200683317", "rs200725328", "rs62639117", "rs767396087", "rs62639119", "rs201544877", "rs200542524",\
    #                         "rs554693544", "rs377659925", "rs373366992", "rs756558238", "rs576317820", "rs62428730", "rs199492794", "rs199491011", "rs61767340",\
    #                         "rs192687139", "rs61767343", "rs374895564", "rs373722660", "rs149807327", "rs375649215", "rs371573906", "rs144672821", "rs201044673",\
    #                         "rs201679166", "rs199511579", "rs61769279", "rs4120955", "rs78894760", "rs150334593", "rs137905425", "rs78216996", "rs1578391",\
    #                         "rs538388406", "rs6594029", "rs55668158", "rs7416152", "rs9283151", "rs6421780", "rs78150957", "rs373437560", "rs200850333",\
    #                         "rs6650105", "rs6594033", "rs796410684", "rs6594035", "rs374132538", "rs368347679", "rs369475281", "rs140845213", "rs143332664",\
    #                         "rs148361140", "rs200919600", "rs770009340", "rs200380870", "rs3857304", "rs4030160", "rs61769323", "rs200021096", "rs760939004",\
    #                         "rs199959440", "rs377611181", "rs760310201", "rs202070913", "rs377354243", "rs61769340", "rs201764515", "rs79279787", "rs530588710",\
    #                         "rs371072735", "rs4951925", "rs200073030", "rs201987734", "rs9768594", "rs71221966", "rs200806899", "rs74564050", "rs368813658",\
    #                         "rs773974007", "rs146067153", "rs139961686", "rs377059754", "rs61769348", "rs200318145", "rs200877292", "rs61769353", "rs9727436",\
    #                         "rs72631875", "rs12029736", "rs113462541", "rs372626744", "rs12028261", "rs3131984", "rs4614203", "rs10900602", "rs10751453",\
    #                         "rs3121393", "rs3115846", "rs34882115", "rs374587598", "rs12735556", "rs12735557", "rs62653102", "rs200283906", "rs201300314",\
    #                         "rs199596027", "rs111203391", "rs747843478", "rs200087344", "rs200420414", "rs3131980", "rs762353340", "rs566981575", "rs3131978",\
    #                         "rs3131977", "rs541848210", "rs61770164", "rs61770165", "rs12564807", "rs7518433", "rs761625058", "rs61770166", "rs199800910",\
    #                         "rs4951928", "rs3131974", "rs2427926", "rs2427917", "rs3131973", "rs3094317", "rs552659822", "rs3115865", "rs191696195",\
    #                         "rs2488505", "rs2519055", "rs3925106", "rs2722509", "rs2519043", "rs199734026", "rs4399101", "rs12045157", "rs12409679",\
    #                         "rs12410209", "rs12409695", "rs11240770", "rs11240771", "rs10793768", "rs3094315", "rs3131972", "rs3131971", "rs61770172",\
    #                         "rs3131970", "rs2073814", "rs59507194", "rs58324164", "rs3131969", "rs3131968", "rs3131967", "rs3115859", "rs3131966",\
    #                         "rs3131964", "rs3115858", "rs3115857", "rs12567639", "rs61768167", "rs61768168", "rs3131963", "rs3131962", "rs3115853",\
    #                         "rs70949520", "rs4951862", "rs3131956", "rs3131955", "rs3131954", "rs3115852", "rs3115851", "rs1048488", "rs3115850",\
    #                         "rs1057213", "rs1064272", "rs1057212", "rs59038458", "rs3115849", "rs3115848", "rs3131950", "rs3131949", "rs3131948",\
    #                         "rs372986968", "rs2519015", "rs2905042", "rs376645387", "rs2977608", "rs2977607", "rs2977606", "rs3095826", "rs2519006",\
    #                         "rs2905039", "rs28830877", "rs747980890", "rs778519663", "rs2905037", "rs2905035", "rs757138879", "rs12124819", "rs2980319",\
    #                         "rs112119688", "rs2977615", "rs2977613", "rs2977612", "rs2980295", "rs2905062", "rs2980300", "rs59319764", "rs2905061",\
    #                         "rs28753393", "rs2905058", "rs2905057", "rs2905056", "rs2905055", "rs187283089", "rs2905054", "rs2905053", "rs2980306" \
    #                         "rs10158938", "rs76089329", "rs72558500", "rs1044922", "rs2905036", "rs2519067", "rs72631879", "rs2909612", "rs12076540",\
    #                         "rs4951864", "rs4245756", "rs6681049", "rs4951931", "rs7516866", "rs9700144", "rs12565310", "rs9725068", "rs11240778",\
    #                         "rs4951933", "rs11240779", "rs6594027", "rs11240780", "rs797008946", "rs7541694", "rs7545373", "rs4246500", "rs4246499",\
    #                         "rs199713952", "rs61768235", "rs201418719", "rs4268301", "rs200699159", "rs376363419", "rs201851918", "rs4571914", "rs72888851",\
    #                         "rs150353320", "rs61768236", "rs61768237", "rs75814620", "rs199660488", "rs4970339", "rs4970338", "rs75468428", "rs112557229",\
    #                         "rs4970387", "rs200893860", "rs6605057", "rs149920886", "rs6677354", "rs35660652", "rs6671445", "rs13302980", "rs143626389",\
    #                         "rs10159337", "rs4475692", "rs13303179", "rs4970337", "rs369620271", "rs2879698", "rs4246498", "rs4970386", "rs377185360",\
    #                         "rs9778019", "rs4437820", "rs28444699", "rs6422669", "rs4970385", "rs9697642", "rs9697378", "rs9697380", "rs9697725",\
    #                         "rs34017275", "rs4500250", "rs4553118", "rs72631888", "rs35084403", "rs111389427", "rs7519340", "rs11516185", "rs7366404",\
    #                         "rs4970333", "rs35331099", "rs7416129", "rs13303222", "rs6664536", "rs6679046", "rs6657440", "rs72631889", "rs4970465",\
    #                         "rs4970463", "rs28436996", "rs7518702", "rs13303369", "rs4970461", "rs1806509", "rs13303019", "rs13303057", "rs6673914",\
    #                         "rs34628185", "rs7418179", "rs71509444", "rs71509445", "rs71509446", "rs112703963", "rs61464428", "rs57465118", "rs57924093",\
    #                         "rs60837925", "rs61338526", "rs773119229", "rs200244572", "rs57816555", "rs28521172", "rs2879816", "rs13302982", "rs13303291",\
    #                         "rs13303101", "rs6680268", "rs6693546", "rs3892970", "rs4040604", "rs28626846", "rs7410998", "rs7417972", "rs7417994",\
    #                         "rs2340589", "rs2340588", "rs2340587", "rs29938", "rs552943457", "rs3079263", "rs2203258", "rs72352263", "rs142286746",\
    #                         "rs2340590", "rs4970464", "rs6689107", "rs2905060", "rs10157303", "rs367950410", "rs140363434", "rs4970388", "rs2977605"], 1500)# just here for testing.  
    
    
    # snp_df = [  {'snpID': 'rs1799971', 'allele': 'G', 'sequence': 'TCCTGGGTCAACTTGTCCCACTTAGATGGCGACCTGTCCGACCCATGCGGTCCGAACCGCA', 'position': 30}, 
    #             {'snpID': 'rs12184297', 'allele': 'T', 'sequence': 'CTTTAAACCTCAACACATTATCAAGCATAATACTGTATATAATAAGTACTCAATACTGAAT', 'position': 30}, 
    #             {'snpID': 'rs116801199', 'allele': 'G', 'sequence': 'TAAAAAATGAATCTAATAATGAGGAAACATGAGAAAAAACCAAACTGAGGGATATTCTACA', 'position': 30}, 
    #             {'snpID': 'rs116801199', 'allele': 'T', 'sequence': 'TAAAAAATGAATCTAATAATGAGGAAACATTAGAAAAAACCAAACTGAGGGATATTCTACA', 'position': 30}, 
    #             {'snpID': 'rs12565286', 'allele': 'G', 'sequence': 'GGAAGCATCCTTCACTATCTTCTACCAAGGGCTTCCTCCTTTGGTGCTTCAAAATTTTTTA', 'position': 30}, 
    #             {'snpID': 'rs12565286', 'allele': 'C', 'sequence': 'GGAAGCATCCTTCACTATCTTCTACCAAGGCCTTCCTCCTTTGGTGCTTCAAAATTTTTTA', 'position': 30}, 
    #             {'snpID': 'rs2977670', 'allele': 'G', 'sequence': 'AACCTTGGAGGACCTATTGCTTAAGGTGTGGGCCAAAGAAAGTAAGTTAGGGCAAGAGACT', 'position': 30}, 
    #             {'snpID': 'rs2977670', 'allele': 'A', 'sequence': 'AACCTTGGAGGACCTATTGCTTAAGGTGTGAGCCAAAGAAAGTAAGTTAGGGCAAGAGACT', 'position': 30}, 
    #             {'snpID': 'rs2977670', 'allele': 'C', 'sequence': 'AACCTTGGAGGACCTATTGCTTAAGGTGTGCGCCAAAGAAAGTAAGTTAGGGCAAGAGACT', 'position': 30}, 
    #             {'snpID': 'rs2977670', 'allele': 'T', 'sequence': 'AACCTTGGAGGACCTATTGCTTAAGGTGTGTGCCAAAGAAAGTAAGTTAGGGCAAGAGACT', 'position': 30}, 
    #             {'snpID': 'rs28454925', 'allele': 'C', 'sequence': 'GGATTCGAATGGAAAGACATGGAATGGACTCGATTGGAATGGGTTGGGATGGAATGATCTA', 'position': 30}, 
    #             {'snpID': 'rs28454925', 'allele': 'G', 'sequence': 'GGATTCGAATGGAAAGACATGGAATGGACTGGATTGGAATGGGTTGGGATGGAATGATCTA', 'position': 30}, 
    #             {'snpID': 'rs28454925', 'allele': 'T', 'sequence': 'GGATTCGAATGGAAAGACATGGAATGGACTTGATTGGAATGGGTTGGGATGGAATGATCTA', 'position': 30}, 
    #             {'snpID': 'rs116587930', 'allele': 'G', 'sequence': 'ATTTTCAACTTTTGTAAATCTCTGTTTTAGGTGGGCTTCTTACGTACAACTTGGAGTTGGG', 'position': 30}, 
    #             {'snpID': 'rs116587930', 'allele': 'A', 'sequence': 'ATTTTCAACTTTTGTAAATCTCTGTTTTAGATGGGCTTCTTACGTACAACTTGGAGTTGGG', 'position': 30}, 
    #             {'snpID': 'rs116720794', 'allele': 'C', 'sequence': 'CTCTAACAGGCATTTCAGAGTGAGGTGGGACGTTCTAGGGCACCTGTTTTGCAGATGCCCT', 'position': 30}, 
    #             {'snpID': 'rs116720794', 'allele': 'T', 'sequence': 'CTCTAACAGGCATTTCAGAGTGAGGTGGGATGTTCTAGGGCACCTGTTTTGCAGATGCCCT', 'position': 30}, 
    #             {'snpID': 'rs4951859', 'allele': 'C', 'sequence': 'TTTGCAGATGCCCTCAGGGTGGGGGAAGGGCAGCTTCCAGCCTTCCCAGTTCCAGCACTCT', 'position': 30}, 
    #             {'snpID': 'rs4951859', 'allele': 'G', 'sequence': 'TTTGCAGATGCCCTCAGGGTGGGGGAAGGGGAGCTTCCAGCCTTCCCAGTTCCAGCACTCT', 'position': 30}, 
    #             {'snpID': 'rs4951859', 'allele': 'T', 'sequence': 'TTTGCAGATGCCCTCAGGGTGGGGGAAGGGTAGCTTCCAGCCTTCCCAGTTCCAGCACTCT', 'position': 30}, 
    #             {'snpID': 'rs148120343', 'allele': 'T', 'sequence': 'TCCCTTCCTTCCAATTCTCCTTCCAGCCTTTCTTGATTTCCAGAATGAGAAATCATTAAGT', 'position': 30}, 
    #             {'snpID': 'rs148120343', 'allele': 'C', 'sequence': 'TCCCTTCCTTCCAATTCTCCTTCCAGCCTTCCTTGATTTCCAGAATGAGAAATCATTAAGT', 'position': 30}, 
    #             {'snpID': 'rs529266287', 'allele': 'TA', 'sequence': 'AGAGAAAGTCCAGTCAATTTTATATAAGTTTAAAAAAAAGATGTGAAACCTATTTTCAGAAT', 'position': 30}, 
    #             {'snpID': 'rs79643588', 'allele': 'A', 'sequence': 'CATGGTCTAGGGAAGGAGAATGAAACATCAAAAATAACTGCAATTCCCCACAGTACGTGTC', 'position': 30}, 
    #             {'snpID': 'rs17396518', 'allele': 'T', 'sequence': 'GCAAAAATTTACAGAGAAGGAAATAGAGCTTCTCCCAAAATGTTAATAAAATTCTTAAAGG', 'position': 30}, 
    #             {'snpID': 'rs17396518', 'allele': 'A', 'sequence': 'GCAAAAATTTACAGAGAAGGAAATAGAGCTACTCCCAAAATGTTAATAAAATTCTTAAAGG', 'position': 30}, 
    #             {'snpID': 'rs17396518', 'allele': 'C', 'sequence': 'GCAAAAATTTACAGAGAAGGAAATAGAGCTCCTCCCAAAATGTTAATAAAATTCTTAAAGG', 'position': 30}, 
    #             {'snpID': 'rs983166', 'allele': 'C', 'sequence': 'GGCTAATATAATACTTATGGAACACTACCACTGTGCCAGATACTACTGATAAATGTTATAT', 'position': 30}, 
    #             {'snpID': 'rs983166', 'allele': 'G', 'sequence': 'GGCTAATATAATACTTATGGAACACTACCAGTGTGCCAGATACTACTGATAAATGTTATAT', 'position': 30}, 
    #             {'snpID': 'rs983166', 'allele': 'T', 'sequence': 'GGCTAATATAATACTTATGGAACACTACCATTGTGCCAGATACTACTGATAAATGTTATAT', 'position': 30}, 
    #             {'snpID': 'rs28842593', 'allele': 'A', 'sequence': 'GATATGTTTTGCATATGATACTCCATTGTAAAGCAGCAACAGCTAGAACTAAGCTGTTGTA', 'position': 30}, 
    #             {'snpID': 'rs28842593', 'allele': 'C', 'sequence': 'GATATGTTTTGCATATGATACTCCATTGTACAGCAGCAACAGCTAGAACTAAGCTGTTGTA', 'position': 30}, 
    #             {'snpID': 'rs7014597', 'allele': 'A', 'sequence': 'TGTCAAGGCCACCCTGGGCTTGAAGGGACCAGCCATGCCTCCAAGCCTTGCCCAGAGAGGG', 'position': 30}, 
    #             {'snpID': 'rs7014597', 'allele': 'C', 'sequence': 'TGTCAAGGCCACCCTGGGCTTGAAGGGACCCGCCATGCCTCCAAGCCTTGCCCAGAGAGGG', 'position': 30}] 
    # snp_df = Fetch_SNP_Data(['rs1044922', 'rs6664536', 'rs367950410'], 30)
    # print(snp_df)
    snp_end = datetime.now()

    primers = generate_allele_specific_primers(snp_df, 18, 24)
    primer_close_end = datetime.now()

    filt = filter_all_list(primers)
    filter_end = datetime.now()

    best_primers, fights = multiplex_list(filt)
    print("fights")
    print(fights)
    multi_end = datetime.now()

    far = generate_matching_primers(best_primers[0], snp_df) 
    # print("far")
    # print(far)
    far_end = datetime.now()

    end = datetime.now()
    
    print(f"fetch_snp took : {snp_end - start}")
    print(f"generate_allele_specific_primer took : {primer_close_end - snp_end}")
    print(f"filter_primers took : {filter_end - primer_close_end}")
    print(f"multiplexing close took : {multi_end - filter_end}")
    print(f"generate far took : {far_end - multi_end}")
    print(f"total time : {end-start}")

if(__name__ == "__main__"):
    Main()

