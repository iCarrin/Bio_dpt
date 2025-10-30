# from Multiplex import *
from Output import *
from Primer_functions import *
import primer3
import re # run 'pip install regex' if not already installed
import time # to handle rate limiting
import requests
from typing import List # for type hinting

# Ensembl REST API base URL
ENSEMBL_REST = "https://rest.ensembl.org"

def Fetch_SNP_Data(rsids: List[str], flank_length: int = 800) -> list[dict]:
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
            ancestral = mapping.get("ancestral_allele")

            #Isaiah change: I dropped any ancestral Alleles that we knew were normal
            #In the future we should pull the list, split it, and have the user pick which ones they want.
            alleles = allele_str.split("/") if allele_str else []
            if ancestral in ['A', 'C', 'G', 'T']:
                alleles.remove(ancestral)


            # Ensure there are at least two alleles to work with
            if len(alleles) < 1:
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
            for allele in alleles:
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
    
    # real fetch snp

    # snp_df = Fetch_SNP_Data(["rs1799971", "rs12184297", "rs116801199", "rs12565286", "rs2977670", "rs28454925"], 30)# just here for testing.  , "rs599839"


    snp_df = [{'snpID': 'rs1799971', 'allele': 'G', 'sequence': 'TCCTGGGTCAACTTGTCCCACTTAGATGGCGACCTGTCCGACCCATGCGGTCCGAACCGCA', 'position': 30}, 
                {'snpID': 'rs12184297', 'allele': 'T', 'sequence': 'CTTTAAACCTCAACACATTATCAAGCATAATACTGTATATAATAAGTACTCAATACTGAAT', 'position': 30}, 
                {'snpID': 'rs116801199', 'allele': 'G', 'sequence': 'TAAAAAATGAATCTAATAATGAGGAAACATGAGAAAAAACCAAACTGAGGGATATTCTACA', 'position': 30}, 
                {'snpID': 'rs116801199', 'allele': 'T', 'sequence': 'TAAAAAATGAATCTAATAATGAGGAAACATTAGAAAAAACCAAACTGAGGGATATTCTACA', 'position': 30}, 
                {'snpID': 'rs12565286', 'allele': 'G', 'sequence': 'GGAAGCATCCTTCACTATCTTCTACCAAGGGCTTCCTCCTTTGGTGCTTCAAAATTTTTTA', 'position': 30}, 
                {'snpID': 'rs12565286', 'allele': 'C', 'sequence': 'GGAAGCATCCTTCACTATCTTCTACCAAGGCCTTCCTCCTTTGGTGCTTCAAAATTTTTTA', 'position': 30}, 
                {'snpID': 'rs2977670', 'allele': 'G', 'sequence': 'AACCTTGGAGGACCTATTGCTTAAGGTGTGGGCCAAAGAAAGTAAGTTAGGGCAAGAGACT', 'position': 30}, 
                {'snpID': 'rs2977670', 'allele': 'A', 'sequence': 'AACCTTGGAGGACCTATTGCTTAAGGTGTGAGCCAAAGAAAGTAAGTTAGGGCAAGAGACT', 'position': 30}, 
                {'snpID': 'rs2977670', 'allele': 'C', 'sequence': 'AACCTTGGAGGACCTATTGCTTAAGGTGTGCGCCAAAGAAAGTAAGTTAGGGCAAGAGACT', 'position': 30}, 
                {'snpID': 'rs2977670', 'allele': 'T', 'sequence': 'AACCTTGGAGGACCTATTGCTTAAGGTGTGTGCCAAAGAAAGTAAGTTAGGGCAAGAGACT', 'position': 30}, 
                {'snpID': 'rs28454925', 'allele': 'C', 'sequence': 'GGATTCGAATGGAAAGACATGGAATGGACTCGATTGGAATGGGTTGGGATGGAATGATCTA', 'position': 30}, 
                {'snpID': 'rs28454925', 'allele': 'G', 'sequence': 'GGATTCGAATGGAAAGACATGGAATGGACTGGATTGGAATGGGTTGGGATGGAATGATCTA', 'position': 30}, 
                {'snpID': 'rs28454925', 'allele': 'T', 'sequence': 'GGATTCGAATGGAAAGACATGGAATGGACTTGATTGGAATGGGTTGGGATGGAATGATCTA', 'position': 30}]
    # print(snp_df)

    primers = generate_allele_specific_primers(snp_df, 24, 30)
    # for prime_list in primers:
    #     for primer in prime_list:
    #         print(primer)
    #     print()
    # print(f"number of snps {len(primers)}")
    # print(f"number of directions {len(primers[0])}")
    # print(f"number of primers in a direction {len(primers[0][0])}")
    # print(len(snp_df))
    low = 0
    high = 0
    dimer = 0
    hairpin = 0
    for allele in primers:
        
        allele_list, fail_ints = filter_one_list_soft(allele, diff = 5.0)
        low += fail_ints[0]
        high += fail_ints[1]
        dimer += fail_ints[2]
        hairpin += fail_ints[3]
        # print(allele_list)
    print(low)
    print(high)
    print(dimer)
    print(hairpin)

if(__name__ == "__main__"):
    Main()

