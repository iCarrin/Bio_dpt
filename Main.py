import primer3,logging
from Primer_functions import generate_allele_specific_probes
from Multiplex import multiplex_close, multiplex_far
from pdfoutput import create_output_json
from datetime import datetime 
from fetch_snp_data import get_data



"""
How the app works:
step 1. Fetch_SNP_Data takes in a list of rsIDs and the length to either side of the SNP location.
    it then makes an api call for the SNP and gets all alleles associated with it. It then makes 
    another call to get the DNA strand and puts the different alleles in the SNP location and it 
    count(allele) number of strands

"""

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
    
    logging.basicConfig(
        filename="primer_info.log",
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )


    snp_df=get_data()
    # primer_3.set_thermo_args(
    #     mv_conc=50.0, 
    #     dv_conc=3, 
    #     dntp_conc=0.8, 
    #     dna_conc=200, 
    #     dmso_conc=0.0, 
    #     dmso_fact=0.0, 
    #     formamide_conc=0.0, 
    #     annealing_temp_c=-10.0, 
    #     temp_c=37.0, 
    #     tm_method='santalucia',
    #     salt_corrections_method='owczarzy')
    # print(thermo_settings.calc_heterodimer(snp_df[0]["sequence"],snp_df[-1]["sequence"]))
    # print(primer3.calc_heterodimer(snp_df[0]["sequence"],snp_df[-1]["sequence"]))
    # exit()


    logger = logging.getLogger(__name__)
    snp_end = datetime.now()
    primers,bad = generate_allele_specific_probes(snp_df, 28, 32)
    primer_close_end = datetime.now()
    best_primers, fights = multiplex_close(primers)
    multi_end = datetime.now()
    poasitive, negative, center, center_reject = multiplex_far(best_primers, snp_df)
    far_end = datetime.now()
    end = datetime.now()
    logger.info(f"allele's run {len(snp_df)}")
    logger.info(f"fetch_snp took : {snp_end - start}")
    logger.info(f"generate_allele_specific_primer took : {primer_close_end - snp_end}")
    # logger.info(f"filter_primers took : {filter_end - primer_close_end}")
    logger.info(f"multiplexing close took : {multi_end - snp_end}")
    logger.info(f"generate far took : {far_end - multi_end}")
    logger.info(f"total time : {end-start}")

    create_output_json(best_primers,"final_primer.pdf",(1000,1000))
    create_output_json(bad,"final_bad.pdf",(12000,1000))
    create_output_json(poasitive,"final_positive.pdf",(1000,1000))
    create_output_json(negative,"final_negitive.pdf",(1000,1000))
    create_output_json(center,"final_center.pdf",(1000,1000))
    create_output_json(center_reject,"final_center_reject.pdf",(900,1000))
    
    
if(__name__ == "__main__"):
    Main()
