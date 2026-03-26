from Multiplexer_class import Multiplexer
from fetch_snp_data import get_data

"""
How the app works:
step 1. Fetch_SNP_Data takes in a list of rsIDs and the length to either side of the SNP location.
    it then makes an api call for the SNP and gets all alleles associated with it. It then makes 
    another call to get the DNA strand and puts the different alleles in the SNP location and it 
    count(allele) number of strands

"""

def Main():
    snp_df=get_data()
    if len(snp_df)==0:
        raise Exception(f"{len(snp_df)=}")
    d1={}
    print("enter values for the things below press 'enter' for default:")
    for i,ty in [("mv_conc (default 50.0): ",float), ("dv_conc (default 3): ",int), ("dntp_conc (default 0.8): ",float), ("dna_conc (default 200): ",int), ("dmso_conc (default 0.0): ",float), ("dmso_fact (default 0.0): ",float), ("formamide_conc (default 0.0): ",float), ("annealing_temp_c (default -10.0): ",float), ("temp_c (default 37.0): ",float), ("tm_method (default 'santalucia'): ",str),("salt_corrections_method (default 'owczarzy'): ",str)]:
        if (t:=input(i)):
            d1[i.split()[0]]=ty(t)
    
    Multiplexer(snp_df=snp_df,**d1).main()
    
if(__name__ == "__main__"):
    Main()
