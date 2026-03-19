# SNP multiplexer tool

* This program finds primer pair and SNP specific probe combinations for all of SNPs entered into the list
* It multiplexes all of the primer and probes to extract all that can sequence in the same run
* Any that fail will be exported with reason that they failed to be used again in with different parameters 

## What it does

This program takes in a list of SNP id's, generates probes at the SNP site, multiplexes the probes so they all can work together, the far primers are generated and then multiplexed against the probes to get the final result of a list of probe primer combinations for all the ones that work and a result file of all the ones that couldn't work. We (will) use Locked Nucleic Acids to force high specificity in the probes so they only match with the SNP they were designed for.



## How the app works:
step 1. 
    Fetch_SNP_Data takes in a list of rsIDs and the length to either side of the SNP location.
    it then makes an api call for the SNP and gets all alleles associated with it. It then makes 
    another call to get the DNA strand and puts the different alleles in the SNP location making a
    new entry for each allele.
step 2. 
    We pass the list to generate_allele_specific_probes which makes a list of probes for each allele returning a list of lists of Probe instances. They Probe class filters at creation so if it can't pass the user defined parameters it will fail before going any further.
step 3. 
    The list of list of good probes is passed to multiplex_close which will check the heterodimer of each close probe each other. It keeps trying combinations until all alleles' probes work, or have been exhausted and scrapped.
step 4. 
    Each allele then calls the primer generator, multiplexes the two far primers against all probes and primers that have already been created. If none of the far primers can multiplex then the primer generator is called again until a combination is found that works or the allele gets scrapped.
step 5.
    The results of the multiplexed alleles and scrapped alleles is printed as two PDFs.



[Primer step testing](https://docs.google.com/spreadsheets/d/10p1sGPLq-MyVW8i0swzkHR8p_JaZrSPe09oy_DVx0Qs/edit?usp=sharing)

[Allele Sequencing time O(N^3) complexity](https://www.desmos.com/calculator/5hqu4aeon9)

