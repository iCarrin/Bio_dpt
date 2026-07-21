# Bio_dpt
The main (first unzips a bunch of settings. That wasn't my work so I don't quite get it) takes the provided SNP ID's and makes a API call to a database. The information on the SNP comes back (what alleles have been seen, and where it's located in which genome) and a new API call is made to get the DNA around the SNP. The DNA strands are copied and each allele is spliced into a copy of that DNA strand.
From there the Probes are generated in multiplex_probes. The Probe has to have it's middle located on the SNP sight which limits what we can make. The Probe DNA has it's middle three bp turned into LNA's (locked nuecalic acids) to increase specificity. We try the largest and if that doesn't work we move to the next largest until we get down to the minimum side. 
For primers and Probes we test it for correct TM, Hairpining and Homodimering durring its creation. If it fails, it doesn't finish getting created. If no probe can be generated the SNP fails and no longer is processed.
As the SNP Allele combos have their probes generated they are tested against all previously generated probes to make sure they don't heterodime with anything else. 
With the list of passing probes, we move on to primer generation. The same steps are followed. For each primer, one side at a time, primers are tested at generation, then they have their heterodimer checked for all previous probes and primers. If a forward and reverse primer can't be found for a probe then the SNP Allele combo will be removed from the list. 
The final list of SNP Allele combos with their associated probes and primers are passed out in a PDF

TODO
add prime blasting
update LNA.
    currently we hard code the middle triplet as LNA but this might not be the best thing.
    Add Zack's code to find the largest TM delta instead of trusting the locked tripplet will always be the best.
Dummy testing 
    do with website because CLI is obsolete.


