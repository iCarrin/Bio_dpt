# Bio_dpt
The main (first unzips a bunch of settings. That wasn't my work so I don't quite get it) takes the provided SNP ID's and makes a API call to a database. The information on the SNP comes back (what alleles have been seen, and where it's located in which genome) and a new API call is made to get the DNA around the SNP. The DNA strands are copied and each allele is spliced into a copy of that DNA strand.
From there the Probes are generated in multiplex_probes. The Probe has to have it's middle located on the SNP sight which limits what we can make. The Probe DNA has it's middle three bp turned into LNA's (locked nuecalic acids) to increase specificity. We try the largest and if that doesn't work we move to the next largest until we get down to the minimum side. 
For primers and Probes we test it for correct TM, Hairpining and Homodimering durring its creation. If it fails, it doesn't finish getting created. If no probe can be generated the SNP fails and no longer is processed.
As the SNP Allele combos have their probes generated they are tested against all previously generated probes to make sure they don't heterodime with anything else. 
With the list of passing probes, we move on to primer generation. The same steps are followed. For each primer, one side at a time, primers are tested at generation, then they have their heterodimer checked for all previous probes and primers. If a forward and reverse primer can't be found for a probe then the SNP Allele combo will be removed from the list. 
The final list of SNP Allele combos with their associated probes and primers are passed out in a PDF

TODO
add prime blasting * this is low priority
update LNA.
    currently we hard code the middle triplet as LNA but this might not be the best thing.
    Add Bio Info Zack's code to find the largest TM delta instead of trusting the locked tripplet will always be the best. * top priority
Dummy testing 
    do with website because CLI is obsolete.
add in shifting DNA string


News
The DNA can slide around the probe as long as no overlap
total length 
max 150
optimal 100
min 75
forward primer minimum gap distance is 5
reverse probe is hanging down below so it doesn't really have to worry about even over lapping the probe as long as it applifies the dna for future probes and it doesn't shoot the total length down to < 75
        This will not apply for right now however because probes for the same SNP might be flipped we'll just insure there is no overlap
reverse stays clear of revers probes and visa versa
60 - 61 - 62 diff of a 1 degree the diff range (60-61 or 61-62) applies to individual targets. The vat can still fall within the range of 60-62
probes at the same melting temp
the primer blasting should just print out the results instead of removing it


Make sure the probes of the same SNP all fall in the same range 

(link to the study on LNA's)[https://pubs.acs.org/bichaw/article-pdf/50/43/9352/7881128/bi200904e.pdf]

