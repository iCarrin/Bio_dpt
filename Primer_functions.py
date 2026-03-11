from Primer_Classes import *
from Bio.Seq import Seq
import logging

logger = logging.getLogger(__name__)


def generate_allele_specific_probes(snp_json: list[dict], min_len: int = 28, max_len: int = 32) -> list[list[Probe]]:
    all_probes = []
    def make_allele_probes_list(snp_id, allele, sequence, snp_pos):
        #add a check here for length so that make primers doesn't have to.
        flank_len = max_len - max_len//2#this gives the longer half
                            #this gives us the longer half so in case we need to drop a g form the 5' end it balances better
        forward = sequence[snp_pos - flank_len : snp_pos+(flank_len)+1]#this gets the largest segment.   

        reverse = str(Seq(sequence[snp_pos-flank_len : snp_pos+flank_len+1]).reverse_complement()) #creates a Biopython sequence, gets the reverse complement, and converts is back to a string

        return (make_probes(forward, min_len, snp_id, allele))\
                + (make_probes(reverse, min_len, snp_id, allele, "reverse"))


    for snp in snp_json:
        allele_probes = make_allele_probes_list(snp["snpID"], snp["allele"], snp["sequence"], snp["position"])
        if allele_probes:
            all_probes.append(allele_probes)
        else:
            print(f"snp: {snp["snpID"]} allele: {snp["allele"]} didn't make the cut")
    
    return all_probes

def make_probes(seq, min_len, snp_id, allele, direction="forward") -> list[Probe]: 
    
    probes = []


    if len(seq) >= min_len:   #len(seq) should just be max_len. It only won't be if the seq length is less than min_len (if the sequence is only 10 long then we'll trigger this)
        len_of_the_flank = (len(seq)-min_len)//2
        for length in range(len_of_the_flank+1):#possible bug if the forward mismatch is smaller than the minimum length
                                            #length is 0-max_len
            trimmed = seq[length: -length if length != 0 else None]#this is assuming that the seq given is already the maximum length. 
            # If given a crazy long string it will start at "length" and give the rest of the string
            #take this part out of the loop, so we can have one dictionary that says the SNP ID and ALLELE and Direction, 
            #and then a list in that dictionary of sequence and lengths. Storing the name over and over seems redundant ID

            try:                                                            #these need to be user controlled inputs
                probes.append(Probe(snp_id, allele, trimmed, direction, 70.0, 3.0, -3.0, -3.0))
            except FilterFail as e:
                # print(e)
                pass
            #These take one off the left and right to see if it gets us anywhere out of 34 it saved an extra 3
            try:                                                            #these need to be user controlled inputs
                probes.append(Probe(snp_id, allele, trimmed[:-1], direction, 70.0, 3.0, -3.0, -3.0))
            except FilterFail as e:
                # print(e)
                pass
            try:                                                            #these need to be user controlled inputs
                probes.append(Probe(snp_id, allele, trimmed[1:], direction, 70.0, 3.0, -3.0, -3.0))
            except FilterFail as e:
                # print(e)
                pass
               

    else:
        print(f"The length of your {direction} primer wasn't long enough. \nYou needed one at least {min_len} long and it ended up only being {len(seq)}")
        logger.warning(f"The length of your {direction} primer {snp_id} allele {allele} wasn't long enough. \nYou needed one at least {min_len} long and it ended up only being {len(seq)}")
    
    probes.sort(key=lambda x: x.rank)
    
    return probes

def make_primers(seq, min_len, max_len, snp_id, allele, direction="forward") -> list[dict]: 
    
    primers = []
    if len(seq) >= min_len:   #len(seq) should just be max_len. It only wont be if the seq length is less than max_len (if the sequence is only 10 long then we'll trigger this)
        for length in range(max_len-min_len):#possible bug if the forward mismatch is smaller than the minimum length
                                            #length is 0-max_len
            trimmed = seq[length:]#this is assuming that the seq given is already the maximum length. 
            # If given a crazy long string it will start at "length" and give the rest of the string
            #take this part out of the loop, so we can have one dictionary that says the SNP ID and ALLELE and Direction, 
            #and then a list in that dictionary of sequence and lengths. Storing the name over and over seems redundant ID
        
            try:                                                            #these need to be user controlled inputs
                primers.append(Primer(snp_id, allele, trimmed, direction, 60.0, 3.0, -3.0, -3.0))
            except FilterFail as e:
                # print(e)
                pass
                # logging.error(f"{snp_id} allele: {allele} had no primers that passed the filtering")

    else:
        print(f"The length of your {direction} primer wasn't long enough. \nYou needed one at least {min_len} long and it ended up only being {len(seq)}")
        logger.warning(f"The length of your {direction} primer {snp_id} allele {allele} wasn't long enough. \nYou needed one at least {min_len} long and it ended up only being {len(seq)}")
    return primers




def generate_matching_primers(primer_king, snp_json, primer_start = 0, direction = 'positive', min_len = 18, max_len = 24, min_dist: int = 50, max_dist: int = 250): 
    """
        Generate matching primers for top  allele-specific primers.
        TODO: Optimize primer pairing.
        - Use primer3-py's designPrimers for more efficient pairing.
        - Add checks for primer pair compatibility (e.g., Tm difference < 5Â°C).
    """
    if len(snp_json[0]['sequence']) < (min_dist+17)/2:
        raise Exception("your sequence is so short it won't allow for even 1 primer to be maid. " \
        "Lower your min distance to have the API call for a longer string")
    snp_dict = {}
    
    for snp in snp_json:
        if  snp['snpID'] == primer_king.snpID: # we only search for a matching SNP because we aren't even using the allele section    
            snp_dict = snp
            break   
    else:
        raise Exception(f"there is no matching entry in the Json list for {primer_king['snpID']}")
   
    
    middle = snp_dict['position']
    whole_sequence = snp_dict['sequence']
    

    #get the far sequence and reverse complement is if necessary
    if direction == "positive":
        seq_pre_rev_comp = whole_sequence[middle+min_dist:middle+max_dist]#everything from snp (middle) plus start dist cutting off at max if necessary (end is exclusive so +1)
        far_sequence = str(Seq(seq_pre_rev_comp).reverse_complement()) #change to sequence object, reverse complement it, and change it back to string
        #                                         _______  
        #                      _________          \______\_
        #  What we're given _/__________|           \_______\ what we make 
        #_________________/_____________|=====================("===" means reverse complement DNA)
        
    elif direction == "negative":
        far_sequence = whole_sequence[middle-max_dist:middle-min_dist-1]
        #___________________________
        #     ___/______/           ====================== ("===" means reverse complement DNA)
        #    /_______/ What makin'  |_____________/ what we have (reverse close)
        #                           |___________/     
    else:
        raise Exception("no direction given. How did we get here?")
    


    #get some far primers for primer_king, but only the best (use strict mode). if the filter fails the while loop should try again farther down the line
    #if the ones that pass don't pass the heterodimer then primer_start should be updated and the whole thing tried again.
    while(True):
        # the reverse complement flips the sequence each time so we can iterate over each one the same way. Researchers expect DNA to come in the forward
        # format even if it was reversed in real life, so no need to flip it back. The note that it should be reversed is enough.
        # each one will walk back from the right side 
        trial_snip = far_sequence[-(max_len+primer_start) : -primer_start or None] # this takes a chunk to feed into a primer generator
        
      
        try:#filter strict mode will throw an error so we use try except
            far_primers = make_primers(trial_snip, min_len, max_len, primer_king.snpID, primer_king.allele, direction)
            # filt_far, _ = filter_one_list(far_primers, temp, 2, strict_mode=True)\
            if far_primers:
                break
            else:
                primer_start += 6
        except ValueError as e:
            print(e)
            logger.warning(e)
            primer_start += 6 #extensive testing shows that 6 is the fastest step to run (See readme for link)
            if primer_start > max_dist-min_len:
                logger.critical(f'{primer_king.snpID} allele {primer_king.allele} had no useable far primers')
                raise Exception("you've tried every possible primer. What have you done??")
        
    
    return far_primers, primer_start

