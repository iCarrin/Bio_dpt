
import re 
from Bio.Seq import Seq
import primer3
import logging

logger = logging.getLogger(__name__)
def introduce_mismatch(primer_sequence: str, fall_back=False) -> str:
    """
    Introduces a base mismatch at the antepenultimate position (3rd from last).
    """
    # Ensure valid string input
    if not primer_sequence or not isinstance(primer_sequence, str):
        print("Warning: Invalid primer input.")
        return primer_sequence

    primer_sequence = primer_sequence.upper().strip()

    # Must only contain A, C, G, T
    if not re.match("^[ACGT]+$", primer_sequence):
        print(f"Warning: Invalid characters in primer: {primer_sequence}")
        return primer_sequence

    # Must be long enough to have a 3rd-to-last base
    if len(primer_sequence) < 3:
        print(f"Warning: Primer too short for mismatch: {primer_sequence}")
        return primer_sequence


    if fall_back:
        mismatch_rules = {#medium rule set
            #technically A->T is a strong mismatch but it was too much work to implement two A options 
            #in the strong rule set and none in the medium, so it sits in the medium rule set
            "A": "T", "G": "C",
            "C": "A", "T": "A"
        }
    else:# Simple mismatch rules (purine↔pyrimidine, purine↔pyrimidine)
        mismatch_rules = { #strong rule set
            "A": "C", "G": "T",
            "C": "A", "T": "G"
        }


    pos = len(primer_sequence) - 3  # Antepenultimate index
    base = primer_sequence[pos]
    mismatch = mismatch_rules.get(base)

    if mismatch is None:
        print(f"Warning: No mismatch rule for base '{base}'")
        return primer_sequence

    # Replace the base with its mismatch
    return primer_sequence[:pos] + mismatch + primer_sequence[pos + 1:]


def calc_gc_content(sequence: str):
    gc_total = 0
    for nucleotide in sequence:
        if nucleotide == 'G' or nucleotide == 'C':
            gc_total += 1
    return gc_total/len(sequence)


def rank_primers(primers: list[dict], target_tm = 60, target_gc = 50) -> list[dict]:
    """
        Rank primers based on Tm proximity to 62.5Â°C and GC content.
        TODO: Refine ranking criteria.
        - Consider weighting Tm vs. GC scores.
        - Add user-configurable ranking metrics.
        """
    for primer in primers:
        primer["tm_score"] = abs(primer["tm"] - target_tm)
        primer["gc_score"] = abs(primer["gc_content"] - target_gc)
        primer["score"] = primer["tm_score"] + primer["gc_score"] 

    primers.sort(key=lambda x: x["score"])

    return primers



def filter_little(filter_name: str, old_list: list[dict], filter_function, strict):
    # fail_count is only 0 or 1, but we'll add it to a total to see what's bugging out
    fail_count = 0
    new_list = list(filter(filter_function, old_list))
    #if after the filter we don't have anything left
    if not new_list and not strict:
        new_list = old_list
        #the fail count goes up
        fail_count = 1
        log_info = []
        if  filter_name =="tm Max":
            for primer in old_list:
                log_info.append(f"{primer["tm"]} ")
            log_info.sort()
        elif filter_name == "tm Min":
            for primer in old_list:
                log_info.append(f"{primer["tm"]} ")
            log_info.sort(reverse=True)
        elif filter_name == "homodimer":
            for primer in old_list:
                log_info.append(f"{primer["homodimer"]} ")
            log_info.sort(reverse=True)
        else:
            for primer in old_list:
                log_info.append(f"{primer["hairpin"]} ")
            log_info.sort(reverse=True)
        #log the fail
        logger.warning(f"{old_list[0]["snpID"]} allele {old_list[0]["allele"]} failed {filter_name}. Results: {log_info}")

    elif not new_list and strict:
        raise ValueError("nothing passed in strict mode")
        
    #return the list that hopefully was able to filter and weather or not it failed
    return(new_list, fail_count)


def filter_one_list(allele_list: list[dict],
                         desired_tm: float = 60.0,
                         diff: float = 3.0,
                         homodimer_goal: float = -3.0,
                         hairpin_goal: float = -3.0,
                         strict_mode = False) -> (list[dict], list[int]):
    
    """
    Soft filter a single candidate list such as the stage1_filter behavior
    Order is homodimer, hairpin, tm > lower, tm < upper
    
    Applies the 4 checks and then returns filtered list of primer string
    for that one candidate list

    The checks are homodimer, hairpin, tm, 

    Used in the R stage1_filter:
        homodimer < homodimer_goal
        hairpin   < hairpin_goal
        Tm        < desired_tm + diff      # "above upper" trim
        Tm        > desired_tm - diff      # "below lower" trim
    """
    #make sure we're getting an actuall list of dictionaries
    if not allele_list:
        raise Exception("There was not list of dictionaries passed in")
    
    # tm > min
    allele_pltm, ltm_fail_count = filter_little("tm Min", allele_list, lambda x : x["tm"] >= (desired_tm - diff), strict_mode)
    
    # tm < max
    allele_phtm, htm_fail_count = filter_little("tm Max", allele_pltm, lambda x : x["tm"] <= (desired_tm + diff), strict_mode)
    # min < homodimer < max
    allele_phomo, homo_fail_count = filter_little("homodimer", allele_phtm, lambda x : ( x["homomdimer"].dg > homodimer_goal* 1000 or x["homomdimer"].tm < 40), strict_mode)
    #min < hairpin < max
    allele_phair, hair_fail_count = filter_little("hairpin", allele_phomo, lambda x : ( x["hairpin"].dg > hairpin_goal* 1000 or x["homomdimer"].tm < 40), strict_mode)
    
    # total_fails = [f"low temp fails: {ltm_fail_count}", f"high temp fails: {htm_fail_count}", f"homodimer fails: {homo_fail_count}", f"hairpin fails: {hair_fail_count}"]
  
    total_fails_ints = [ltm_fail_count, htm_fail_count, homo_fail_count, hair_fail_count]
    result = rank_primers(allele_phair)

    return (result, total_fails_ints)

def filter_all_list(all_snp_primers, desired_tm: float = 60.0, diff: float = 3.0, homodimer_goal: float = -3.0, hairpin_goal: float = -3.0, strict_mode = False) -> (list[dict], list[int]):

    filt_all_primers = []
    tm_min_fails = 0
    tm_max_fails = 0
    homodimer_fails = 0
    hairpin_fails = 0
    for allele_primer in all_snp_primers:
        good_primers, fails = filter_one_list(allele_primer, desired_tm, diff, homodimer_goal, hairpin_goal, strict_mode)
        filt_all_primers.append(good_primers)
        tm_min_fails += fails[0]
        tm_max_fails += fails[1]
        homodimer_fails += fails[2]
        hairpin_fails += fails[3]
    logger.info(f"lower TM fails: {tm_min_fails}, all tm too high: {tm_max_fails}, low homodimer with > 40C tm: {homodimer_fails}, low haripin with > 40C tm: {hairpin_fails}")
    return filt_all_primers
def generate_allele_specific_primers(snps_list: list[dict], min_len: int = 18, max_len: int = 24, fallback=False) -> list[list[dict]]:
    # make primers makes a dictionary for every length of one direction of an allele for a SNP.
    # those dictionaries are stored in a list, so a list for forward and a list for backward
    # those lists are stored in another list, one for each SNP. 
    # then those are stored in one big list, a list of all SNP given this session.
    # list_of_every_SNP[list_of_forward_and_reverse_directions[list_of_dictionaries[dictionary_of_particular_length]]]
    """
        Generate allele-specific primers (forward/reverse) ending at the SNP.
        TODO: Optimize for large SNP sets.
        - Use parallel processing (e.g., multiprocessing) for many SNPs.
        - Add validation for sequence length and SNP position.
        """
    all_primers = []
    min_len -= 2 #don't know why but 2 and 1 have to be removed from the inputs to get the desired lengths
    max_len -= 1
    def make_allele_list():
            forward = sequence[snp_pos - max_len :snp_pos+1]#this gets the largest segment.   
            forward_mismatch = introduce_mismatch(forward, fallback)
        
            reverse = str(Seq(sequence[snp_pos:snp_pos+max_len+1]).reverse_complement()) #creates a Biopython sequence, gets the reverse complement, and converts is back to a string
            reverse_mismatch = introduce_mismatch(reverse, fallback)

            return (make_primers(forward_mismatch, min_len, max_len, snp_id, allele))\
                    + (make_primers(reverse_mismatch, min_len, max_len, snp_id, allele, "reverse"))# this make one list of 
    for snp_dict in snps_list:
        #why pass this in seperately when it's already in the dict?
        snp_id = snp_dict["snpID"]
        allele = snp_dict["allele"]
        sequence = snp_dict["sequence"]
        snp_pos = snp_dict["position"]
        allele_primers = make_allele_list()
    
        all_primers.append(allele_primers) #this adds this list to the larger list
    return all_primers # this will return a list of lists of dictionaries. Each allele is a list. 

   

def make_primers(seq, min_len, max_len, snp_id, allele, direction="forward") -> list[dict]: 
    seq_length = len(seq)
    
    primers = []
    if seq_length >= min_len:
        for length in range(max_len-min_len):#possible bug if the forward mismatch is smaller than the minimum length

            trimmed = seq[length:]#this is assuming that the seq given is already the maximum length. 
            # If given a crazy long string it will start at "length" and give the rest of the string
            #take this part out of the loop, so we can have one dictionary that says the SNP ID and ALLELE and Direction, 
            #and then a list in that dictionary of sequence and lengths. Storing the name over and over seems redundant IDK
            primers.append({
                "snpID": snp_id,
                "allele": allele,
                "primer_sequence": seq[length:], #this is the trimmed length
                "direction": direction,
                "length": seq_length-length,
                "tm" : primer3.bindings.calc_tm(trimmed),
                "gc_content" : calc_gc_content(trimmed),
                "hairpin" : primer3.bindings.calc_hairpin(trimmed),
                "homomdimer" : primer3.bindings.calc_homodimer(trimmed)
            })
    else:
        print(f"The length of your {direction} primer wasn't long enough. \nYou needed one at least {min_len} long and it ended up only being {seq_length}")
        logger.warning(f"The length of your {direction} primer {snp_id} allele {allele} wasn't long enough. \nYou needed one at least {min_len} long and it ended up only being {seq_length}")
    return primers




def generate_matching_primers(primer_king, snp_json, primer_start = 0, min_len = 18, max_len = 24, min_dist: int = 800, max_dist: int = 1200): 
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
    found = False
    for snp in snp_json:
        if  snp['snpID'] == primer_king['snpID']: # we only search for a matching SNP because we aren't even using the allele section    
            snp_dict = snp
            found = True
            break
        
    if found == False:
        raise Exception(f"there is no matching entry in the Json list for {primer_king['snpID']}")
   
    
    middle = snp_dict['position']
    whole_sequence = snp_dict['sequence']
    direction = primer_king['direction']
    far_sequence = ""
    temp = primer_king['tm']

    #get the far sequence and reverse complement is if necessary
    if direction == "forward":
        seq_pre_rev_comp = whole_sequence[middle+min_dist:middle+max_dist]#everything from snp (middle) plus start dist cutting off at max if necessary (end is exclusive so +1)
        far_sequence = str(Seq(seq_pre_rev_comp).reverse_complement()) #change to sequence object, reverse complement it, and change it back to string
        #                                         _______  
        #                      _________          \______\_
        #  What we're given _/__________|           \_______\ what we make 
        #_________________/_____________|=====================("===" means reverse complement DNA)
        
    elif direction == "reverse":
        far_sequence = whole_sequence[middle-max_dist:middle-min_dist-1]
        #___________________________
        #     ___/______/           ====================== ("===" means reverse complement DNA)
        #    /_______/ What makin'  |_____________/ what we have (reverse close)
        #                           |___________/     
    else:
        raise Exception("no direction given. How did we get here?")
    
    passes = False
    start = primer_start

    #get some far primers for primer_king, but only the best (use strict mode). if the filter fails the while loop should try again farther down the line
    #if the ones that pass don't pass the heterodimer then primer_start should be updated and the whole thing tried again.
    while(not passes):
        # the reverse complement flips the sequence each time so we can iterate over each one the same way. Researchers expect DNA to come in the forward
        # format even if it was reversed in real life, so no need to flip it back. The note that it should be reversed is enough.
        # each one will walk back from the right side 

        trial_snip = far_sequence[-(max_len+start) : -start or None] # this takes a chunk to feed into a primer generator
        
        far_primers = make_primers(trial_snip, min_len, max_len, primer_king['snpID'], primer_king['allele'], direction)
      
        try:#filter strict mode will throw an error so we use try except
            filt_far, _ = filter_one_list(far_primers, temp, 2, strict_mode=True)
            passes = True
        except ValueError as e:
            print(e)
            logger.warning(e)
            start += 3 #extensive testing shows that 6 is the fastest step to run (See google sheet in the git hub)
            if start > max_dist-max_len:
                logger.critical(f'{primer_king['snpID']} allele {primer_king['allele']} had no useable far primers')
                raise Exception("you've tried every possible primer. What have you done??")
        
    
    return filt_far, start

# def generate_probe (dna_dict, far_start, probe_start, min_len = 18, max_len = 24)
#     if not probe_start:
#         probe_start = dna_dict["position"] - far_start

#         1 - 1500 - 3000