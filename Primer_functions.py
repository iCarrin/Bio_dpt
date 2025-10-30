
import re 
from Bio.Seq import Seq
import logging
from typing import Dict
import primer3
from collections.abc import Callable
# from itertools import filter

# the order is 
# generate allele specific (and evaluate)
# generate matching (and evaluate)
# filter
# rank


# # Sequence -> Sequence
# def Reverse_Complement(sequence: str): #  -> str
    
#     return str(Seq(sequence).reverse_complement())
# This function was unnecessary


def introduce_mismatch(primer_sequence: str) -> str:
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

    # Simple mismatch rules (purine↔purine, pyrimidine↔pyrimidine)
    mismatch_rules = {
        "A": "G", "G": "A",
        "C": "T", "T": "C"
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


def rank_primers(primers: list[dict], target_tm = 62.5, target_gc = 50, optimism = 5) -> list[dict]:
    """
        Rank primers based on Tm proximity to 62.5Â°C and GC content.
        TODO: Refine ranking criteria.
        - Consider weighting Tm vs. GC scores.
        - Add user-configurable ranking metrics.
        """
    for primer in primers:
        primer["tm_score"] = abs(primer["tm"] - target_tm)
        primer["gc_score"] = abs(primer["gc_content"] - target_gc)
        primer["score"] = primer["tm_score"] + primer["hairpin_dg"] + primer["homodimer_dg"] + primer["gc_score"] 


    primers.sort(key=lambda x: x["score"])

    
    return primers[:optimism]


def filter_little(filter_name: str, old_list: list[dict], filter_function):
    # fail_count is only 0 or 1, but we'll add it to a total to see what's bugging out
    fail_count = 0
    #the new list off of a filter of the passed in list. The filter will pick any allele that worked
    new_list = list(filter(filter_function, old_list))
    #if after the filter we don't have anything left
    if not new_list:
        #just use the old list
        new_list = old_list
        #print a statement for debugging purposes
        print(f"{old_list[0]["snpID"]}: {filter_name} filter failed; using previous list.")
        #the fail count goes up
        fail_count = 1
    #return the list that hopefully was able to filter and weather or not it failed
    return(new_list, fail_count)


def filter_one_list_soft(allele_list: list[dict],
                         desired_tm: float = 60.0,
                         diff: float = 3.0,
                         homodimer_goal: float = 3.0,
                         hairpin_goal: float = 3.0) -> (list[dict], list[int]):
    
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
    allele_pltm, ltm_fail_count = filter_little("tm Min", allele_list, lambda x : x["tm"] >= (desired_tm - diff))
    # tm < max
    allele_phtm, htm_fail_count = filter_little("tm Max", allele_pltm, lambda x : x["tm"] <= (desired_tm + diff))
    # min < homodimer < max
    allele_phomo, homo_fail_count = filter_little("homodimer", allele_phtm, lambda x : (homodimer_goal*-1) < x["homodimer_dg"] < homodimer_goal)
    #min < hairpin < max
    allele_phair, hair_fail_count = filter_little("hairpin", allele_phomo, lambda x : (hairpin_goal*-1) < x["hairpin_dg"] < hairpin_goal)
    
    # total_fails = [f"low temp fails: {ltm_fail_count}", f"high temp fails: {htm_fail_count}", f"homodimer fails: {homo_fail_count}", f"hairpin fails: {hair_fail_count}"]
    # print(total_fails)
    total_fails_ints = [ltm_fail_count, htm_fail_count, homo_fail_count, hair_fail_count]
    
    return (allele_phair, total_fails_ints)


def generate_allele_specific_primers(snps_list: list[dict], min_len: int = 18, max_len: int = 28) -> list[list[list[dict]]]:
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

    for snp_dict in snps_list:
        #why pass this in seperately when it's already in the dict?
        snp_id = snp_dict["snpID"]
        allele = snp_dict["allele"]
        sequence = snp_dict["sequence"]
        snp_pos = snp_dict["position"]

        forward = sequence[snp_pos - max_len :snp_pos+1]#this gets the largest segment.   
        forward_mismatch = introduce_mismatch(forward)
    
        reverse = str(Seq(sequence[snp_pos:snp_pos+max_len+1]).reverse_complement()) #creates a Biopython sequence, gets the reverse complement, and converts is back to a string
        reverse_mismatch = introduce_mismatch(reverse)

        this_allele_primers = (make_primers(forward_mismatch, min_len, max_len, snp_id, allele))\
                            + (make_primers(reverse_mismatch, min_len, max_len, snp_id, allele, "reverse"))# this make one list of dictionaries for both flanking directions
        all_primers.append(this_allele_primers) #this adds this list to the larger list
    return all_primers # this will return a list of lists of dictionaries. Each allele is a list. 
#[[snp1 allele1 dictionaries],[snp1 allele2 dictionaries],[snp2 allele2 dictionaries],[snp2 allele2 dictionaries]] each snp and allele are on the same level.
   

def make_primers(seq, min_len, max_len, snp_id, allele, direction="forward") -> list[dict]: 
    seq_length = len(seq)
    primers = []
    if seq_length >= min_len:
        for length in range(max_len-min_len):#possible bug if the forward mismatch is smaller than the minimum length
            
            primers.append({
                "snpID": snp_id,
                "allele": allele,
                "primer_sequence": seq[length:], #this is the trimmed length
                "direction": direction,
                "length": seq_length-length,
                "tm" : primer3.bindings.calc_tm(length),
                "gc_content" : calc_gc_content(length),
                "hairpin_dg" : primer3.bindings.calc_hairpin(length).dg,
                "homodimer_dg" : primer3.bindings.calc_homodimer(length).dg
            })
            
    else:
        print(f"The length of your forward primer wasn't long enough. \nYou needed one at least {min_len} long and it ended up only being {seq_length}")
    return primers
