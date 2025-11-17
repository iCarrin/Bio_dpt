import primer3
from itertools import combinations


#this is basically a glorified heterodimer filter. Glorified because it has to check all options against all others 
def multiplex_list(big_list: list[list[dict]], heterodimer_max = 50.0):
    """
    This is the heterodimer close primer filter. 
    The thought was that if we filter the close primers to where they like each other than the far primers will have
    the slimmest possible list to filter against. All close primers are have their heterodimer checked against all other 
    close primers. If there is a heterodimer one of the primers is cycled through it's primer list until either there
    are no more problems or the end of the primer list is reached, in which case the primer with the leaset problems is
    saved. Once all problem primers have been cycled through the 

    """
    # time saved for keeping track of improvements
 

    #instead of looping through a list of lists (N^2) we find all combinations. This avoids looking at combinations already tried ((n(n-1))/2)
    allele_combos = combinations(range(len(big_list)),2)
    #this is called enough in for loops that we save it as a variable
    list_size = len(big_list)
    
    #this keeps track of how many problems each allele has
    alleles_prob_count = [0] * list_size
    #this keeps track of the current primer that is working for the allele corresponding to the indices 
    golden_primers = [0] * list_size


    def get_primer(allele_index, primer_index=None):
        """
        this function was made just to cut down on the noise
        that comes from referencing a dict in a list in a list
        """
        if primer_index is None:
            primer_index = allele_index
        return big_list[allele_index][golden_primers[primer_index]]
         
        

    def get_heterodimer(left, right, leftPrimer = None):
        """
        This fucntion was made to cut down on the noise that comes from calling the primer calc_heterodimer function
        """
        return primer3.calc_heterodimer(get_primer(left, leftPrimer)['primer_sequence'], get_primer(right)['primer_sequence']).dg


    def find_best_primer(allele):
        """
        This function takes a problem allele and cycles through it's primers until it finds one that isn't a problem or runs out
        of primers in which case it uses the least problematic one it could find
        """
        #this is to make sure we don't run off the end of the list of dictionaries
        num_primes = len(big_list[allele])

        #we start at the original primer even though it's failings are the reason we're in this function in the first place.
        #The reason for that iscase the other primer that formed a hetero dimer was fixed before coming to this one. In that case
        #we check the first primer again, see that the other primer it fought with was straitened out, and end the loop
        for primer in range(num_primes):
            #AI gave me this. I wanted to clean up a if not statement and it up classed me out of town.
                        #add every index from the list that isn't the allele it's self, and that make a hetero dimer
            fight_list = [i for i in range(list_size) if i != allele and get_heterodimer(allele, i, primer) > heterodimer_max]
            probs_found = len(fight_list)                                                           # primer is only used in this function 
                                                                                                    # to make sure the primers are looping
         
            #if it's better
            if probs_found < alleles_prob_count[allele]:
                #update who the current best primer is 
                golden_primers[allele] = primer 
                # update it's lower problem count
                alleles_prob_count[allele] = probs_found
                #and if those problems were 0 we can go home early
                if probs_found == 0:
                    break
            #at the end (will get skipped if probs was 0) we check to make sure that we didn't make any new problems
            for i in fight_list:
                if alleles_prob_count[i] == 0:
                    alleles_prob_count[i] = 1
    """that's the end of the internal functions"""

    #loop the whole list (using the combo list to cut the N^2 time in half)
    for left, right in allele_combos:
        #log all of the problems
        if get_heterodimer(left, right) > heterodimer_max:
            #since we're using the combo list we update both locations when finding a problem
            alleles_prob_count[left] += 1                     
            alleles_prob_count[right] += 1

 
    #every time we go over something we will change it's remaining problems to be negative so we it won't trigger the while loop
    while(max(alleles_prob_count)>0):
        #find the allele that's causing the most problems and start with it first.
        worst = alleles_prob_count.index(max(alleles_prob_count))
        find_best_primer(worst)                   
        alleles_prob_count[worst] = -alleles_prob_count[worst]
        
    fighting_alleles = []
    #this makes a list of where the failures still are
    remaining_stubborn = [i for i in range(list_size) if alleles_prob_count[i] < 0]
    #this makes that into a combo list so that the final loop only has to look at those locations
    fight_combos = combinations(remaining_stubborn,2)
    

    for left, right in fight_combos:#
        fighting_alleles.append((f"{get_primer(left)['snpID']} : {get_primer(left)['allele']}", 
                                f"{get_primer(right)['snpID']} : {get_primer(right)['allele']}"))
        
 

    
    out_list = [get_primer(i) for i in range(list_size)]
  

   
    # multiplexing took : 0:00:05.508837
    return out_list, fighting_alleles









# def add(A, B):
#     return A + B

# def Check_Multiplex_Compatibility(primer_pairs: pd.DataFrame, heterodimer_max: float = 50.0): # -> [{score: 32.2, combination: (P1, P2, P3, P4)}, {...}]:
#     """
#     Check primer sets for multiplex compatibility.
#     TODO: Enhance for multiple SNPs.
#     - Extend to check cross-SNP interactions (current checks within SNP).
#     - Optimize for large primer sets using batch dimer calculations.
#     """

#     def _get_snp_ids(primer_pairs: pd.DataFrame) -> list:
#         """
#         Extract unique SNP IDs from primer pairs DataFrame.
#         """
#         # Assuming 'snpId' is a column in the primer_pairs DataFrame
#         return primer_pairs['snpId'].unique().tolist()
    
#     def _get_primer_pair_list(primer_pairs: pd.DataFrame) -> list:
#         """
#         Extract primer pairs from DataFrame.
#         """
#         # return [tuple(item) for item in primer_pairs[['Forward', 'Reverse']].values]
#         # TODO: add index from 0 to len(primer_pairs)
#         return []
    
#     def _calculate_hetero_compatibility_score(heterodimer_results) -> float:
#         """
#         Calculate compatibility score based on hairpin, homodimer, and heterodimer results.
#         """
#         # Example scoring logic (to be replaced with actual logic)
#         # score = 0.0
#         # for result in hairpin_results + homodimer_results + heterodimer_results:
#         #     if result['tm'] < heterodimer_max:
#         #         score += 1.0
#         # return score / len(hairpin_results + homodimer_results + heterodimer_results)
#         return random.uniform(0, 100)  # Placeholder for actual score calculation
    
#     def _calculate_homo_compatibility_score(hairpin_results, homodimer_results, heterodimer_results) -> float:
#         """
#         Calculate compatibility score based on hairpin, homodimer, and heterodimer results.
#         """
#         # Example scoring logic (to be replaced with actual logic)
#         # score = 0.0
#         # for result in hairpin_results + homodimer_results + heterodimer_results:
#         #     if result['tm'] < heterodimer_max:
#         #         score += 1.0
#         # return score / len(hairpin_results + homodimer_results + heterodimer_results)
#         return random.uniform(0, 100)  # Placeholder for actual score calculation

#     ### 1. Prepare the data
#     pp_list = _get_primer_pair_list(primer_pairs) # [(primer_forward, primer_reverse), ...]

#     ### 2. Define compatibility score calculation for all primer pairs
#     homo_table = [0] * len(pp_list) # score of homodimer&hairpin. as score goes higher, it makes worse effect on pcr
#     hetero_table = [[0] * len(pp_list)] * len(pp_list) # score of heterodimer. as score goes higher, it makes worse effect on pcr

#     for i in range(len(pp_list)):
#         primer1_f, primer1_r = pp_list[i]

#         # calculate hairpin&homodimer score
#         hairpin_results = []
#         hairpin_results.append(primer3.calcHairpin(primer1_f.sequence))
#         hairpin_results.append(primer3.calcHairpin(primer1_r.sequence))
#         hairpin_results.append(primer3.calcHairpin(primer2_f.sequence))
#         hairpin_results.append(primer3.calcHairpin(primer2_r.sequence))

#         homodimer_results = []
#         homodimer_results.append(primer3.calcHomodimer(primer1_f.sequence))
#         homodimer_results.append(primer3.calcHomodimer(primer1_r.sequence))
#         homodimer_results.append(primer3.calcHomodimer(primer2_f.sequence))
#         homodimer_results.append(primer3.calcHomodimer(primer2_r.sequence))

#         homo_table[i] = _calculate_homo_compatibility_score(hairpin_results, homodimer_results)

#         for allele in range(i + 1, len(pp_list)):
#             if(i == allele):
#                 continue
#             primer2_f, primer2_r = pp_list[allele]
#             # Calculate heterodimer score
#             heterodimer_results = []
#             heterodimer_results.append((primer1_f.sequence, primer1_r.sequence))
#             heterodimer_results.append(primer3.calcHeterodimer(primer1_f.sequence, primer2_f.sequence))
#             heterodimer_results.append(primer3.calcHeterodimer(primer1_f.sequence, primer2_r.sequence))
#             heterodimer_results.append(primer3.calcHeterodimer(primer1_r.sequence, primer2_f.sequence))
#             heterodimer_results.append(primer3.calcHeterodimer(primer1_r.sequence, primer2_r.sequence))
#             heterodimer_results.append(primer3.calcHeterodimer(primer2_f.sequence, primer2_r.sequence))

#             result = _calculate_hetero_compatibility_score(heterodimer_results)
#             hetero_table[i][allele] = result
#             hetero_table[allele][i] = result

#     ### 3. calculate compatibility scores for all **combinations** of primer pairs
#     # prepare the all combinations of primer pairs
#     from itertools import product
#     snp_ids = _get_snp_ids(primer_pairs)  # Extract unique SNP IDs from primer pairs
#     grouped = []
#     for snp in snp_ids:
#         grouped.append([pp for pp in pp_list if pp['snpId'] == snp])

#     combinations = list(product(*grouped)) # [(P1, P2, P3, P4), (P5, P2, P3, P4)] if len(snp_ids) == 4 

#     # calculate compatibility scores for all combinations
#     scores = []
#     for comb in combinations: # comb is (P1, P2, P3, P4)
#         score = 0
#         for i in range(len(comb)):
#             pp1 = comb[i]
#             score += homo_table[pp1.index]
#             for allele in range(i + 1, len(comb)):
#                 pp2 = comb[allele]
#                 score += hetero_table[pp1.index][pp2.index]
#         scores.append(score)
    
#     # sort combinations by their score
#     zipped = list(zip(scores, combinations))
#     zipped_sorted = sorted(zipped, key=lambda x: x[0])
#     sorted_scores, sorted_combinations = zip(*zipped_sorted)
    
#     ### Packing result
#     result = []
#     for i in range(len(sorted_scores)):
#         dict = {}
#         dict['score'] = sorted_scores[i]
#         dict['combination'] = sorted_combinations[i]
#         result.append(dict)
#     return result

    