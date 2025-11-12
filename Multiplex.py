import pandas as pd
import primer3
import random



def loop_part(primer, allele, big_list, best_primers, Homodimer_Max):
    curr_prime_probs = 0

    for i in big_list:
        if big_list[allele] == big_list[i]:
            pass
        elif (primer3.calcHeterodimer(big_list[allele][primer], big_list[i][best_primers[i]]).dg > Homodimer_Max):
            curr_prime_probs += 1

    return curr_prime_probs
        
def loop_whole(big_list, best_primers, allele_probs_count, Homodimer_Max):
    for i in range(len(big_list)):
        allele_probs_count[i] = loop_part(best_primers[i], big_list[i], big_list, best_primers, Homodimer_Max)

def find_best_primer(allele, big_list, best_primers, allele_probs_count, Homodimer_Max):
    num_primes_allele_has = best_primers[allele].count()

    for primer in range(num_primes_allele_has):
        probs_found = loop_part(primer, allele, big_list, best_primers, Homodimer_Max)

        if probs_found == 0:
            best_primers[allele] = primer
            return
        elif (probs_found < allele_probs_count[allele]):
            allele_probs_count[allele] = probs_found




def multiplex_list(big_list: list[list[dict]], Homodimer_Max = 5):
    list_size = big_list.count()
    allele_probs_count = [0] * list_size
    best_primers = [0] * list_size
    loop_whole(big_list, best_primers, allele_probs_count, Homodimer_Max)

    times_through = 0
    while(max(allele_probs_count) > 1 or times_through < list_size):
        worst = allele_probs_count.index(max(allele_probs_count))
        find_best_primer(worst, big_list, best_primers, allele_probs_count, Homodimer_Max)
        loop_whole(big_list, best_primers, allele_probs_count, Homodimer_Max)
        allele_probs_count[worst] = allele_probs_count[worst] * -1 # I invert it so if won't trigger the while condition, but we still can see who many probblems it had
        times_through +=1
     













def add(A, B):
    return A + B

def Check_Multiplex_Compatibility(primer_pairs: pd.DataFrame, heterodimer_max: float = 50.0): # -> [{score: 32.2, combination: (P1, P2, P3, P4)}, {...}]:
    """
    Check primer sets for multiplex compatibility.
    TODO: Enhance for multiple SNPs.
    - Extend to check cross-SNP interactions (current checks within SNP).
    - Optimize for large primer sets using batch dimer calculations.
    """

    def _get_snp_ids(primer_pairs: pd.DataFrame) -> list:
        """
        Extract unique SNP IDs from primer pairs DataFrame.
        """
        # Assuming 'snpId' is a column in the primer_pairs DataFrame
        return primer_pairs['snpId'].unique().tolist()
    
    def _get_primer_pair_list(primer_pairs: pd.DataFrame) -> list:
        """
        Extract primer pairs from DataFrame.
        """
        # return [tuple(item) for item in primer_pairs[['Forward', 'Reverse']].values]
        # TODO: add index from 0 to len(primer_pairs)
        return []
    
    def _calculate_hetero_compatibility_score(heterodimer_results) -> float:
        """
        Calculate compatibility score based on hairpin, homodimer, and heterodimer results.
        """
        # Example scoring logic (to be replaced with actual logic)
        # score = 0.0
        # for result in hairpin_results + homodimer_results + heterodimer_results:
        #     if result['tm'] < heterodimer_max:
        #         score += 1.0
        # return score / len(hairpin_results + homodimer_results + heterodimer_results)
        return random.uniform(0, 100)  # Placeholder for actual score calculation
    
    def _calculate_homo_compatibility_score(hairpin_results, homodimer_results, heterodimer_results) -> float:
        """
        Calculate compatibility score based on hairpin, homodimer, and heterodimer results.
        """
        # Example scoring logic (to be replaced with actual logic)
        # score = 0.0
        # for result in hairpin_results + homodimer_results + heterodimer_results:
        #     if result['tm'] < heterodimer_max:
        #         score += 1.0
        # return score / len(hairpin_results + homodimer_results + heterodimer_results)
        return random.uniform(0, 100)  # Placeholder for actual score calculation

    ### 1. Prepare the data
    pp_list = _get_primer_pair_list(primer_pairs) # [(primer_forward, primer_reverse), ...]

    ### 2. Define compatibility score calculation for all primer pairs
    homo_table = [0] * len(pp_list) # score of homodimer&hairpin. as score goes higher, it makes worse effect on pcr
    hetero_table = [[0] * len(pp_list)] * len(pp_list) # score of heterodimer. as score goes higher, it makes worse effect on pcr

    for i in range(len(pp_list)):
        primer1_f, primer1_r = pp_list[i]

        # calculate hairpin&homodimer score
        hairpin_results = []
        hairpin_results.append(primer3.calcHairpin(primer1_f.sequence))
        hairpin_results.append(primer3.calcHairpin(primer1_r.sequence))
        hairpin_results.append(primer3.calcHairpin(primer2_f.sequence))
        hairpin_results.append(primer3.calcHairpin(primer2_r.sequence))

        homodimer_results = []
        homodimer_results.append(primer3.calcHomodimer(primer1_f.sequence))
        homodimer_results.append(primer3.calcHomodimer(primer1_r.sequence))
        homodimer_results.append(primer3.calcHomodimer(primer2_f.sequence))
        homodimer_results.append(primer3.calcHomodimer(primer2_r.sequence))

        homo_table[i] = _calculate_homo_compatibility_score(hairpin_results, homodimer_results)

        for allele in range(i + 1, len(pp_list)):
            if(i == allele):
                continue
            primer2_f, primer2_r = pp_list[allele]
            # Calculate heterodimer score
            heterodimer_results = []
            heterodimer_results.append((primer1_f.sequence, primer1_r.sequence))
            heterodimer_results.append(primer3.calcHeterodimer(primer1_f.sequence, primer2_f.sequence))
            heterodimer_results.append(primer3.calcHeterodimer(primer1_f.sequence, primer2_r.sequence))
            heterodimer_results.append(primer3.calcHeterodimer(primer1_r.sequence, primer2_f.sequence))
            heterodimer_results.append(primer3.calcHeterodimer(primer1_r.sequence, primer2_r.sequence))
            heterodimer_results.append(primer3.calcHeterodimer(primer2_f.sequence, primer2_r.sequence))

            result = _calculate_hetero_compatibility_score(heterodimer_results)
            hetero_table[i][allele] = result
            hetero_table[allele][i] = result

    ### 3. calculate compatibility scores for all **combinations** of primer pairs
    # prepare the all combinations of primer pairs
    from itertools import product
    snp_ids = _get_snp_ids(primer_pairs)  # Extract unique SNP IDs from primer pairs
    grouped = []
    for snp in snp_ids:
        grouped.append([pp for pp in pp_list if pp['snpId'] == snp])

    combinations = list(product(*grouped)) # [(P1, P2, P3, P4), (P5, P2, P3, P4)] if len(snp_ids) == 4 

    # calculate compatibility scores for all combinations
    scores = []
    for comb in combinations: # comb is (P1, P2, P3, P4)
        score = 0
        for i in range(len(comb)):
            pp1 = comb[i]
            score += homo_table[pp1.index]
            for allele in range(i + 1, len(comb)):
                pp2 = comb[allele]
                score += hetero_table[pp1.index][pp2.index]
        scores.append(score)
    
    # sort combinations by their score
    zipped = list(zip(scores, combinations))
    zipped_sorted = sorted(zipped, key=lambda x: x[0])
    sorted_scores, sorted_combinations = zip(*zipped_sorted)
    
    ### Packing result
    result = []
    for i in range(len(sorted_scores)):
        dict = {}
        dict['score'] = sorted_scores[i]
        dict['combination'] = sorted_combinations[i]
        result.append(dict)
    return result

    