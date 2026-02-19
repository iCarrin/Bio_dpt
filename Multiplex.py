from Primer_Classes import *
import primer3
from itertools import combinations,product
from Primer_functions import generate_matching_primers

def func(primer1,primer2,heterodimer_max = 50.0,tm_max=40,memo={}):
    if (key:=tuple(sorted((primer1,primer2)))) in memo:
        return memo[key]
    result=primer3.calc_heterodimer(primer1,primer2)
    ans=result.dg > heterodimer_max or result.tm > tm_max
    memo[key]= ans
    return ans

def multiplex_far(close_primers, snp_list):
    '''
    Multiplex the far primers against the close primers and any far primers that have already succeeded
    
    '''
    all_good_fars = [] 
    done_snpid = set()
    for close_primer in close_primers: # find each close primer a far primer match
        if close_primer.snpID in done_snpid:
            continue
        else:
            done_snpid.add(close_primer.snpID)
        for far in generate_matching_primers(close_primer, snp_list):#loop every possible primer given
            for primer in (close_primers + all_good_fars): #compare this far primer against all close and already found far primers
                 # calculate it's heterodimer value every other primer far and close
                if func(far.sequence, primer.sequence):#it only fails if it has a delta gibbs lower than the max AND the dimer will happen at temp that will bother us
                    break #stop checking early
            else: #only if it got through the check everything loop
                all_good_fars.append(far) # we add it to the final list
                # and tell the while loop that this close primer has found it's soul mate
                break #the far primer is found stop searching the far primer list
    return all_good_fars



#this is basically a glorified heterodimer filter. Glorified because it has to check all options against all others 
def multiplex_close(big_list: list[list[Primer]]):
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
            primer_index = golden_primers[allele_index]
        return big_list[allele_index][primer_index]
         
        

    def get_heterodimer(left, right, leftPrimer = None):
        """
        This function was made to cut down on the noise that comes from calling the primer calc_heterodimer function
        """ 
        return func(get_primer(left, leftPrimer).sequence,get_primer(right).sequence)

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
            probs_found = 0                                                           # primer is only used in this function 
            for i in range(list_size):
                if i != allele and get_heterodimer(allele, i, primer):
                    probs_found+=1
                    if alleles_prob_count[i] == 0:
                        alleles_prob_count[i] = 1
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
                
    
    """that's the end of the internal functions"""
    #loop the whole list (using the combo list to cut the N^2 time in half)
    for left, right in allele_combos:
        #log all of the problems
        if get_heterodimer(left, right):
            #since we're using the combo list we update both locations when finding a problem
            alleles_prob_count[left] += 1                     
            alleles_prob_count[right] += 1
    #every time we go over something we will change it's remaining problems to be negative so we it won't trigger the while loop

    while((res:=max(enumerate(alleles_prob_count),key=lambda x: x[1]))[1]>0):
        #find the allele that's causing the most problems and start with it first.
        find_best_primer(res[0])                   
        alleles_prob_count[res[0]] = -alleles_prob_count[res[0]]
        
    fighting_alleles = []
    #this makes a list of where the failures still are
    remaining_stubborn = [i for i in range(list_size) if alleles_prob_count[i] < 0]
    #this makes that into a combo list so that the final loop only has to look at those locations
    fight_combos = combinations(remaining_stubborn,2)
    

    for left, right in fight_combos:#
        fighting_alleles.append((f"{(l:=get_primer(left)).snpID} : {l.allele}", 
                                f"{(r:=get_primer(right)).snpID} : {r.allele}"))
        
    out_list = [get_primer(i) for i in range(list_size)]

    # multiplexing took : 0:00:05.508837
    return out_list, fighting_alleles

def multiplex_close2(big_list: list[list[Primer]], heterodimer_max = 50.0):
    print("start")
    primer_sequences = [
    [p.sequence for p in allele_list]
    for allele_list in big_list
    ]
    print("finish primer seq")
    calc_heterodimer = primer3.calc_heterodimer

    hetero_cache = {}

    for a1, a2 in combinations(range(len(primer_sequences)), 2):
        print(a1,a2)
        primers1 = primer_sequences[a1]
        primers2 = primer_sequences[a2]
        for s1, s2 in product(primers1, primers2):
            key = (s1, s2) if s1 < s2 else (s2, s1)
            if key not in hetero_cache:
                result = calc_heterodimer(s1, s2)
                hetero_cache[key] = (
                    result.dg > heterodimer_max or result.tm > 40
                )
                            
    print("finish hetero_cashe")
    conflicts = {}
    for a1, a2 in combinations(range(len(primer_sequences)), 2):
            print(a1,a2)
            primers1 = primer_sequences[a1]
            primers2 = primer_sequences[a2]
            for (i, s1), (j, s2) in product(enumerate(primers1), enumerate(primers2)):
                    key = (s1, s2) if s1 < s2 else (s2, s1)
                    if hetero_cache[key]:
                        conflicts.setdefault((a1,i), set()).add((a2,j))
                        conflicts.setdefault((a2,j), set()).add((a1,i))

    print("finished conficts")
    selected = [0]*len(big_list)

    for allele in range(len(big_list)):

        best = None
        best_score = float("inf")

        for primer_i in range(len(big_list[allele])):

            score = 0

            for other_allele in range(allele):

                other_primer = selected[other_allele]

                if (allele,primer_i) in conflicts:
                    if (other_allele,other_primer) in conflicts[(allele,primer_i)]:
                        score += 1

            if score < best_score:
                best_score = score
                best = primer_i

        selected[allele] = best
    print("finish selected")
    out_list = [
    big_list[a][p]
    for a,p in enumerate(selected)
    ]   
    return out_list,conflicts

def multiplex_close3(big_list: list[list[Primer]], heterodimer_max = 50.0):
    print("start")
    num_alleles = len(big_list)

    # Create integer masks for each primer
    # conflicts_mask[allele][primer_index] = integer where bit j set if conflicts with primer j of other alleles
    conflicts_mask = [
        [0]*len(allele) for allele in big_list
    ]

    # Flatten all primers to assign unique IDs for bits
    primer_id_map = {}  # (allele, primer_index) -> unique id
    id_counter = 0
    for a, allele in enumerate(big_list):
        for i in range(len(allele)):
            primer_id_map[(a,i)] = id_counter
            id_counter += 1

    # print("start hetero cashe")
    # Precompute sequences
    primer_sequences = [[p.sequence for p in allele] for allele in big_list]
    calc_heterodimer = primer3.calc_heterodimer
    hetero_cache = {}

    for a1, a2 in combinations(range(num_alleles), 2):
        print(a1,a2)
        primers1 = primer_sequences[a1]
        primers2 = primer_sequences[a2]

        for (i, s1), (j, s2) in product(enumerate(primers1), enumerate(primers2)):
            key = (s1, s2) if s1 < s2 else (s2, s1)
            if key not in hetero_cache:
                result = calc_heterodimer(s1, s2)
                hetero_cache[key] = (result.dg > heterodimer_max or result.tm > 40)

            if hetero_cache[key]:
                id1 = primer_id_map[(a1,i)]
                id2 = primer_id_map[(a2,j)]
                conflicts_mask[a1][i] |= 1 << id2
                conflicts_mask[a2][j] |= 1 << id1
    # print("finish hetero cashe")
    selected = [0]*num_alleles
    current_mask = 0

    for a, allele in enumerate(big_list):
        best = None
        min_conflicts = float("inf")
        
        for i in range(len(allele)):
            conflicts = (current_mask & conflicts_mask[a][i]).bit_count()
            if conflicts < min_conflicts:
                min_conflicts = conflicts
                best = i
                
        selected[a] = best
        # Update current mask with the newly selected primer
        current_mask |= 1 << primer_id_map[(a,best)]
    
    out_list = [big_list[a][p] for a, p in enumerate(selected)]

    fighting_alleles = []
    for a1, a2 in combinations(range(num_alleles), 2):
        p1 = selected[a1]
        p2 = selected[a2]
        id2 = primer_id_map[(a2,p2)]
        if conflicts_mask[a1][p1] & (1 << id2):
            fighting_alleles.append((
                f"{big_list[a1][p1].snpID} : {big_list[a1][p1].allele}",
                f"{big_list[a2][p2].snpID} : {big_list[a2][p2].allele}"
            ))
    return out_list,fighting_alleles