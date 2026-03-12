from Primer_Classes import Primer, Probe
import primer3
from itertools import combinations,chain
from Primer_functions import generate_matching_primers

def check_heterodimer(primer1,primer2,heterodimer_max = 50.0,tm_max=40,memo=None):
    if memo is None: memo = {}
    if (key:=tuple(sorted((primer1,primer2)))) in memo:
        return memo[key]
    result=primer3.calc_heterodimer(primer1,primer2)
    ans=result.dg > heterodimer_max * 1000 and result.tm > tm_max
    memo[key]= ans
    return ans

def multiplex_far(close_primers, snp_list):
    '''
    Multiplex the far primers against the close primers and any far primers that have already succeeded
    
    '''
    pos_far_primers = [] 
    neg_far_primers = [] 
    temp_pos = None
    temp_neg = None
    done_snpid = set()
    snps_to_remove = set()
    for close_primer in close_primers: # find each close primer a far primer match
        if close_primer.snpID in done_snpid:
            continue
        done_snpid.add(close_primer.snpID)

        temp_pos = None
        temp_neg = None
        flipped = False
        count = 0
        while True:
            for direction in ['forward', 'reverse']:
                count +=1
                if count >2:
                    print(f"{count} times through!!!!!!!")
               
                for far in generate_matching_primers(close_primer, snp_list, direction, flipped):#loop every possible primer given 
                    #loop 4
                    for primer in chain(close_primers, pos_far_primers, neg_far_primers): #compare this far primer against all close and already found far primers
                        # calculate it's heterodimer value every other primer far and close
                        if check_heterodimer(far.sequence, primer.sequence):#it only fails if it has a delta gibbs lower than the max AND the dimer will happen at temp that will bother us
                            break #stop checking early
                    else: #only if it got through the check everything loop

                        if direction == 'forward': 
                            temp_pos = far # we add it to the final list
                        else: 
                            temp_neg = far # we add it to the final list
                                # and tell the while loop that this close primer has found it's soul mate
                        break #the far primer is found stop searching the far primer list
                 
            if temp_pos and temp_neg :# we've run through every primer we can try
                break
            elif not flipped:
                print("we're in the end game")
                flipped = True
            else:
                print("all hope is lost")
                snps_to_remove.add(close_primer.snpID)
                print(f"snpID: {close_primer.snpID} (and all associated alleles) didn't make the list. It couldn't find far primers to that could work")
                break

        pos_far_primers.append(temp_pos) 
        neg_far_primers.append(temp_neg)
    
    done_snpid = [p for p in done_snpid if p not in snps_to_remove]

    if len(done_snpid) == len(neg_far_primers) == len(pos_far_primers):
        pass
    else:
        print(f'close primers after removal {len(close_primers)}')
        print(f'neg fars {len(neg_far_primers)}')
        print(f'pos fars {len(pos_far_primers)}')
    return (neg_far_primers, pos_far_primers, done_snpid, snps_to_remove)



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
        return check_heterodimer(get_primer(left, leftPrimer).sequence,get_primer(right).sequence)

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

