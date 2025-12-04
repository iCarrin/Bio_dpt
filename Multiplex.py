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


def multiplex_far (Generate_Matching_Primers, get_heterodimer):
        # compare each heterodimer, if that primer plays well with other primers then save. Try new one if it doesn't works.