import primer3
import logging
from itertools import combinations,chain
from Bio.Seq import Seq
from typing import Generator


from Primer_Classes import Primer, Probe, FilterFail
from pdfoutput import create_output_json

class Multiplexer():
    def __init__(self,snp_df,alleles,**kwarg):
        self.snp_df=snp_df
        self.known_alleles = alleles
        self.primer3=primer3.thermoanalysis.ThermoAnalysis()
        self.primer3.set_thermo_args(
            mv_conc=kwarg.get("mv_conc",50), dv_conc=kwarg.get("dv_conc",3), dntp_conc=kwarg.get("dntp_conc",0.8), 
            dna_conc=kwarg.get("dna_conc",200), dmso_conc=kwarg.get("dmso_conc",0.0), 
            dmso_fact=kwarg.get("dmso_fact",0.0), formamide_conc=kwarg.get("formamide_conc",0.0),
            annealing_temp_c=kwarg.get("annealing_temp_c",-10.0), temp_c=kwarg.get("temp_c",37), 
            tm_method=kwarg.get("tm_method","santalucia"), salt_corrections_method=kwarg.get("salt_corrections_method","owczarzy")) 

        self.desired_tm=kwarg.get("desired_tm",60.0)
        self.diff=kwarg.get("diff",2)
        self.homodimer_goal=kwarg.get("homodimer_goal",-3.0)
        self.hairpin_goal=kwarg.get("hairpin_goal",-3.0)
        self.target_gc=kwarg.get("target_gc",50.0)
        self.heterodimer_max = kwarg.get("heterodimer_max",50.0)
        self.primer_dimer_distance=kwarg.get("primer_dimer_distance",20)
        self.min_probe_len=kwarg.get("min_probe_len",12)
        self.max_probe_len=kwarg.get("max_probe_len",20)
        self.min_primer_len=kwarg.get("min_primer_len",18)
        self.max_primer_len=kwarg.get("max_primer_len",24)
        self.min_primer_dist=kwarg.get("min_primer_dist",50)
        self.max_primer_dist=kwarg.get("max_primer_dist",75)
        self.pdfoutput_precision=kwarg.get("pdfoutput_precision",2)
        self.bad_alleles = []
        self.logger = logging.getLogger(__name__)

    def _check_heterodimer(self,primer1,primer2,memo={}):
        if (key:=tuple(sorted((primer1,primer2)))) in memo:
            return memo[key]
        result=self.primer3.calc_heterodimer(primer1,primer2)
        ans=result.dg > self.heterodimer_max * 1000 and result.tm > (self.desired_tm + self.diff - self.primer_dimer_distance )
        memo[key]= ans
        return ans 
    
    def _multiplex_primers(self,snp_site_probes):
        '''
        Multiplex the far primers against the close primers and any far primers that have already succeeded
        
        '''
        pos_far_primers = [] 
        neg_far_primers = [] 
        temp_pos = None
        temp_neg = None
        done_snpid = set()

        for probe in snp_site_probes: # find each probe a set of primers
            if probe.snpID in done_snpid:
                continue
            
            temp_pos = None
            temp_neg = None

            
            while True:
                for direction in ['forward', 'reverse']:
                    for far in self.generate_matching_primers(probe, direction):#loop every possible primer given 
                        for entry in chain(snp_site_probes, pos_far_primers, neg_far_primers): #compare this far primer against all close and already found far primers
                            # calculate it's heterodimer value every other primer far and close
                            if self._check_heterodimer(far.sequence, entry.sequence):#it only fails if it has a delta gibbs lower than the max AND the dimer will happen at temp that will bother us
                                break #stop checking early
                        else: #only if it got through the check everything loop
                            
                            if direction == 'forward': 
                                temp_pos = far # we add it to the final list
                            else: 
                                temp_neg = far # we add it to the final list
                                    # and tell the while loop that this close primer has found it's soul mate
                            done_snpid.add(probe.snpID)
                            break #the far primer is found stop searching the far primer list
                    
                if temp_pos and temp_neg :# we've run through every primer we can try
                    break
                else:
                    print(f"all hope is lost for probe: {probe.snpID} allele: {probe.allele}")

                    for entry in snp_site_probes:
                        if entry.snpID == probe.snpID:
                            self.bad_alleles.append({"snpid":probe.snpID, "allele":entry.allele, "reason":"probe couldn't find primers"})
                    print(f"snpID: {probe.snpID} (and all associated alleles) didn't make the list. It couldn't find far primers to that could work")
                    break
            if temp_pos is not None and temp_neg is not None:
                pos_far_primers.append(temp_pos) 
                neg_far_primers.append(temp_neg)
        

        if not (len(done_snpid) == len(neg_far_primers) == len(pos_far_primers)):
            print(f'probes after removal {len(snp_site_probes)}')
            print(f'neg fars {len(neg_far_primers)}')
            print(f'pos fars {len(pos_far_primers)}')
        


        return (neg_far_primers, pos_far_primers)
    
    def _multiplex_probes(self,big_list: list[list[Primer]]):
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

        def _get_next_probe(allele_index, primer_index=None):
            """
            this function was made just to cut down on the noise
            that comes from referencing a dict in a list in a list
            """
            if primer_index is None:
                primer_index = golden_primers[allele_index]
            return big_list[allele_index][primer_index]
            
        def _get_heterodimer(left, right, leftPrimer = None):
            """
            This function was made to cut down on the noise that comes from calling the primer calc_heterodimer function
            """ 
            return self._check_heterodimer(_get_next_probe(left, leftPrimer).sequence,_get_next_probe(right).sequence)

        def find_best_probe(allele):
            """
            This function takes a problem allele and cycles through it's primers until it finds one that isn't a problem or runs out
            of primers in which case it uses the least problematic one it could find
            """
            #this is to make sure we don't run off the end of the list of dictionaries
            num_probes = len(big_list[allele])

            #we start at the original primer even though it's failings are the reason we're in this function in the first place.
            #The reason for that iscase the other primer that formed a hetero dimer was fixed before coming to this one. In that case
            #we check the first primer again, see that the other primer it fought with was straitened out, and end the loop
            for probe in range(num_probes):
                #AI gave me this. I wanted to clean up a if not statement and it up classed me out of town.
                            #add every index from the list that isn't the allele it's self, and that make a hetero dimer
                probs_found = 0                                                           # primer is only used in this function 
                for i in range(list_size):
                    if i != allele and _get_heterodimer(allele, i, probe):
                        probs_found+=1
                        if alleles_prob_count[i] == 0:
                            alleles_prob_count[i] = 1
                                                                                        # to make sure the primers are looping
                #if it's better
                if probs_found < alleles_prob_count[allele]:
                    #update who the current best primer is 
                    golden_primers[allele] = probe 
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
            if _get_heterodimer(left, right):
                #since we're using the combo list we update both locations when finding a problem
                alleles_prob_count[left] += 1                     
                alleles_prob_count[right] += 1
        #every time we go over something we will change it's remaining problems to be negative so we it won't trigger the while loop
        while((res:=max(enumerate(alleles_prob_count),key=lambda x: x[1]))[1]>0):
            #find the allele that's causing the most problems and start with it first.
            find_best_probe(res[0])                   
            alleles_prob_count[res[0]] = -alleles_prob_count[res[0]]
        
        #this makes a list of where the failures still are
        remaining_stubborn = [i for i in range(list_size) if alleles_prob_count[i] < 0]

        for i in remaining_stubborn:
            bad_guy = _get_next_probe(i)
            self.bad_alleles.append({"snpid":bad_guy.snpID, "allele":bad_guy.allele, "reason":"probe on probe heterodimer"})

        #this makes that into a combo list so that the final loop only has to look at those locations
        # fight_combos = combinations(remaining_stubborn,2)
        
        # for left, right in fight_combos:#
        #     self.bad_alleles.append({"snpid":(l:=_get_next_probe(left)).snpID, "allele":l.allele}) 
        #     self.bad_alleles.append({"snpid":(r:=_get_next_probe(right)).snpID, "allele":r.allele})
            
        probes_that_work = [_get_next_probe(i) for i in range(list_size) ]

        return probes_that_work
    
    def _make_probes(self, seq, snp_id, allele, direction="forward",index=0) -> list[Probe]:
        probes = []
        if len(seq) >= self.min_probe_len:   #len(seq) should just be max_len. It only won't be if the seq length is less than min_len (if the sequence is only 10 long then we'll trigger this)
            len_of_the_flank = (len(seq)-self.min_probe_len)//2
            cof = [0,0,0,0,0] #Count Of Failure
            rff = ["TM too low", "Tm Too high", "Homodimer", "Hairpin", "Probes"] # Reasons For Failure
            for length in range(len_of_the_flank):#possible bug if the forward mismatch is smaller than the minimum length
                                                #length is 0-max_len
                trimmed = seq[length: -length if length != 0 else None]#this is assuming that the seq given is already the maximum length. 
                # If given a crazy long string it will start at "length" and give the rest of the string
                #take this part out of the loop, so we can have one dictionary that says the SNP ID and ALLELE and Direction, 
                #and then a list in that dictionary of sequence and lengths. Storing the name over and over seems redundant ID   
                #These take one off the left and right to see if it gets us anywhere out of 34 it saved an extra 3
                for segment in [trimmed, trimmed[:-1], trimmed[1:]]:
                    cof[4] += 1 
                    try:                                           
                        probes.append(Probe(snp_id, allele, segment, direction, self.primer3, self.desired_tm, self.diff, self.homodimer_goal, self.hairpin_goal,self.target_gc))
                    except FilterFail as e:
                        match e.fail_type:
                            case "lower Tm":
                                cof[0] += 1
                            case "upper Tm":
                                cof[1] += 1
                            case "homodimer":
                                cof[2] += 1
                            case "hairpin":
                                cof[3] += 1
                        pass

        else:
    
            print(f"The length of your {direction} primer wasn't long enough. \nYou needed one at least {self.min_probe_len} long and it ended up only being {len(seq)}")
            self.logger.warning(f"The length of your {direction} primer {snp_id} allele {allele} wasn't long enough. \nYou needed one at least {self.min_probe_len} long and it ended up only being {len(seq)}")
            
        '''probe error messages (simple and extensive)'''
        # br = cof.index(max(cof[:-1]))
        # print(f'snp: {snp_id} allele: {allele} {cof[4]} {rff[4]}, {cof[0]+cof[1]+cof[2]+cof[3]} misses, {cof[br]} were {rff[br]}')
        # print(f'snp: {snp_id} allele: {allele} {cof[4]} {rff[4]}, {cof[0]+cof[1]+cof[2]+cof[3]} misses, {cof[0]} were {rff[0]}, {cof[1]} were {rff[1]}, {cof[2]} were {rff[2]}, {cof[3]} were {rff[3]}')
        return probes

    def _make_primers(self,seq, snp_id, allele, direction="forward") -> Generator[Primer]: 
    
        if len(seq) >= self.min_primer_len:   #len(seq) should just be max_len. It only wont be if the seq length is less than max_len (if the sequence is only 10 long then we'll trigger this)
            for length in range(self.max_primer_len-self.min_primer_len+1):#possible bug if the forward mismatch is smaller than the minimum length
                                                #length is 0-max_len
                trimmed = seq[length:]#this is assuming that the seq given is already the maximum length. 
                # If given a crazy long string it will start at "length" and give the rest of the string
                #take this part out of the loop, so we can have one dictionary that says the SNP ID and ALLELE and Direction, 
                #and then a list in that dictionary of sequence and lengths. Storing the name over and over seems redundant ID
            
                try:                                     #these need to be user controlled inputs
                    yield Primer(snp_id, allele, trimmed, direction,self.primer3 ,self.desired_tm, self.diff, self.homodimer_goal, self.hairpin_goal,self.target_gc)
                except FilterFail as e:
                    # print(e)
                    pass
                    # logging.error(f"{snp_id} allele: {allele} had no primers that passed the filtering")
        else:
            print(f"The length of your {direction} primer wasn't long enough. \nYou needed one at least {self.min_primer_len} long and it ended up only being {len(seq)}")
            self.logger.warning(f"The length of your {direction} primer {snp_id} allele {allele} wasn't long enough. \nYou needed one at least {self.min_primer_len} long and it ended up only being {len(seq)}")
            raise ValueError
        # return primers

    def _generate_allele_specific_probes(self) -> list[list[Probe]]:
        all_probes = []
        
        def _make_allele_probes_list(snp_id, allele, sequence, snp_pos,i):
            #this is a lazy way to put an LNA triplet into the middle of our probe sequence
            lna_sequence = sequence[:snp_pos-1] + (sequence[snp_pos-1:snp_pos+2]).lower() + sequence[snp_pos+2:]
            known_allele_list = self.known_alleles[snp_id]['allele_str']
      
            #add a check here for length so that make primers doesn't have to.
            flank_len = self.max_probe_len - self.max_probe_len//2#this gives the longer half
            #this gives us the longer half so in case we need to drop a g form the 5' end it balances better
            forward = lna_sequence[snp_pos - flank_len : snp_pos+(flank_len)+1]#this gets the largest segment.   
            forward_probes = self._make_probes(forward, snp_id, allele,index=i) if allele != "G" or "A" not in known_allele_list else []

            reverse = str(Seq(lna_sequence[snp_pos-flank_len : snp_pos+flank_len+1]).reverse_complement()) #creates a Biopython sequence, gets the reverse complement, and converts is back to a string                                                                  #adding a check to see if 
            reverse_probes = self._make_probes(reverse,snp_id, allele, "reverse",index=i) if allele != "C" or "A" not in known_allele_list else []

            probes = (forward_probes + reverse_probes)
            probes.sort(key=lambda x: x.rank)
            return probes
        
        for i,snp in enumerate(self.snp_df):
            allele_probes = _make_allele_probes_list(snp["snpID"], snp["allele"], snp["sequence"], snp["position"],i)
            if allele_probes:
                all_probes.append(allele_probes)
            else:
                # print(f"snp: {snp["snpID"]} allele: {snp["allele"]} didn't make the cut")
                self.bad_alleles.append({"snpid":snp["snpID"], "allele":snp["allele"], "reason":"probe wouldn't generate"})
        
        return all_probes

    def generate_matching_primers(self,primer_king, direction, primer_start = 0): 
        """
            Generate matching primers for best allele-specific probe.
        """

        if len(self.snp_df[0]['sequence']) < 1+(self.min_primer_dist+self.min_primer_len)*2:
            raise Exception("your sequence is so short it won't allow for even the shortest of primers at the minimum distance to be made" \
            "Lower your min distance or have the API call for a longer string")
        
        snp_dict = {}   
        for snp in self.snp_df:
            if snp['snpID'] == primer_king.snpID: # we only search for a matching SNP because we aren't even using the allele section    
                snp_dict = snp
                break   
        else:
            raise Exception(f"there is no matching entry in the Json list for {primer_king.snpID}")
          
        middle = snp_dict['position']
        whole_sequence = snp_dict['sequence']
        #get the far sequence and reverse complement if necessary
        if direction == "forward":
            far_sequence = whole_sequence[middle+self.min_primer_dist:middle+self.max_primer_dist]#everything from snp (middle) plus start dist cutting off at max if necessary (end is exclusive so +1)
            '''                                  _______  
                                                |______\_
                                                |________\ what we make eventually 
            _________________________|___________=====================("===" means reverse complement DNA)
                                    SNP                  All the way to max_dist gets reverse complemented
                                        min_dist away
            
            '''       
        elif direction == "reverse":
            far_sequence = whole_sequence[middle-self.max_primer_dist:middle-self.min_primer_dist-1]

            '''
                                    SNP
            ____________________________|____________________________
                        __/______|           
            What makin'./_________|  
                            Go min_dist away 
            and take up to Max_dist to use in primer making   
            '''
        else:
            raise Exception("no direction given. How did we get here?")
        
        #get some far primers for primer_king, but only the best (use strict mode). if the filter fails the while loop should try again farther down the line
        #if the ones that pass don't pass the heterodimer then primer_start should be updated and the whole thing tried again.
        while(True):
            # the reverse complement flips the sequence each time so we can iterate over each one the same way. Researchers expect DNA to come in the forward
            # format even if it was reversed in real life, so no need to flip it back. The note that it should be reversed is enough.
            # each one will walk back from the right side 
            trial_snip = far_sequence[-(self.max_primer_len+primer_start) : -primer_start or None] # this takes a chunk to feed into a primer generator    
            try:
                yield from self._make_primers(trial_snip,primer_king.snpID, primer_king.allele, direction)
                primer_start += 1
            except ValueError as e:
                self.logger.warning(f'{primer_king.snpID} allele {primer_king.allele} had no useable far primers')
                "you've tried every possible primer. What have you done??"
                break
                
    def main(self):
        """
        Main function to generate, filter, rank, pair, and export primers.
        TODO: Add comprehensive error handling and logging.
        - Log progress and errors to a file.
        - Add input validation for rsIDs and output format.
        TODO: Test with real SNP data.
        - Validate output with biological experts.
        - Benchmark performance for large SNP sets.
        """
    

        primers = self._generate_allele_specific_probes()
        if len(primers)==0:
            raise Exception(f"{len(primers)=}")
        best_probes= self._multiplex_probes(primers)
        if len(best_probes)==0:
            raise Exception(f"{len(best_probes)=}")
        forward, reverse, = self._multiplex_primers(best_probes)
        tre= best_probes + forward+reverse ####this might not be the right course of action to add the probes into the mix
        tre2= best_probes ####this might not be the right course of action to add the probes into the mix
        tre.sort(key=lambda x : x.snpID)
        # tre2.sort(key=lambda x : x.snpID)
        print("these alleles didn't make it")
        print(self.bad_alleles)
        # print(forword,reverse,center,center_reject,sep="\n")
        # create_output_json(best_primers,"final_primer.pdf",(1000,1000),self.pdfoutput_precision)
        # create_output_json(fights,"final_fights.pdf",(1000,1000),self.pdfoutput_precision)
        # create_output_json(self.bad_alleles,"final_bad_combined.pdf",(120,350),self.pdfoutput_precision)
        create_output_json(tre,"final_combined.pdf",(1000,1000),self.pdfoutput_precision)
        # create_output_json(tre2,"final_probe_combined.pdf",(1000,1000),self.pdfoutput_precision)
        # create_output_json(forword,"final_forword.pdf",(1000,1000),self.pdfoutput_precision)
        # create_output_json(reverse,"final_reverse.pdf",(1000,1000),self.pdfoutput_precision)
        # create_output_json(center,"final_center.pdf",(1000,1000),self.pdfoutput_precision)
        # create_output_json(center_reject,"final_center_reject.pdf",(900,1000),self.pdfoutput_precision)


                    


