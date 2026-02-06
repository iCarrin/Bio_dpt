import primer3
def calc_gc_content(sequence: str):
    gc_total = 0
    for nucleotide in sequence:
        if nucleotide == 'G' or nucleotide == 'C':
            gc_total += 1
    return gc_total/len(sequence)

class FilterFail(Exception):
    '''just passes it silently. No need to cause any stink'''
    pass

class Primer():

    def __init__(self, snp_id, allele, primer, direction, desired_tm: float = 60.0, diff: float = 3.0, homodimer_goal: float = -3.0, hairpin_goal: float = -3.0):
        
        self.tm = primer3.bindings.calc_tm(primer)
        if self.tm < (desired_tm-diff):
            raise FilterFail
        if self.tm > (desired_tm+diff):
            raise FilterFail
        self.__homodimer_thermo = primer3.bindings.calc_homodimer(primer)
        if self.__homodimer_thermo.dg > homodimer_goal* 1000 or self.__homodimer_thermo.tm < desired_tm-diff-20:
            raise FilterFail
        self.__hairpin_thermo = primer3.bindings.calc_hairpin(primer)
        if self.__hairpin_thermo.dg > hairpin_goal* 1000 or self.__hairpin_thermo.tm < desired_tm-diff-20:
            raise FilterFail
        self.snpID = snp_id
        self.allele = allele
        self.direction = direction
        self.primer_sequence = primer #this is the primer length
        self.length = len(primer)
        self.gc_content = calc_gc_content(primer)
        
    
    # not sure how to make this better than just calling hairpin_thermo.tm
    
    def get_homomdimer_tm(self):
        return self.__homodimer_thermo.tm
    def get_homomdimer_dg(self):
        return self.__homodimer_thermo.dg
    def get_hairpin_tm(self):
        return self.__hairpin_thermo.tm
    def get_hairpin_tm(self):
        return self.__hairpin_thermo.dg

class Probe(Primer):
    # def __init__(self, snp_id, allele, seq, length, direction):
    #     super().__init__(snp_id, allele, seq, length, direction)
    pass

        