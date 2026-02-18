import primer3
def calc_gc_content(sequence: str):
    gc_total = 0
    for nucleotide in sequence:
        if nucleotide == 'G' or nucleotide == 'C':
            gc_total += 1
    return round(gc_total/len(sequence), 2)

class FilterFail(Exception):
    '''just passes it silently. No need to cause any stink'''
    def __init__(self, fail_type:str):
        self.fail_type = fail_type
        super().__init__(f'filter faild: {fail_type}')
    

class Primer():

    def __init__(self, snp_id, allele, sequence, direction, desired_tm: float = 60.0, diff: float = 3.0, homodimer_goal: float = -3.0, hairpin_goal: float = -3.0):
        
        self.tm = primer3.bindings.calc_tm(sequence)
        if self.tm < (desired_tm-diff):
            raise FilterFail("lower Tm")
        
        if self.tm > (desired_tm+diff):
            raise FilterFail("upper Tm")
        
        __homodimer_thermo = primer3.bindings.calc_homodimer(sequence)
        self.homodimer_tm = __homodimer_thermo.tm
        self.homodimer_dg = __homodimer_thermo.dg
        if self.homodimer_dg < homodimer_goal* 1000 and self.homodimer_tm > desired_tm-diff-20:
            raise FilterFail("homodimer")
        
        __hairpin_thermo = primer3.bindings.calc_hairpin(sequence)
        self.hairpin_tm = __hairpin_thermo.tm
        self.hairpin_dg = __hairpin_thermo.dg
        if self.hairpin_dg < hairpin_goal* 1000 and self.hairpin_tm > desired_tm-diff-20:
            raise FilterFail("hairpin")
        
        self.snpID = snp_id
        self.allele = allele
        self.direction = direction
        self.sequence = sequence #this is the primer length
        self.length = len(sequence)
        self.gc_content = calc_gc_content(sequence)
        self.rank = ""
    

    def set_rank(self, rank):
        self.rank = rank
   

        
class Probe(Primer):
    # self, snp_id, allele, primer,  direction, desired_tm: float = 60.0, diff: float = 3.0, homodimer_goal: float = -3.0, hairpin_goal: float = -3.0
    def __init__(self, snp_id, allele, sequence, direction, desired_tm: float = 70.0, diff: float = 3.0, homodimer_goal: float = -3.0, hairpin_goal: float = -3.0):
        self.sequence = sequence
        if self.sequence[0] == "G":
            raise FilterFail("started with G")
        super().__init__(snp_id, allele, sequence, direction, desired_tm, diff, homodimer_goal, hairpin_goal)
    

        