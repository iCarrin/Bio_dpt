import primer3

def calc_gc_content(sequence: str):
    gc_total = sequence.count('G') + sequence.count('C')
    return round(gc_total/len(sequence), 2)

class FilterFail(Exception):
    def __init__(self, id:str, allele:str, fail_type:str, result = None):
        self.fail_type = fail_type
        super().__init__(f'{id} allele: {allele} failed the {fail_type} {result}')
    

class Primer():

    def __init__(self, snp_id, allele, sequence, direction,thermo, desired_tm: float = 60.0, diff: float = 2.0, homodimer_goal: float = -3.0, hairpin_goal: float = -3.0, target_gc = 50.0):
        
        self.tm = thermo.calc_tm(sequence)
        if self.tm < (desired_tm-diff):
            # print(self.tm)
            raise FilterFail(snp_id, allele, "lower Tm", self.tm)
        if self.tm > (desired_tm+diff):
            raise FilterFail(snp_id, allele, "upper Tm", self.tm)
        
        __homodimer_thermo = thermo.calc_homodimer(sequence)
        self.homodimer_tm = __homodimer_thermo.tm
        self.homodimer_dg = __homodimer_thermo.dg
        if self.homodimer_dg < homodimer_goal* 1000 and self.homodimer_tm > desired_tm-diff-20:
            raise FilterFail(snp_id, allele, "homodimer")
        
        __hairpin_thermo = thermo.calc_hairpin(sequence)
        self.hairpin_tm = __hairpin_thermo.tm
        self.hairpin_dg = __hairpin_thermo.dg
        if self.hairpin_dg < hairpin_goal* 1000 and self.hairpin_tm > desired_tm-diff-20:
            raise FilterFail(snp_id, allele, "hairpin")
        
        self.snpID = snp_id
        self.allele = allele
        self.direction = direction
        self.sequence = sequence #this is the primer length
        self.length = len(sequence)
        self.gc_content = calc_gc_content(sequence)
        self.target_gc = target_gc
        self.rank = abs(self.tm - desired_tm) + abs(self.gc_content - target_gc)
   

        
class Probe(Primer):
    # self, snp_id, allele, primer,  direction, desired_tm: float = 60.0, diff: float = 3.0, homodimer_goal: float = -3.0, hairpin_goal: float = -3.0
    def __init__(self, snp_id, allele, sequence, direction,thermo ,desired_tm: float = 60.0, diff: float = 2.0, homodimer_goal: float = -3.0, hairpin_goal: float = -3.0, target_gc = 50.0):
        super().__init__(snp_id, allele, sequence, direction, thermo,desired_tm, diff, homodimer_goal, hairpin_goal, target_gc)
    def to_list(self,percision):
        return [ val if  type(val:=self.__getattribute__(i))!=float else round(val,percision)  for i in vars(self)]
    

        