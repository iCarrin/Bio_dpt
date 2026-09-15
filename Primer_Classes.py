import primer3
from lna_tm_shiny import calc_tm_with_lna
from calc_gc import calc_gc

class FilterFail(Exception):
    def __init__(self, id:str, allele:str, fail_type:str, result = None):
        self.fail_type = fail_type
        super().__init__(f'{id} allele: {allele} failed the {fail_type} {result}')
    

class Primer():

    def _calc_tm(self, sequence, thermo):
        return thermo.calc_tm(sequence)

    def __init__(self, snp_id, allele, sequence, direction, thermo, desired_tm: float, diff: float, homodimer_goal: float, hairpin_goal: float, target_gc: float, salt_corrections_method):

        self.tm = self._calc_tm(sequence, thermo, salt_corrections_method)
        if self.tm < (desired_tm-diff):
            raise FilterFail(snp_id, allele, "lower Tm", self.tm)
        if self.tm > (desired_tm+diff):
            raise FilterFail(snp_id, allele, "upper Tm", self.tm)
        
        __homodimer_thermo = thermo.calc_homodimer(sequence)
        self.homodimer_tm = __homodimer_thermo.tm
        self.homodimer_dg = __homodimer_thermo.dg
        if self.homodimer_dg < homodimer_goal* 1000 and self.homodimer_tm > desired_tm-diff-20:
            raise FilterFail(snp_id, allele, "homodimer")
        self.homodimer_dg = round(self.homodimer_dg, 2)

        __hairpin_thermo = thermo.calc_hairpin(sequence)
        self.hairpin_tm = __hairpin_thermo.tm
        self.hairpin_dg = __hairpin_thermo.dg
        if self.hairpin_dg < hairpin_goal* 1000 and self.hairpin_tm > desired_tm-diff-20:
            raise FilterFail(snp_id, allele, "hairpin")
        
        self.hairpin_dg = round(self.hairpin_dg, 2)
        self.snpID = snp_id
        self.allele = allele
        self.direction = direction
        self.sequence = sequence #this is the primer length
        self.length = len(sequence)
        self.gc_content = calc_gc(sequence)
        self.target_gc = target_gc
        self.rank = abs(self.tm - desired_tm) + abs(self.gc_content - target_gc)
    def to_list(self,percision):
        return [ val if  type(val:=self.__getattribute__(i))!=float else round(val,percision)  for i in vars(self)]
   

        
class Probe(Primer):
    def _calc_tm(self, sequence, thermo, salt_corrections_method):
                                                                                            
        return calc_tm_with_lna(sequence, thermo.dna_conc, thermo.mv_conc, thermo.dv_conc, thermo.dntp_conc, thermo.dmso_conc, thermo.dmso_fact, thermo.formamide_conc, salt_corrections_method)
    
    def __init__(self, snp_id, allele, sequence, direction, thermo, desired_tm: float, diff: float, homodimer_goal: float, hairpin_goal: float , target_gc: float, salt_corrections_method):
        super().__init__(snp_id, allele, sequence, direction, thermo, desired_tm, diff, homodimer_goal, hairpin_goal, target_gc, salt_corrections_method)
  
# only related to primer3py
# dmso_conc     =
# dmso_fact     = 
# formamide_conc=

# dna_conc_nM = dna_conc = Oligo Conc
# K_mM = mv_conc = Na+ Conc
# divalent_mM = dv_conc = Mg++ Conc
# dntp_mM = dntp_conc = dNTPs Conc
    

        