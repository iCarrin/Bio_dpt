from Bio.Seq import Seq
import primer3
import re
import json
from Primer_Classes import *
from lna_tm_shiny import *

class TempStuff():
    def __init__(self):
        self.dna_conc  = 200
        self.mv_conc = 50
        self.dv_conc = 3
        self.dntp_conc = 0.8
        self.dmso_conc = 0.0
        self.dmso_fact = 0.0
        self.formamide_conc = 0.0
        self.salt_correction_method = "owczarzy"
temp_guy = TempStuff()
billy = primer3.thermoanalysis.ThermoAnalysis()
billy.set_thermo_args(
            mv_conc = 200,
            dv_conc = 50,
            dntp_conc = 3,
            dna_conc = 0.8,
            dmso_conc = 0.0,
            dmso_fact = 0.0,
            formamide_conc = 0.0,
            salt_correction_method = "owczarzy"
)
 

# print(calc_tm_with_lna("AAAAAAAAAAAA", temp_guy))
string = "AGCAACTAGTGACTGACTT"

location = len(string)//2
# print(string[location])
# print(string[:-location-1] + "[" + "Q" + "]" + string[-location:])


asdf = "AGCAACTAGGGACTGACTT"
asdfg = []
for i in range(len(asdf)):
    asdfg.append(asdf[-(i+1)])

asdfgh = "".join(asdfg)
print(asdfgh)

# 0123456789012345678
# AAGTCAGTCACTAGTTGCT
# TTCAGTCAGGGATCAACGA
# 8765432109876543210