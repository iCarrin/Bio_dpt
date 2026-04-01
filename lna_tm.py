import math
from collections import namedtuple
HSG = namedtuple('HSG', ['H','S', 'G'])

dublets = {
            #         cal/mol     cal/k*mol
            'AA': HSG(H=-7900,    S=-22.2,   G=-1.02),
            'TT': HSG(H=-7900,    S=-22.2,   G=-1.02),
            'AT': HSG(H=-7200,    S=-20.4,   G=-0.88),
            'TA': HSG(H=-7200,    S=-21.3,   G=-0.73),
            'CA': HSG(H=-8500,    S=-22.7,   G=-1.60),
            'TG': HSG(H=-8500,    S=-22.7,   G=-1.60),
            'GT': HSG(H=-8400,    S=-22.4,   G=-1.38),
            'AC': HSG(H=-8400,    S=-22.4,   G=-1.38),
            'CT': HSG(H=-7800,    S=-21.0,   G=-1.43),
            'AG': HSG(H=-7800,    S=-21.0,   G=-1.43),
            'GA': HSG(H=-8200,    S=-22.2,   G=-1.16),
            'TC': HSG(H=-8200,    S=-22.2,   G=-1.16),
            'CG': HSG(H=-10600,   S=-27.2,   G=-3.36),
            'GC': HSG(H=-9800,    S=-24.4,   G=-2.09),
            'GG': HSG(H=-8000,    S=-19.9,   G=-2.28),
            'CC': HSG(H=-8000,    S=-19.9,   G=-2.28),
            #locked then norma
            'aA': HSG(H=707,   S=2.477,   G=-0.092),
            'aC': HSG(H=1131,   S=4.064,   G=-0.122),
            'aG': HSG(H=264,   S=2.613,   G=-0.561),
            'aT': HSG(H=2282,   S=7.457,   G=-0.007),
            'cA': HSG(H=1049,   S=4.320,   G=-0.270),
            'cC': HSG(H=2096,   S=7.996,   G=-0.457),
            'cG': HSG(H=785,   S=3.709,   G=-0.332),
            'cT': HSG(H=708,   S=4.175,   G=-0.666),
            'gA': HSG(H=3162,   S=10.544,  G=-0.072),
            'gC': HSG(H=-360,  S=-0.251,  G=-0.414),
            'gG': HSG(H=-2844,  S=-6.680,  G=-0.700),
            'gT': HSG(H=-212,  S=0.073,   G=-0.194),
            'tA': HSG(H=-46,  S=1.562,   G=-0.563),
            'tC': HSG(H=1893,   S=6.685,   G=-0.208),
            'tG': HSG(H=-1540,  S=-3.044,  G=-0.548),
            'tT': HSG(H=1528,   S=5.298,   G=-0.130),
            #Normal then locke
            'Aa': HSG(H=992,   S=4.065,   G=-0.396),
            'Ac': HSG(H=2890,   S=10.576,  G=-0.390),
            'Ag': HSG(H=-1200,  S=-1.826,  G=-0.603),
            'At': HSG(H=1816,   S=6.863,   G=-0.309),
            'Ca': HSG(H=1358,   S=4.367,   G=0.046),
            'Cc': HSG(H=2063,   S=7.565,   G=-0.404),
            'Cg': HSG(H=-276,  S=-0.718,  G=-0.003),
            'Ct': HSG(H=-1671,  S=-4.070,  G=-0.409),
            'Ga': HSG(H=444,   S=2.898,   G=-0.437),
            'Gc': HSG(H=-925,  S=-1.111,  G=-0.535),
            'Gg': HSG(H=-943,  S=-0.933,  G=-0.666),
            'Gt': HSG(H=-635,  S=-0.342,  G=-0.520),
            'Ta': HSG(H=1591,   S=5.281,   G=0.004),
            'Tc': HSG(H=609,   S=3.169,   G=-0.396),
            'Tg': HSG(H=2165,   S=7.163,   G=-0.106),
            'Tt': HSG(H=2326,   S=8.051,   G=-0.212),
            #Locked then Locked         
            'aa': HSG(H=-9991,  S=-27.175, G=-1.57),
            'ac': HSG(H=-11389, S=-28.963, G=-2.44),
            'ag': HSG(H=-12793, S=-31.607, G=-3.07),
            'at': HSG(H=-14703, S=-40.750, G=-2.12),
            'ca': HSG(H=-14177, S=-35.498, G=-3.11),
            'cc': HSG(H=-15399, S=-36.375, G=-4.15),
            'cg': HSG(H=-14558, S=-35.239, G=-3.65),
            'ct': HSG(H=-15737, S=-41.218, G=-2.92),
            'ga': HSG(H=-13959, S=-35.097, G=-3.03),
            'gc': HSG(H=-16109, S=-40.738, G=-3.48),
            'gg': HSG(H=-13022, S=-29.673, G=-3.82),
            'gt': HSG(H=-17361, S=-45.858, G=-3.12),
            'ta': HSG(H=-10318, S=-26.108, G=-2.19),
            'tc': HSG(H=-9166,  S=-21.535, G=-2.48),
            'tg': HSG(H=-10046, S=-22.591, G=-3.08),
            'tt': HSG(H=-10419, S=-27.683, G=-1.83)
        }

def divalent_to_monovalent(divalent_mM, dntp_mM):
    """Convert divalent cation concentration to monovalent equivalent (mM)."""
    if divalent_mM == 0:
        dntp_mM = 0
    if divalent_mM < dntp_mM:
        divalent_mM = dntp_mM  # dNTPs chelate Mg2+
    return 120 * math.sqrt(divalent_mM - dntp_mM)


# dna_conc_nM = dna_conc = Oligo Conc
# K_mM = mv_conc = Na+ Conc
# divalent_mM = dv_conc = Mg++ Conc
# dntp_mM = dntp_conc = dNTPs Conc

# 1 mM = 1000 µM = 1,000,000 nM
# 1 µM = 1000 nM
# 1 nM = 0.001 µM = 0.000001 mM

def calc_lna_tm(seq, dna_conc_nM, K_mM=50, divalent_mM=3, dntp_mM=.8,
                dmso_conc=0, dmso_fact=0, formamide_conc=0):
    H = 0
    S = 0
   


    for i in range(len(seq) - 1):
        window = seq[i:i+2]

        H += dublets[window].H
        S += dublets[window].S

    for end in [0, -1]:
        if seq[end] in ["A", "T"]:
            H += 2300
            S += 4.1
        else:
            H += 100
            S += -2.8

    seq_len = len(seq)
    GC_count = sum(1 for b in seq if b in ('G', 'C'))

    # Add divalent -> monovalent conversion
    K_mM = K_mM + divalent_to_monovalent(divalent_mM, dntp_mM)

    # SantaLucia salt correction modifies entropy
    S = S + 0.368 * (seq_len - 1) * math.log(K_mM / 1000.0)

    # Compute Tm (non-symmetric assumed; add symmetry check if needed)
    ### CHECK INTO THE SELF COMPLEMENTRY 
    Tm = H / (S + 1.987 * math.log((dna_conc_nM / 1e9)/4)) - 273.15

    # DMSO and formamide corrections
    Tm -= dmso_conc * dmso_fact
    Tm += (0.453 * GC_count / seq_len - 2.88) * formamide_conc

    return Tm

# my testing results https://docs.google.com/spreadsheets/d/1vhDBcd0xlyCc8ZkBS0_ncm_Wgm1kFlaCr6CGU6G4WYU/edit?usp=sharing
# the thoughts are that the salting method is Owczazry online and I'm using SantaLucia also my LNA tables could just be out of date































