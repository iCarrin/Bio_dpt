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
            #Locked Then Normal (+XY)
            'aA': HSG(H=-7193,  S=-19.723, G=-1.09),
            'aC': HSG(H=-7269,  S=-18.336, G=-1.56),
            'aG': HSG(H=-7536,  S=-18.387, G=-1.84),
            'aT': HSG(H=-4918,  S=-12.943, G=-0.89),
            'cA': HSG(H=-7451,  S=-18.380, G=-1.72),
            'cC': HSG(H=-5904,  S=-11.904, G=-2.30),
            'cG': HSG(H=-9815,  S=-23.491, G=-2.50),
            'cT': HSG(H=-7092,  S=-16.825, G=-1.95),
            'gA': HSG(H=-5038,  S=-11.656, G=-1.37),
            'gC': HSG(H=-10160, S=-24.651, G=-2.65),
            'gG': HSG(H=-10844, S=-26.580, G=-2.54),
            'gT': HSG(H=-8612,  S=-22.327, G=-1.63),
            'tA': HSG(H=-7246,  S=-19.738, G=-1.14),
            'tC': HSG(H=-6307,  S=-15.515, G=-1.51),
            'tG': HSG(H=-10040, S=-25.744, G=-2.00),
            'tT': HSG(H=-6372,  S=-16.902, G=-1.13),
            #Normal Then Locked (X+Y)
            'Aa': HSG(H=-6908,  S=-18.135, G=-1.40),
            'Ac': HSG(H=-5510,  S=-11.824, G=-1.83),
            'Ag': HSG(H=-9000,  S=-22.826, G=-1.88),
            'At': HSG(H=-5384,  S=-13.537, G=-1.19),
            'Ca': HSG(H=-7142,  S=-18.333, G=-1.40),
            'Cc': HSG(H=-5937,  S=-12.335, G=-2.24),
            'Cg': HSG(H=-10876, S=-27.918, G=-2.17),
            'Ct': HSG(H=-9471,  S=-25.070, G=-1.69),
            'Ga': HSG(H=-7756,  S=-19.302, G=-1.74),
            'Gc': HSG(H=-10725, S=-25.511, G=-2.78),
            'Gg': HSG(H=-8943,  S=-20.833, G=-2.51),
            'Gt': HSG(H=-9035,  S=-22.742, G=-1.96),
            'Ta': HSG(H=-5609,  S=-16.019, G=-0.58),
            'Tc': HSG(H=-7591,  S=-19.031, G=-1.70),
            'Tg': HSG(H=-6335,  S=-15.537, G=-1.56),
            'Tt': HSG(H=-5574,  S=-14.149, G=-1.21),
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


# big = "C"
# plus = f"+{big}"
# small = big.lower() 
# for i in range(12,40,2):
#     first_half = i //2
#     # print(f"{i}	{big*first_half}{plus*3}{big*first_half}")

#     print(round(calc_lna_tm((big*first_half + small*3 + big*first_half),200),1))



# Seq1	ATCGTATGCTGATCTGATCAC
# Seq2	AACGTTTGCTGATCTAATCAC
# Seq3	TTCGCATGCCGATCAGATCAC

    