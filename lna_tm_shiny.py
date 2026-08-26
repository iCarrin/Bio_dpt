import math
import primer3
from calc_gc import calc_gc 
from collections import namedtuple

hs = namedtuple("hs", ["H", "S"])
# calc temp of both sides, add the temp of the LNAs to one
dublets = {
#Standard DNA doublets. No missmatch
'AA/TT': hs(-7.9, -22.2), 'AT/TA': hs(-7.2, -20.4), 'TA/AT': hs(-7.2, -21.3), 'TT/AA': hs(-7.9, -22.2),
'CA/GT': hs(-8.5, -22.7), 'GT/CA': hs(-8.4, -22.4), 'CT/GA': hs(-7.8, -21.0), 'GA/CT': hs(-8.2, -22.2),
'CG/GC': hs(-10.6, -27.2), 'GC/CG': hs(-9.8, -24.4), 'GG/CC': hs(-8.0, -19.9), 'CC/GG': hs(-8.0, -19.9),
'AC/TG': hs(-8.4, -22.4), 'TG/AC': hs(-8.4, -22.4),'AG/TC': hs(-7.8, -21.0), 'TC/AG': hs(-7.8, -21.0),
#Locked pairs to DNA Pairs no missmatch (Table 2, "Consecutive LNA Modifications")
'aa/TT': hs(-9.991, -27.175), 'ac/TG': hs(-11.389, -28.963), 'ag/TC': hs(-12.793, -31.607), 'at/TA': hs(-14.703, -40.750),
'ca/GT': hs(-14.177, -35.498), 'cc/GG': hs(-15.399, -36.375), 'cg/GC': hs(-14.558, -35.239), 'ct/GA': hs(-15.737, -41.218),
'ga/CT': hs(-13.959, -35.097), 'gc/CG': hs(-16.109, -40.738), 'gg/CC': hs(-13.022, -29.673), 'gt/CA': hs(-17.361, -45.858),
'ta/AT': hs(-10.318, -26.108), 'tc/AG': hs(-9.166, -21.535), 'tg/AC': hs(-10.046, -22.591), 'tt/AA': hs(-10.419, -27.683),
#+A·A Mismatch (Table 3, full parameters)
'aa/AT': hs(-3.826, -13.109), 'ac/AG': hs(-2.367, -7.322), 'ag/AC': hs(-4.849, -13.007), 'at/AA': hs(-5.049, -17.514),
'aa/TA': hs(-4.229, -15.160), 'ca/GA': hs(-5.878, -17.663), 'ga/CA': hs(-8.558, -23.976), 'ta/AA': hs(2.074, 3.446),
#+C·C Mismatch
'ca/CT': hs(2.218, 4.750), 'cc/CG': hs(1.127, 1.826), 'cg/CC': hs(-10.903, -32.025), 'ct/CA': hs(-2.053, -10.517),
'ac/TC': hs(1.065, -1.403), 'cc/GC': hs(-9.522, -27.024), 'gc/CC': hs(-4.767, -14.897), 'tc/AC': hs(4.114, 9.258),
#+G·G Mismatch
'ga/GT': hs(-2.920, -9.387), 'gc/GG': hs(-8.139, -21.784), 'gg/GC': hs(-5.149, -12.508), 'gt/GA': hs(-8.991, -27.311),
'ag/TG': hs(-4.980, -15.426), 'cg/GG': hs(-4.441, -12.158), 'gg/CG': hs(-13.505, -36.021), 'tg/AG': hs(-2.775, -9.286),
#+T·T Mismatch
'ta/TT': hs(-3.744, -12.149), 'tc/TG': hs(-4.387, -13.520), 'tg/TC': hs(-6.346, -16.629), 'tt/TA': hs(-7.697, -25.049),
'at/TT': hs(-4.207, -14.307), 'ct/GT': hs(-8.176, -22.962), 'gt/CT': hs(-7.241, -20.622), 'tt/AT': hs(-2.051, -7.055),
#+A·C Mismatch
'aa/CT': hs(-1.362, -5.551), 'ac/CG': hs(-1.759, -6.511), 'ag/CC': hs(-6.549, -18.073), 'at/CA': hs(-3.563, -14.105),
'aa/TC': hs(-2.078, -10.088), 'ca/GC': hs(-5.868, -16.952), 'ga/CC': hs(-8.477, -24.565), 'ta/AC': hs(2.690, 4.965),
#+C·A Mismatch
'ca/AT': hs(-9.844, -29.673), 'cc/AG': hs(-3.761, -11.204), 'cg/AC': hs(-9.845, -27.316), 'ct/AA': hs(-3.389, -12.517),
'ac/TA': hs(0.753, -0.503), 'cc/GA': hs(-12.714, -35.555), 'gc/CA': hs(-12.658, -35.729), 'tc/AA': hs(-1.719, -7.023),
#+A·G Mismatch
'aa/GT': hs(2.193, 4.374), 'ac/GG': hs(-8.453, -22.672), 'ag/GC': hs(-1.164, -2.532), 'at/GA': hs(-7.418, -24.066),
'aa/TG': hs(-1.963, -9.013), 'ca/GG': hs(-8.712, -23.779), 'ga/CG': hs(-7.875, -21.661), 'ta/AG': hs(3.207, 7.156),
#+G·A Mismatch
'ga/AT': hs(-2.914, -9.402), 'gc/AG': hs(-9.131, -25.347), 'gg/AC': hs(-2.154, -3.871), 'gt/AA': hs(-8.515, -26.313),
'ag/TA': hs(-6.691, -21.148), 'cg/GA': hs(-3.960, -10.588), 'gg/CA': hs(-12.898, -34.656), 'tg/AA': hs(0.334, -0.440),
#+C·T Mismatch
'ca/TT': hs(0.382, -0.579), 'cc/TG': hs(-2.716, -8.000), 'cg/TC': hs(-10.363, -29.315), 'ct/TA': hs(-5.783, -20.173),
'ac/TT': hs(-0.692, -5.278), 'cc/GT': hs(-10.299, -28.503), 'gc/CT': hs(-9.062, -26.356), 'tc/AT': hs(2.073, 3.968),
#+T·C Mismatch
'ta/CT': hs(-5.485, -17.347), 'tc/CG': hs(1.451, 1.556), 'tg/CC': hs(-7.213, -20.128), 'tt/CA': hs(-2.397, -11.371),
'at/TC': hs(-0.633, -5.801), 'ct/GC': hs(-6.868, -21.000), 'gt/CC': hs(-5.853, -16.643), 'tt/AC': hs(0.211, -1.446),
#+G·T Mismatch
'ga/TT': hs(-5.551, -15.398), 'gc/TG': hs(-14.943, -40.148), 'gg/TC': hs(-8.110, -18.349), 'gt/TA': hs(-14.213, -40.041),
'ag/TT': hs(-7.130, -20.786), 'cg/GT': hs(-14.862, -39.430), 'gg/CT': hs(-14.622, -37.510), 'tg/AT': hs(-6.703, -18.111),
#+T·G Mismatch
'ta/GT': hs(-4.612, -14.039), 'tc/GG': hs(-9.798, -26.406), 'tg/GC': hs(-4.519, -11.065), 'tt/GA': hs(-4.523, -15.693),
'at/TG': hs(-2.364, -8.834), 'ct/GG': hs(-11.396, -30.732), 'gt/CG': hs(-6.233, -15.933), 'tt/AG': hs(-2.960, -9.305),
# LNA on top only
'aA/TT': hs(-7.193, -19.723), 'aC/TG': hs(-7.269, -18.336), 'aG/TC': hs(-7.536, -18.387), 'aT/TA': hs(-4.918, -12.943),
'cA/GT': hs(-7.451, -18.380), 'cC/GG': hs(-5.904, -11.904), 'cG/GC': hs(-9.815, -23.491), 'cT/GA': hs(-7.092, -16.825),
'gA/CT': hs(-5.038, -11.656), 'gC/CG': hs(-10.160, -24.651), 'gG/CC': hs(-10.844, -26.580), 'gT/CA': hs(-8.612, -22.327),
'tA/AT': hs(-7.246, -19.738), 'tC/AG': hs(-6.307, -15.515), 'tG/AC': hs(-10.040, -25.744), 'tT/AA': hs(-6.372, -16.902),
# LNA on bottom only
'Aa/TT': hs(-6.908, -18.135), 'Ac/TG': hs(-5.510, -11.824), 'Ag/TC': hs(-9.000, -22.826), 'At/TA': hs(-5.384, -13.537),
'Ca/GT': hs(-7.142, -18.333), 'Cc/GG': hs(-5.937, -12.335), 'Cg/GC': hs(-10.876, -27.918), 'Ct/GA': hs(-9.471, -25.070),
'Ga/CT': hs(-7.756, -19.302), 'Gc/CG': hs(-10.725, -25.511), 'Gg/CC': hs(-8.943, -20.833), 'Gt/CA': hs(-9.035, -22.742),
'Ta/AT': hs(-5.609, -16.019), 'Tc/AG': hs(-7.591, -19.031), 'Tg/AC': hs(-6.335, -15.537), 'Tt/AA': hs(-5.574, -14.149)
}

counter_char = {
    "A" : "T", "T" : "A", "C" : "G", "G" : "C", "a" : "t", "t" : "a", "c" : "g", "g" : "c"
}


def ow_salt_correction(divalent_mM, dntp_mM, K_mM, gc_content, seq_len):
   # depending on the value of div_monov_ratio respect to value of crossover_point 
   # Eq 16 (divalent corr, Owczarzy et al., 2008) or 
   # Eq 22 (monovalent corr, Owczarzy et al., 2004) should be used
   crossover_point = .22 
   if(dntp_mM >= divalent_mM):
      free_divalent = 0.00000000001
   else:
      free_divalent = (divalent_mM - dntp_mM)/1000.0
   a = b = c = d = f = g = 0
   if K_mM==0:
      div_monov_ratio = 6.0
   else: 
      div_monov_ratio = (math.sqrt(free_divalent))/(K_mM/1000)
   if div_monov_ratio < crossover_point:
      # use only monovalent salt correction, Eq 22 (Owczarzy et al., 2004) */
      return (((4.29 * gc_content) - 3.95) * 10**-5 * math.log(K_mM / 1000.0)) + (9.40e-6 * math.log(K_mM / 1000.0)**2)
   else:
       # magnesium effects are dominant, Eq 16 (Owczarzy et al., 2008) is used */
      b =- 9.11e-6
      c = 6.26e-5
      e =- 4.82e-4
      f = 5.25e-4
      a = 3.92e-5
      d = 1.42e-5
      g = 8.31e-5
      if(div_monov_ratio < 6.0):
         # in particular ratio of conc of monov and div cations some parameters of Eq 16 must be corrected (a,d,g) 
         a = 3.92e-5 * (0.843 - 0.352 * math.sqrt(K_mM/1000.0) * math.log(K_mM/1000.0))
         d = 1.42e-5 * (1.279 - 4.03e-3 * math.log(K_mM/1000.0) - 8.03e-3 * math.log(K_mM/1000.0)**2)
         g = 8.31e-5 * (0.486 - 0.258 * math.log(K_mM/1000.0) + 5.25e-3 * math.log(K_mM/1000.0)**3)

      return a + (b * math.log(free_divalent)) + gc_content * (c + (d * math.log(free_divalent))) + (1/(2 * (seq_len - 1))) * (e + (f * math.log(free_divalent)) + g * math.log(free_divalent)**2)

def divalent_to_monovalent(divalent_mM, dntp_mM):
    """Convert divalent cation concentration to monovalent equivalent (mM)."""
    if divalent_mM == 0:
        dntp_mM = 0
    if divalent_mM < dntp_mM:
        divalent_mM = dntp_mM  # dNTPs chelate Mg2+
    return 120 * math.sqrt(divalent_mM - dntp_mM)

def santa_salt_correction(K_mM, divalent_mM, dntp_mM, seq_len, S):
    K_mM = K_mM + divalent_to_monovalent(divalent_mM, dntp_mM)
    return S + 0.368 * (seq_len - 1) * math.log(K_mM / 1000.0 )

def parse_seq(seq):
    chars = []
    it = iter(seq)
    for char in it:
        if char == "+":
            char = next(it).lower()
        chars.append(char)
    return "".join(chars)

def check_sym(seq, comp_seq):
    for i in range(len(seq)):
        if counter_char[seq[i]] != comp_seq[-(i+1)]:
            return False
    else:
        return True

def calc_tm_with_lna(seq, comp_seq, dna_conc_nM=200, K_mM=50, divalent_mM=3, dntp_mM=.8,
                dmso_conc=0, dmso_fact=0, formamide_conc=0, salt="ow"):
    seq = parse_seq(seq)
    comp_seq = parse_seq(comp_seq)
    seq_len = len(seq)
    correction = 0
    sym = check_sym(seq, comp_seq)
    gc_content = calc_gc(seq)
    if salt.lower() != "ow":
        salt = "santa"

        
    H = 0
    S = 0

    for i in range(len(seq) - 1):
        window = seq[i:i+2] + '/' + comp_seq[-i-1:-i-3:-1]
        try:
            H += dublets[window].H
            S += dublets[window].S
        except:
           print("there was no match for you DNA doublet in the dictionary. Likely a 'N' or something snuck in")

    for end in [0, -1]:
        if seq[end] in ["A", "T", "a", "t"]:
            H += 2.3
            S += 4.1
        else:
            H += .1
            S += -2.8

    if salt == "santa":
        S = santa_salt_correction(K_mM, divalent_mM, dntp_mM, seq_len, S)
    elif salt =="ow":
        correction = ow_salt_correction(divalent_mM, dntp_mM, K_mM, gc_content, seq_len)


    # Compute Tm 
    Tm = (H*1000) / (S + 1.9872 * math.log(dna_conc_nM / (4e9 if sym == False else 1e9)))
    if salt == "ow":
        Tm = 1/((1/Tm)+correction) - 273.15
    else:
        Tm -= 273.15
    # DMSO and formamide corrections
    Tm -= dmso_conc * dmso_fact
    Tm += (0.453 * gc_content - 2.88) * formamide_conc

    return Tm

print(calc_tm_with_lna("ATACGGGC+T+C+GATAGCTAC", "GTAGCTATCGAGCCCGTAT"))

# only related to primer3py
# dmso_conc     =
# dmso_fact     = 
# formamide_conc=

# dna_conc_nM = dna_conc = Oligo Conc
# K_mM = mv_conc = Na+ Conc
# divalent_mM = dv_conc = Mg++ Conc
# dntp_mM = dntp_conc = dNTPs Conc


# ## References

# # 1. Owczarzy, R., You, Y., Groth, C. L., et al. (2011). "Stability and Mismatch Discrimination of Locked Nucleic Acid−DNA Duplexes." *Biochemistry* 50(40), 9988-10009.
# # 2. SantaLucia, J. Jr. (1998). "A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics." *Proc. Natl. Acad. Sci. U.S.A.* 95(4), 1460-1465.  
# from primer3-py scr/libprimer3/oligo.c line 371


