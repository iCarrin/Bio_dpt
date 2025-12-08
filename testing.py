from Bio.Seq import Seq
import primer3
import re
# import numpy as np
# import prc_lib as plib
# # snp = 'ATGCAATTGGCCAAATTTGGGCCCAAAATTTTGGGGCCCCAAAAATTTTTGGGGGCCCCCAAAAAATTTTTTGGGGGGCCCCCC'
print(primer3.bindings.calc_hairpin('AAAATAACTTCAATTCTACACAGTACGTGTC'))

# snp = 'ABCDEFGHIJKLmNOPQRSTUVWXYZ'
# pessimism = 5
# start = 0
# middle = 13
# max_dist = 12
# min_dist = 2
# # for i in range (3):
# #     for length in range(10):
# # trimmed = snp[middle+min_dist:middle+max_dist]
# # print(trimmed)
# # start += 1
# # print("try again")
#     # print(length)
# myString = "Cot"

# print (myString)
# myString = myString[:1] + "A" + myString[2:]
# print (myString)

# snp_data = [{'snpID': 'rs1799971', 'allele': 'A', 'sequence': 'TCCTGGGTCAACTTGTCCCACTTAGATGGCAACCTGTCCGACCCATGCGGTCCGAACCGCA', 'position': 30}, 
#             {'snpID': 'rs1799971', 'allele': 'G', 'sequence': 'TCCTGGGTCAACTTGTCCCACTTAGATGGCGACCTGTCCGACCCATGCGGTCCGAACCGCA', 'position': 30}, 
#             {'snpID': 'rs599839', 'allele': 'G', 'sequence': 'AAAAAAAGAGAAAGAAATAGGAGCAGGATCGACTTCCAGATATACAGAGAATATAAAAATA', 'position': 30}, 
#             {'snpID': 'rs599839', 'allele': 'A', 'sequence': 'AAAAAAAGAGAAAGAAATAGGAGCAGGATCAACTTCCAGATATACAGAGAATATAAAAATA', 'position': 30}, 
#             {'snpID': 'rs599839', 'allele': 'C', 'sequence': 'AAAAAAAGAGAAAGAAATAGGAGCAGGATCCACTTCCAGATATACAGAGAATATAAAAATA', 'position': 30}, 
#             {'snpID': 'rs599839', 'allele': 'T', 'sequence': 'AAAAAAAGAGAAAGAAATAGGAGCAGGATCTACTTCCAGATATACAGAGAATATAAAAATA', 'position': 30}]

# snp_short = snp[0:10]
# print(f"Normal string: {snp_short}")
# seq = Seq.Seq(snp_short)
# print(f"Now a sequence in BIO Python: {seq}")
# reverse = seq.reverse_complement()
# print(f"Now a reverse sequence in Bio Python: {reverse}")
# rev_string = str(reverse)
# print(f"Now a string again: {rev_string}")
# reverse

# np_frame = np.array(snp_data)

# for i in np_frame:
#     print(i['sequence'])
# plib.generate_allele_specific_primers(snp_data, 24, 26)
test = "TCCTGGGTCAACTTGTCCCACTTANATGGCYACCTGNCCGACCCATGCGGTCCGAACCGCA"
# print(test[30+7:])
# print(test[:30+1-6])# the end is exclusive 
print(test[-(6+2):-(0)])# the end is exclusive 
# reverse_test = str(Seq(test[30+7:]).reverse_complement())
# print(reverse_test)

# results =
# YACCTGNCCGACCCATGCGGTCCGAACCGCA
# TCCTGGGTCAACTTGTCCCACTTANATGGCY

# CCGACCCATGCGGTCCGAACCGCA
# TGCGGTTCGGACCGCATGGGTCGG
# inverted and spun around



