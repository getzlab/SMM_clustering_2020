"""
SMM study features.
-----------------------
This file includes features and their classifications
for the SMM study. They are broken into feature type
for easy featurization for downstream clustering with annotations to
distinguish them.

"""


# --------------------------------------
# Numerical Entries
#
# SNVs, CNVs, Translocations
# --------------------------------------
snvs = [
    'KRAS',   #
    'NRAS',   #
    'BRAF',   #
    'PTPN11', # Phosphatase; tumor suppressor; should de-phosphorylate RAS
    'NF1',    # Tumor suppressor gene; controls the RAS genes
    'TP53',   #
    'ATM',    #
    'DIS3',   # very myeloma specific (see FAM46C); Found on 13q --> higher m-spke
    'FAM46C', # very myeloma specific; tumor supressor; exonuclease, ribosomasl proteins that cleaves protein complex
              # at certain time; mutations here are like loss of function related to the protein processing
              # altered protein regulation; myeloma is very sensitive to this because it produces a lot of immunoglobulins
              # Found on 1P
    'XBP1',   # can delete
    'ZNF292', # tumor suppressor; TF repressor
    'CCND1',  # CD-1 gene; activation; occurs more frequently with t(11;14)
    'RB1',    # Cell-cycle gene, tumor suppressor; 13q
    'CDKN1A',
    'CDKN1B',
    'CDKN2A', # Cyclin; chr9; loss of function mutation
    'CDKN2C', # Cyclin; 1p; loss of function mutation
    'LTB',    # Lymphotoxin; chr6; NFkB pathway
    'NFKBIA', # NFkB
    'NFKB2',  # ''
    'TRAF2',  # ''
    'TRAF3',  # ''
    'CYLD',   # ''; chr16
    'IRF4',   # Unique to myeloma; many mutations in the UTR --> they don't know what affects 3' or 5' and if it will
              # affect protein
    'PRDM1',  # Transcription factor for plasma cell maturation
    'SP140',  # Tumor necrosis factor / NFkB
    'MAX',    # Myc activator; tumor suppressor through a separate pathway
    'EGR1',
    'ACTG1',
    'HIST1H1B',
    'HIST1H1E',
    'ARID1A',
    'EP300',
    'DUSP2',
    'MAML2',
    'DTX1',
    'UBR5',  # Not sure of the impact but found commonly mutated
    'KLHL6', # Not sure of the impact but found commonly mutated
    'MET',
    'FGFR3', # t(4;14)
    'MAF',   # t(14;16), t(14;20)
    'MAFB',  # t(14;16), t(14;20)
    'NCOR1',
    'BTG1',
    'SAMHD1',
    'KMT2B',
    'KMT2C',
    'CREBBP',
    'SETD2',
    'SF3B1',
    'MAN2C1',
    'TCL1A',
    'PIM1',
    'PRKD2'  # Related to protein processing;
]

filt_snvs = [
    'KRAS',   #
    'NRAS',   #
    'BRAF',   #
    'PTPN11', # Phosphatase; tumor suppressor; should de-phosphorylate RAS
    'NF1',    # Tumor suppressor gene; controls the RAS genes
    'TP53',   #
    'ATM',    #
    'DIS3',   # very myeloma specific (see FAM46C); Found on 13q --> higher m-spke
    'FAM46C', # very myeloma specific; tumor supressor; exonuclease, ribosomasl proteins that cleaves protein complex
              # at certain time; mutations here are like loss of function related to the protein processing
              # altered protein regulation; myeloma is very sensitive to this because it produces a lot of immunoglobulins
              # Found on 1P
    'ZNF292', # tumor suppressor; TF repressor
    'CCND1',  # CD-1 gene; activation; occurs more frequently with t(11;14)
    'RB1',    # Cell-cycle gene, tumor suppressor; 13q
    'CDKN2A', # Cyclin; chr9; loss of function mutation
    'CDKN2C', # Cyclin; 1p; loss of function mutation
    'LTB',    # Lymphotoxin; chr6; NFkB pathway
    'NFKBIA', # NFkB
    'NFKB2',  # ''
    'TRAF2',  # ''
    'TRAF3',  # ''
    'CYLD',   # ''; chr16
    'IRF4',   # Unique to myeloma; many mutations in the UTR --> they don't know what affects 3' or 5' and if it will
              # affect protein
    'PRDM1',  # Transcription factor for plasma cell maturation
    'SP140',  # Tumor necrosis factor / NFkB
    'MAX',    # Myc activator; tumor suppressor through a separate pathway
    'UBR5',  # Not sure of the impact but found commonly mutated
    'KLHL6', # Not sure of the impact but found commonly mutated
    'FGFR3', # t(4;14)
    'MAF',   # t(14;16), t(14;20)
    'MAFB',  # t(14;16), t(14;20)
    'PRKD2'  # Related to protein processing;
]

deletions = [
    '6q_del',
    '17p_del',
    '11q_del',
    '16q_del',
    '12p_del',
    '14q_del',
    'del_4q',
    'del_21',
    '1p_del',
    'del_8p',
    'del_22q',
    'del_20q',
    'del_10p',
    'del_18q',
    '13q del'
]

amplifications = ['8q24 amp', 'amp_11_isolated', 'amp_9_isolated', '2p_amp', '1q_gain']

translocations =[
    't(11;14)',
    't(14;16)',
    't(14;20)',
    't(4;14)',
    't(6;14)'
]

hyperdiploidy = ['HRD']

trisomies = [
 'Tri_2',
 'Tri_3',
 'Tri_4',
 'Tri_5',
 'Tri_6',
 'Tri_7',
 'Tri_9',
 'Tri_10',
 'Tri_11',
 'Tri_12',
 'Tri_15',
 'Tri_17',
 'Tri_18',
 'Tri_19',
 'Tri_20',
 'Tri_21'
 ]

# --------------------------------------
# Metadata & Pathway Annotation
# --------------------------------------
pathways = [
    'Cell Cycle pathway (CCND1, RB1, CDKN1B, CDKN2A)',
    'NFkB pathway (NFKBIA, TRAF3, CYLD, NFKB2, and LTB)',
    'Protein processing and folding pathway (DIS3, FAM46C, 1p del)',
    'Other pathway (SP140,KDM6A,MET,MAX, ALK, IDH1, NOTHCH1)',
    'B cell_differentiation_pathway',
    'MAPK Pathway (KRAS, NRAS, BRAF, PTPN11)',
    # Sort of a pathway
    'DNA repair (TP53, ATM, 17p del, 11q del)'
]

biallelic = [
    'TRAF3_biallelic',
    'ZNF292_biallelic',
    'CYLD_biallelic',
    'DIS3_biallelic_alteration',
    'TP53_biallelic_inactivation',
    'FAM46C_biallelic',
    'Rb1_biallelic_inactivation'
]

myc = ['MYC _aberrations_all', 'MYC_translocation']

patient_metadata = [
    'AGE',
    'SEX',
    'RACE'
]

clinical_metadata = [
    'Albumin',
    'B2_microglobulin',
    'BM Involvement_per_ Bx',
    'Ca+2',
    'Creatinine',
    'Globulin',
    'Hgb',
    'IgH_type',
    'Inv/Uninv LC_ratio',
    'KLC',
    'LC_MM',
    'LDH',
    'LLC',
    'Light_chain_type',
    'M_spike',
    'Management_protocol',
    'Other_IgH_translocations',
    'Sequencing',
    'Total_protein',
    'serum_IgH_ level',
    # Irrelevant metadata
    'SNV',
    'SNV.1',
    'Translocations',
    'CNV'
]
