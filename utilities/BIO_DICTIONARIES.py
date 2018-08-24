#####################################################
#AMINOACID SYMBOLS
#####################################################
AA= {}
AA['Ala']='A'
AA['Arg']='R'
AA['Asn']='N'
AA['Asp']='D'
AA['Cys']='C'
AA['Glu']='E'
AA['Gln']='Q'
AA['Gly']='G'
AA['His']='H'
AA['Ile']='I'
AA['Leu']='L'
AA['Lys']='K'
AA['Met']='M'
AA['Phe']='F'
AA['Pro']='P'
AA['Ser']='S'
AA['Thr']='T'
AA['Trp']='W'
AA['Tyr']='Y'
AA['Val']='V'

#####################################################
#STANDARD TRANSLATION CODE
#####################################################
CODE= {}
CODE['ATG']='M' #translation start
CODE['TGG']='W'
for triplet in ['TTT','TTC']:
	CODE[triplet]='F'
for triplet in ['TAT','TAC']:
	CODE[triplet]='Y'
for triplet in ['TGT','TGC']:
	CODE[triplet]='C'
for triplet in ['AAT','AAC']:
	CODE[triplet]='N'
for triplet in ['CAA','CAG']:
	CODE[triplet]='Q'
for triplet in ['GAT','GAC']:
	CODE[triplet]='D'
for triplet in ['GAA','GAG']:
	CODE[triplet]='E'
for triplet in ['AAA','AAG']:
	CODE[triplet]='K'
for triplet in ['CAT','CAC']:
	CODE[triplet]='H'
for triplet in ['ATT','ATC','ATA']:
	CODE[triplet]='I'
for triplet in ['GTT','GTC','GTA','GTG']:
	CODE[triplet]='V'
for triplet in ['ACT','ACC','ACA','ACG']:
	CODE[triplet]='T'
for triplet in ['GCT','GCC','GCA','GCG']:
	CODE[triplet]='A'
for triplet in ['CCT','CCC','CCA','CCG']:
	CODE[triplet]='P'
for triplet in ['GGT','GGC','GGA','GGG']:
	CODE[triplet]='G'
for triplet in ['TTA','TTG','CTT','CTC','CTA','CTG']:
	CODE[triplet]='L'
for triplet in ['TCT','TCC','TCA','TCG','AGT','AGC']:
	CODE[triplet]='S'
for triplet in ['CGT','CGC','CGA','CGG','AGA','AGG']:
	CODE[triplet]='R'
for triplet in ['TAA','TAG','TGA']:
	CODE[triplet]='*' #translation stop
