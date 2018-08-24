#sequence from Kern&Hupalo 2013, structure predicted by RNAfold

RNAfold < rna_seq.txt  > rna_RNAfold.out
cat rna_RNAfold.out
>araTha9_evofold_vf_1_16548 range=chr1:12451826-12451841 5'pad=0 3'pad=0 strand=+ repeatMasking=none
AUUUGAAGAGGAAGAU
................ (  0.00)
>araTha9_evofold_vf_1_49268 range=chr1:19660798-19660813 5'pad=0 3'pad=0 strand=+ repeatMasking=none
GUAUCUUGUUGGAUGC
((((((....)))))) ( -4.40)
>araTha9_evofold_vf_1_131750 range=chr1:29491193-29491202 5'pad=0 3'pad=0 strand=+ repeatMasking=none
GGCCACGGUU
.......... (  0.00)
>araTha9_evofold_vf_1_21575 range=chr1:7219-7228 5'pad=0 3'pad=0 strand=+ repeatMasking=none
GGGAAUUCCC
(((....))) ( -2.50)

RNAplot < rna_RNAfold.out
RNAdistance -Xm < rna_RNAfold.out  > rna_RNAdist.out
> f   4
24 
6 18 
18 6 12 


#both sequence and structure from Kern&Hupalo 2013
cat Kern\&Hupalo_RNAs_seq+str.txt 
>araTha9_evofold_vf_1_16548 range=chr1:12451826-12451841 5'pad=0 3'pad=0 strand=+ repeatMasking=none
ATTTGAAGAGGAAGAT
((((........))))
>araTha9_evofold_vf_1_49268 range=chr1:19660798-19660813 5'pad=0 3'pad=0 strand=+ repeatMasking=none
GTATCTTGTTGGATGC
(((((.....)).)))
>araTha9_evofold_vf_1_131750 range=chr1:29491193-29491202 5'pad=0 3'pad=0 strand=+ repeatMasking=none
GGCCACGGTT
(((....)))
>araTha9_evofold_vf_1_21575 range=chr1:7219-7228 5'pad=0 3'pad=0 strand=+ repeatMasking=none
GGGAATTCCC
((......))

RNAplot < Kern&Hupalo_RNAs_seq+str.txt
RNAdistance -Xm < Kern&Hupalo_RNAs_seq+str.txt > rnas_RNAdist.out
> f   4
6 
6 6 
6 6 4 

