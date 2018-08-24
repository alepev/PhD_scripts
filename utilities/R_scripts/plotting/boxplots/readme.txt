to generate input data for the boxplot scripts (.tmp files), 
run codon_usage_per-gene_FULL_STATISTICS_PDF.py or codon_usage_per-gene_GC+LEN_only_PDF.py
scripts by providing one list at the time (instead of 2);
this way, besides performing normality checks, scripts do not delete .tmp files.
(behavior can be changed any moment by restoring .tmp removal)

.tmp are just text files containing a list of values, one per row (can be easily made from whatever source)

WHAT IS EACH STATISTICS SCRIPT FOR?

../01a_scripts_length+GC:
	> codon_usage_per-gene_GC+LEN_only_PDF.py 			-> GC+LEN for full transcript, 5'UTR, 3'UTR, CDS
																		-> with one list, can be used to assess distributions (+ create .tmp file)
																		-> with 2 lists, 1st is background and 2nd is considered sample from same population
	> codon_usage_per-gene_GC+LEN_only_PDF_2lists.py
																		-> can only be used with 2 lists, which are considered samples from same bkg population

../01b_scripts_length+GC+codon_bias_CDSonly:
	> codon_usage_per-gene_FULL_STATISTICS_PDF.py		-> GC,LEN,GC3,Nc for CDS ONLY (more filtering included, see scripts for details)
																			(when computing Nc, set MIN=100 for meaningful results)
																		-> with one list, can be used to assess distributions (+ create .tmp file)
																		-> with 2 lists, 1st is background and 2nd is considered sample from same population
	> codon_usage_per-gene_FULL_STATISTICS_PDF_2lists.py
																		-> can only be used with 2 lists, which are considered samples from same bkg population
	
lists are actual lists of gene identifiers (1 per row), 
other necessary files are listed in the script and should be present in the same directory (together with R scripts)