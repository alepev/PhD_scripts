# PhD_scripts
A collection of Python&R scripts from my PhD research (mostly plant phylogenetics)

## Intro

This is a partial collection of scripts I made during my PhD,
usable with publicly available datasets and software for (plant) genomics & phylogenetics.
Some of the original datasets consisted of huge genomic files and annotations,
so instead of including them directly in the repository, 
I opted for a reference list with download links whenever possible.

Had I known about version control earlier I would have structured things differently,
likely into separate repositories for each project, 
and spent more time in creating mock datasets for testing
(this is something I started doing more recently, see for instance test data in
01_FASTA_operations/test_files_for_fasta_EXTRACT-EXTEND).

..Better late than never I guess ;)

## Contents

**Folders 01-03** and **Utilities** contain scripts I used for general cross-project operations,
such as handling of FASTA files and gene functional annotation.

**Folders 04-05** are specific to two publications, respectively

  *Bai, Peviani, van der Horst, et al. (2017). 
  Extensive translational regulation during seed germination revealed by polysomal profiling. 
  The New Phytologist, 214(1), 233–244.*

and

  *Peviani, Lastdrager, Hanson, & Snel (2016). 
  The phylogeny of C/S1 bZIP transcription factors reveals a shared algal ancestry 
  and the pre-angiosperm translational regulation of S1 transcripts. 
  Scientific Reports, 6, 30444.*

**Folder 05** in particular is an actual pipeline (with steps described in the HOWTO file)
I used to resolve a complex gene family phylogeny, combining several different tools 
and manual inspection of the resulting alignments iteratively
(reason why it would be extremely difficult to fully automate).
