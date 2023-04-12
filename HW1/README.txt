Compile and run with python 3.
Recieves only one argument that is the filepath of the fasta containing two sequences to align.
EX: python3 trana_align.py filepath.fasta

The traceback function will check for left gaps first, then up gaps, and if neither could result in the alignment table, then assumes there was a substitution.
The above matters only for breaking ties. Different best alignments will be generated depending on how the traceback function breaks ties.