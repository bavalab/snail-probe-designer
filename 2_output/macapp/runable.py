#!/usr/bin/env python
# runable script for non-coding usage

import sys
import os

sys.path.append(os.path.abspath("~/Documents/\"Coding Projects\"/STARMAp_Probe_Designer/src/"))
from snail_probe_designer import snail_probe_designer

spd = snail_probe_designer() # uses the default settings for gc, tm, size, and separation
file_or_sequence = raw_input("Input the sequence or the path to a fasta file, then press Enter.\n")
gene_name = raw_input("Input the name of the gene, then press Enter.\n")
spd.prime(file_or_sequence, gene_name)
spd.get_kmers()
spd.score_kmers()
print "{} probe pairs found!".format(len(spd.probe_pairs))
output_filename = raw_input("Input the path/filename to save your output, then press Enter.\n")
spd.write_probes_to_csv(output_filename)
print "File has been saved. Please check the output folder you indicated."