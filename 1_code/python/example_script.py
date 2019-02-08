# Example usage of the SNAIL Probe Designer

from snail_probe_designer import snail_probe_designer

spd = snail_probe_designer() # uses the default settings for gc, tm, size, and separation
spd.prime("../data/mll-sequence.fasta", "MLL")
spd.get_kmers()
spd.score_kmers()
spd.write_probes_to_csv("../tmp/test_mll_probes.csv")
print "{} probe pairs found!".format(len(spd.probe_pairs))

