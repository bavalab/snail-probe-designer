# Oligo PLA Designer Class
# Author: Eric Cramer <eric.cramer@curie.fr>

import os#, primer3

class oligo_pla_designer:
    """

    Oligonucleotide Proximity Ligation Assay Designer (OPLAD)

    Defines an object for designing oligonucleotide probes for a given DNA/RNA sequence. 
    Sequences can be from fasta files or strings.
    Probes are designed in pairs with 5' and 3' modifications provided by the user.
    For example, in a SNAIL design, one padlock and one splint will be created with the padlock detection 
    readout/barcode, padlock phosphorylated end, and splint arm provided by the user.

    The OPLAD class uses python wrapper functions for the primer3 library to calculate thermodynamic stability 
    of the probe products.

    OPLAD is capable of batch processing

    """

    def __init__(self, size=(18, 22), tm=(50.0, 65.0), gc=(30.0, 70.0), sep=(0, 4)):
        """
        Class Constructor
        Declares the object and sets its internal values

        Input:
        self = the object being declared, standard python syntax
        size = the length of the probe in base pairs. Default is a tuple of (18, 22) bp 
        tm = the desired range of melting temperatures of the probes in C. Default is a tuple, (50.0, 65.0)
        gc = the "G" and "C" content of the probes as a percent. Default is a tuple (40.0, 60.0)
        sep = the desired separation between the two probes in base pairs. Default is a tuple of (4)
        padlock_barcode = the barcode to be attached to the padlock probe. Check code for default value
        padlock_leader = the phosphoylated end of the padlock probe that binds the splint. check code for default value
        splint_connector = the arm of the splint which binds the padlock probe. Check code for default value

        Usage:
        spd = snail_probe_designer()
        """
        self.probe_size = size 
        self.probe_tm = tm
        self.probe_gc = gc 
        self.probe_sep = sep
        self._sequences = {} # empty dictionary to store the sequences
        self._top_probes = None # init a variable to store the best probes

    def set_modifications(self, padlock_5prime, padlock_3prime, splint_5prime, splint_3prime):
        """
        Sets the modification tags attached to each end of the probe target sequences.

        Input:
        padlock_5prime = the modification for the 5' (leading) end of the padlock
        padlock_3prime = the modification for the 3' (trailing) end of the padlock
        splint_5prime = the modification for the 5' (leading) end of the splint
        splint_3prime = the modification for the 3' (trailing) end of the splint

        Usage:
        spd.set_modifications(padlock_5prime='NNN', padlock_3prime='NNN', splint_5prime='NNN', splint_3prime='NNN')

        """
        self._padlock_5prime = padlock_5prime
        self._padlock_3prime = padlock_3prime
        self._splint_5prime = splint_5prime
        self._splint_3prime = splint_3prime

    def set_sequence(self, seq_or_file_or_folder, gene_name=None):
        """
        Prime the probe designer with a sequence or folder of fasta files to find probes for

        Input:
        seq_or_file_or_folder = either the sequence as a string, a fasta file containing sequence(s), or a folder of fasta files

        Usage:
        spd.set_sequence('fasta_file')
        spd.set_sequence("SEQUENCE")
        spd.set_sequence('directory_of_fasta_files')

        """
        # first check to see if the user passes a directory
        if os.path.isdir(seq_or_file_or_folder):
            # if so process all of the files in the directory
            for file in os.listdir(seq_or_file_or_folder):
                if file.endswith('fasta'):
                    with open(seq_or_file_or_folder + file, 'r') as  fasta_file:
                        for line in fasta_file:
                            line = line.strip()
                            if not line: continue
                            if line.startswith('>'):
                                seq_name = line[1:-1] # remove the leading '>'
                                if seq_name not in self._sequences:
                                    self._sequences[seq_name] = ""
                                continue
                            self._sequences[seq_name] = line
        # if the user passed a fasta file
        if seq_or_file_or_folder.endswith("fasta"):
            # open and import the fasta file
            with open(seq_or_file_or_folder, 'r') as fasta_file:
                for line in fasta_file:
                    line = line.strip()
                    if not line: continue
                    if line.startswith('>'):
                        seq_name = line[1:-1] # remove the leading '>'
                        if seq_name not in self._sequences:
                            self._sequences[seq_name] = ""
                        continue
                    self._sequences[seq_name] = line
        # otherwise treat the input as a sequence
        else:
            if gene_name is not None:
                self._sequences[gene_name] = seq_or_file_or_folder.strip()
            else:
                print("Please provide a gene name for the sequence.")

def main():
	"""
	Currently used for testing. Otherwise will provide CLI support in future versions.
	"""
    oplad = oligo_pla_designer()
    print("From directory")
    oplad.set_sequence("../../0_data/raw/")
    print(oplad._sequences)
    oplad = oligo_pla_designer()
    print("From file")
    oplad.set_sequence("../../0_data/raw/actb.fasta")
    print(oplad._sequences)
    print("From sequence")
    oplad = oligo_pla_designer()
    oplad.set_sequence("ATGGATGACGATATCGCTGCGCTGGTCGTCGACAACGGCTCCGGCATGTGCAAAGCCGGCTTCGCGGGCGACGATGCTCCCCGGGCTGTATTCCCCTCCATCGTGGGCCGCCCTAGGCACCAGGGTGTGATGGTGGGAATGGGTCAGAAGGACTCCTATGTGGGTGACGAGGCCCAGAGCAAGAGAGGTATCCTGACCCTGAAGTACCCCATTGAACATGGCATTGTTACCAACTGGGACGACATGGAGAAGATCTGGCACCACACCTTCTACAATGAGCTGCGTGTGGCCCCTGAGGAGCACCCTGTGCTGCTCACCGAGGCCCCCCTGAACCCTAAGGCCAACCGTGAAAAGATGACCCAGATCATGTTTGAGACCTTCAACACCCCAGCCATGTACGTAGCCATCCAGGCTGTGCTGTCCCTGTATGCCTCTGGTCGTACCACAGGCATTGTGATGGACTCCGGAGACGGGGTCACCCACACTGTGCCCATCTACGAGGGCTATGCTCTCCCTCACGCCATCCTGCGTCTGGACCTGGCTGGCCGGGACCTGACAGACTACCTCATGAAGATCCTGACCGAGCGTGGCTACAGCTTCACCACCACAGCTGAGAGGGAAATCGTGCGTGACATCAAAGAGAAGCTGTGCTATGTTGCTCTAGACTTCGAGCAGGAGATGGCCACTGCCGCATCCTCTTCCTCCCTGGAGAAGAGCTATGAGCTGCCTGACGGCCAGGTCATCACTATTGGCAACGAGCGGTTCCGATGCCCTGAGGCTCTTTTCCAGCCTTCCTTCTTGGGTATGGAATCCTGTGGCATCCATGAAACTACATTCAATTCCATCATGAAGTGTGACGTTGACATCCGTAAAGACCTCTATGCCAACACAGTGCTGTCTGGTGGTACCACCATGTACCCAGGCATTGCTGACAGGATGCAGAAGGAGATTACTGCTCTGGCTCCTAGCACCATGAAGATCAAGATCATTGCTCCTCCTGAGCGCAAGTACTCTGTGTGGATCGGTGGCTCCATCCTGGCCTCACTGTCCACCTTCCAGCAGATGTGGATCAGCAAGCAGGAGTACGATGAGTCCGGCCCCTCCATCGTGCACCGCAAGTGCTTCTAG", gene_name="ACTB")
    print(oplad._sequences)

if __name__ == "__main__":
    main()