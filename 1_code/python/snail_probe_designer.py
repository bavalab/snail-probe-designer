"""

SNAIL Probe Designer class and associated functions
Author: Eric Cramer <eric.cramer@curie.fr>

"""

# Import necessary dependencies
import math, os
from itertools import groupby

# Utility functions

def rev_comp(in_str, d = 0):
    """
    Reverse complement input strings, default is to convert between dna sequences

    Input:
    in_str = the string or sequence to be reverse complemented
    d = an integer indicating the type of conversion to do
        if d = 0 then converting dna to dna
        if d = 1 then converting rna to dna
        if d = 2 then converting dna to rna
        if d = 3 then converting rna to rna

    Output:
    A string which is translated into the complement, and reversed with respect to the input string

    Usage:
    spd.rev_comp(in_str)
    """
    if "U" in in_str and "T" in in_str: return "Invalid input. Is it RNA or DNA? Make up your mind."
    conversion = {"A":"T", "T":"A", "G":"C", "C":"G"}
    if d == 1:
        conversion = {"A":"T", "U":"A", "G":"C", "C":"G"} # rna to dna
    elif d == 2: 
        conversion = {"A":"U", "T":"A", "G":"C", "C":"G"} # dna to rna
    elif d == 3:
        conversion = {"A":"U", "U":"A", "G":"C", "C":"G"} # rna to rna
    return "".join(conversion.get(base, base) for base in reversed(in_str.upper()))

def calc_gc(p):
    """
    calculates the GC content of a primer p

    Input:
    p = the primer sequence (a string)

    Output:
    The GC content as a percentage (a float)

    Usage:
    calc_gc(primer_sequence)
    """
    return(100*((p.count("G") + p.count("C"))/float(len(p))))

def check_gc(primer_gc, probe_seq):
    """
    Checks the G and C content of the string to make sure it is within parameters for annealing

    Input:
    primer_gc = the GC content parameters for the primer as a tuple
    probe_seq = the sequence of a given probe

    Output:
    A boolean indicating if the GC content is within parameters

    Usage:
    gc_content = spd.check_gc(primer_gc, seq)

    """ 
    gc_content = calc_gc(probe_seq)
    if gc_content > primer_gc[0] and gc_content < primer_gc[1]:
        return(True)
    return(False)

def calc_dhds(p):
    """
    Calculates the delta H (entropy) and delta S (enthalpy) of a given primer/probe
    Uses the pre-set nearest neighbor dynamics tables given in Kampke et al (2000)
    Input:
    p = the probe sequence

    Output:
    A dictionary containing the entropy and enthalpy of the probe sequence

    Usage:
    dh_and_ds = spd.calc_dhds(probe_sequence)

    """
    dh, ds = 0, 0
    # init the lookup tables
    nnt_dh = {"AA":9.1, "TT":9.1, "AT":8.6, "TA":6.0, "CA":5.8, "TG":5.8, "GT":6.5, "AC":6.5, "CT":7.8, "AG":7.8, "GA":5.6, "TC":5.6, "CG":11.9, "GC":11.1, "GG":11.0, "CC":11.0}
    nnt_ds = {"AA":24.0, "TT":24.0, "AT":23.9, "TA":16.9, "CA":12.9, "TG":12.9, "GT":17.3, "AC":17.3, "CT":20.8, "AG":20.8, "GA":13.5, "TC":13.5, "CG":27.8, "GC":26.7, "GG":26.6, "CC":26.6}
    
    # calculate dh and ds
    for i in range(len(p)):
        if i+2 > len(p): break # don't walk off the end of the sequence
        doublet = p[i:i+2]
        dh += nnt_dh[doublet]
        ds += nnt_ds[doublet]
    return({"dH":dh, "dS":ds})

def calc_tm_energy(p):
    """
    Calculates the melting temperature of a primer/probe sequence
    Uses the forumla suggested by Freier et al (1986)

    Input:
    p = the probe sequence

    Output:
    a float indicating the melting temperature for the probe in C

    Usage:
    tm = spd.calc_tm_energy(p)
    """
    # init constants for the equation
    R = 1.987 # gas constant
    gamma = 50*math.pow(10, -9) # molar concentration of primer
    T0 = -273.15 # temp in C
    t = -21.6 # temp correction in C
    dh_ds = calc_dhds(p)
    tm = (dh_ds["dH"]/(dh_ds["dS"] + R + math.log(gamma/4))) + T0 + t
    return(tm)

def calc_tm_simple(p):
    """
    Calculates the melting temperature of a primer/probe sequence
    Uses the simplified formula proposed by Kampke et al (2000) for probes and primers 15 - 28 base pairs long

    Input: 
    p = the probe sequence

    Output:
    a float indicating the melting temperature for the probe in C

    Usage:
    tm = spd.calc_tm_simple(p)
    """
    return((4*p.count("G"))+(4*p.count("C"))+(2*p.count("A"))+(2*p.count("T"))
    )

def calc_alignment_score(p1, p2):
    """ 
    Calculates exact character alignment between 2 DNA sequences using sliding windows to align

    Input:
    p1 = the sequence as a string for the first probe
    p2 = the sequence as a string for the second probe

    Output:
    Returns the score of the alignment with 2 points given for every A-T alignment and 4 points given for every G-C alignment

    Usage:
    alignment_score = calc_alignment_score(seq1, seq2)
    """
    score = 0
    base_matches = {"A":"T", "T":"A", "C":"G", "G":"C"}
    for i in range(len(p1)):
        if base_matches[p1[i]] == p2[i]:
            if p1[i] == "A" or p1[i] == "T": score += 2
            else: score += 4
    return(score)

def get_max_alignment_score(p1, p2):
    """
    Checks the alignment between two sequences
    Use to find the two probes least likely to hybridize to themselves and each other
    To check self-hybridization, reverse the probe to simulate 5' to 3' matching and pass to function

    Input:
    p1 = the sequence as a string for the first probe
    p2 = the sequence as a string for the second probe

    Output:
    Returns the maximum score all alignments between p1 and p2 using the sliding windows. Calls calc_alignment_score for each window.

    Usage:
    max_alignment_score = get_max_alignment_score(seq1, seq2)
    """
    # check to make sure the strings are the same length, otherwise reject
    if len(p1) != len(p2): 
        print("The probes should be the same length.")
        return(False)
    # init list to store the scores for each alignment
    scores = []
    # use a sliding window to align one string to the other with an "overhang"
    sliding, stationary = "", ""
    # subscript one string from the beginning and the other from the end, then check how many characters match/align
    for i in range(1, 2*len(p1)):
        if i <= len(p1):
            sliding = p1[:i]
            stationary = p2[-i:]
        else:
            # we need to subtract the length of the string to slice it properly
            sliding = p1[i-len(p1):]
            stationary = p2[:-(i-len(p1))]
        scores.append(calc_alignment_score(sliding, stationary))
    return(max(scores))

# Definition of the snail_probe_designer class

class snail_probe_designer:
    """

    SNAIL Probe Designer Class

    Defines an object for designing SNAIL probes for a given DNA/RNA sequence. 
    Sequences can be from fasta files or strings.
    Probes are designed in pairs, one padlock and one splint.

    Author: Eric Cramer <eric.cramer@curie.fr>

    """

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
    def __init__(self, size=(18, 22), tm=(50.0, 65.0), gc=(30.0, 70.0), sep=(0, 4), padlock_barcode="AATTATTACTGAAACATACACTAAAGATA", padlock_leader="PACATTA", splint_connector="TAATGTTATCTT"):
        self.probe_size = size 
        self.probe_tm = tm 
        self.probe_gc = gc 
        self.probe_sep = sep 
        self.padlock_barcode = padlock_barcode
        self.padlock_leader = padlock_leader
        self.splint_connector = splint_connector
        self.top_probes = None

    def prime(self, seq_or_file, gene_name=None):
        """
        Prime the probe designer with a sequence to find probes for
        This function was designed with single genes in mind (e.g. one probe designer object per gene)

        Input:
        seq_or_file = either the sequence as a string or a fasta file containing a single sequence
        gene_name = the name of the gene as a string. Defaults to none

        Usage:
        spd.prime(fasta_file, gene_name)
        spd.prime("SEQUENCE", gene_name)

        """
        self.gene_name = gene_name
        if seq_or_file.endswith("fasta"):
            # open and import the fasta file
            fasta_file = open(seq_or_file)
            faiter = (x[1] for x in groupby(fasta_file, lambda line: line[0] == ">"))
            for header in faiter:
                self.fasta_header = next(header)[1:].strip() # removes the leading ">"
                self.sequence = "".join(s.strip() for s in next(faiter))
        else:
            self.sequence = seq_or_file

    def get_kmers(self):
        """
        Finds all pairs of probes for a given size, temperature, gc content, and separation combination
        
        Output:
        sets an internal list of probe pairs as a tuple

        Usage:
        spd = snail_probe_designer()
        spd.prime(filename, gene_name)
        spd.get_kmers()

        Pseudocode:
        ---
        init a list of probe pairs stored as tuples
        for each potential probe separation
            for each potential probe length
                for each character index i in the sequence
                    if i+2(primer_size)+primer_sep > sequence length
                        stop and return the list of probe pairs because we hit the end
                    p1 = substring from i to i+primer_size
                    p2 = substring from i+primer_sep+primer_size to i+primer_sep+(2*primer_size)
                    if gc of p1 or p2 > primer_gc
                        continue to next iteration
                    anneal_temp = calculate the probe annealing temperatures
                    if anneal_temp for p1 or p2 > primer_tm
                        continue to next iteration
                    add (p1,p2) to the list of probe pairs
            
        """
        self.probe_pairs = []
        for sep in range(self.probe_sep[0], self.probe_sep[1]):
            for p in range(min(self.probe_size), max(self.probe_size)+1):
                for i in range(len(self.sequence)):
                    if i+2*(p)+sep > len(self.sequence): 
                        break
                    # reverse complement the probes and check their separation
                    p1 = rev_comp(self.sequence[i:i+p])
                    p2 = rev_comp(self.sequence[i+sep+p:i+sep+(2*p)])
                    # is the GC content within parameters?
                    if not check_gc(self.probe_gc, p1) and not check_gc(self.probe_gc, p2):
                        continue
                    # calculate the primer temperatures
                    temp_p1 = calc_tm_simple(p1)
                    temp_p2 = calc_tm_simple(p2)
                    # are the two primers close enough in melting temperature?
                    if math.fabs(temp_p1 - temp_p2) >= 5.0:
                        continue
                    # are the primers within the desired melting temperature range?
                    if min(temp_p1, temp_p2) < self.probe_tm[0] or max(temp_p1, temp_p2) > self.probe_tm[1]:
                        continue
                    self.probe_pairs.append((p1, p2, math.fabs(temp_p1 - temp_p2)))

    def score_kmers(self):
        """
        Score the found kmers to find which are optimal
        Optimal probes do not self hybridize, or hybridize to the other probe in the pair
        
        Output:
        sets an internal list with the probe pairs, and a mixed score as a 3-part tuple
        e.g. (probe1, probe2, equally weighted average of: probe1 self-score, probe2 self-score, and inter-score)
        
        Usage:
        spd = snail_probe_designer()
        spd.prime(filename, gene_name)
        spd.get_kmers()
        spd.score_kmers()
        """
        self.scored_probe_pairs = []
        if len(self.probe_pairs) > 0:
            probe_pairs_scored = []
            for probe_pair in self.probe_pairs:
                self_score_0 = get_max_alignment_score(probe_pair[0], probe_pair[0][::-1])
                self_score_1 = get_max_alignment_score(probe_pair[1], probe_pair[1][::-1])
                inter_score = get_max_alignment_score(probe_pair[0], probe_pair[1][::-1])
                mixed_score = ((1/float(3))*self_score_0 + (1/float(3))*self_score_1 + (1/float(3))*inter_score)/float(3)
                probe_pairs_scored.append((probe_pair[0], probe_pair[1], probe_pair[2], mixed_score))
            self.scored_probe_pairs = sorted(probe_pairs_scored, key=lambda tup: (tup[2], tup[3]))

    def get_top_probes(self, top_x = 5):
        """
        Return the top X kmers (the X number of kmers with the lowest scores)

        Input:
        top_x = the desired number of probe pairs (e.g. top_x = 5 means 5 probe pairs or 10 total oligo probes). Default is 5.
        no_overlap = do not show overlapping probe pairs. Default is True
        
        Output:
        returns a list of the top_x non-overlapping probes (unless otherwise specified)

        Usage:
        spd = snail_probe_designer()
        spd.prime(filename, gene_name)
        spd.get_kmers()
        spd.score_kmers()
        spd.get_top_probess()
        """
        if len(self.scored_probe_pairs) > 0:
            # init a list of probes that are in the top x and meet additional criteria
            self.top_probes = [self.scored_probe_pairs[0]]

            # get the number of probes in the designer
            n_probes = len(self.scored_probe_pairs)

            # for each probe pair in the self contained list of pairs (until top_x is met)
            for i in range(1, n_probes):
                new_probe = self.scored_probe_pairs[i]
                overlap_flag = 0
                # check to see if it overlaps with a top probe pair
                for probe in self.top_probes:
                    if not self.check_overlaps(new_probe, probe):
                        overlap_flag += 1
                # if the new probe does not overlap with a top probe, add it
                if overlap_flag is len(self.top_probes):
                    self.top_probes.append(new_probe)
                # if the list is long enough break and return
                if len(self.top_probes) is top_x:
                    break
        else:
            print("No probe pairs found.")
        return(self.top_probes)

    def write_probes_to_csv(self, filename, num_probes=5):
        """
        Writes a given number of probes to a csv file in the order of their hybridization score

        The file is formatted as a csv table (for easy opening in Excel). Even number probes are padlocks and odd number probes are splints. Example:

        | Gene Name | Probe Number | Probe Name | Padlock Leader | Probe Sequence | Probe Connector or Barcode| Full Probe Sequence|
        |---|---|---|---|---|---|---|
        | ALB | 0 | ALB_0 | PACATTA    | GGGGGAGGTTTGGGTTGTCA    | AATTATTACTGAAACATACACTAAAGATA    | PACATTAGGGGGAGGTTTGGGTTGTCAAATTATTACTGAAACATACACTAAAGATA |
        | ALB | 1 | ALB_1 | | CATCAACCTCTGGTCTCACC | TAATGTTATCTT | CATCAACCTCTGGTCTCACCTAATGTTATCTT |
        ...

        Input:
        filename = the name and path of the file to write the probe information to as a string 
        num_probes = the number of probe pairs to write save

        Output:
        the file containing the probe information is written to the given filename and folder as a csv

        Usage:
        spd = snail_probe_designer()
        spd.prime(filename, gene_name)
        spd.get_kmers()
        spd.score_kmers()
        spd.write_probes_to_csv(path+filename)
        """
        structured_output = ""
        counter = 0
        probes_to_write = self.top_probes if self.top_probes is not None else self.get_top_probes(top_x = num_probes)
        for probe_pair in probes_to_write:
            structured_output += self.gene_name + "," + str(counter) + "," + self.gene_name + "_" + str(counter) + "," + self.padlock_leader + "," + probe_pair[0] + "," + self.padlock_barcode + "," + self.padlock_leader + probe_pair[0] + self.padlock_barcode + "\n"
            counter += 1
            structured_output += self.gene_name + "," + str(counter) + "," + self.gene_name + "_" + str(counter) + "," + "," + probe_pair[1] + "," + self.splint_connector + "," + probe_pair[1] + self.splint_connector + "\n"
            counter += 1

        outfile = open(filename, 'w')
        outfile.write(structured_output)
        outfile.close()

    def write_probes_to_eurogentec(self, filename, num_probes=5, product='Custom Oligos', syn_scale='40 nmol (10-99 bases)', pur='SePOP Desalted', form='Dried'):
        """
        Writes probes to an excel file in the order of their hybridization score, ready for upload to Eurogentec (https://secure.eurogentec.com/products/custom-oligonucleotides.html)
        ---
        Input:
        filename = the name and path of the file to write the probe information to as a string 
        num_probes = the number of probe pairs to write save. Default is 5
        product = the type of product being requested from Eurogentec. Default is 'Custom Oligos'
        syn_scale = the synthesis scale of the oligos. Default is '40 nmol (10-99 bases)'
        pur = the purity of the oligos. Default is 'SePOP Desalted'
        form = The form that the oligos will be sent in. Default is 'Dried'

        Output:
        the file containing the probe information is written to the given filename and folder as a .xlsx

        Usage:
        spd = snail_probe_designer()
        spd.prime(filename, gene_name)
        spd.get_kmers()
        spd.score_kmers()
        spd.write_probes_to_eurogentec(path+filename)
        """

        # check to see if the excel writer is installed
        try:
            import xlsxwriter
        except ImportError:
            print("Error, module not available. Cannot write to Excel file.")
            return(False)

        workbook = xlsxwriter.Workbook(filename)
        worksheet = workbook.add_worksheet()

        # write the column headers
        excel_headers = ["PRODUCT", "NAME", "SEQUENCE (5' > 3')", "LENGTH", "SYNTHESIS SCALE/(MIN) DELIVERED QUANTITY", "PURIFICATION", "FORMAT", "5' MODIFICATION", "3' MODIFICATION", "INTERNAL MODIFICATION 1 (indicate with 1: AAA1AAA)", "INTERNAL MODIFICATION 2 (indicate with 2: AAA2AAA)", "QC/ADDITIONAL QC"]
        col = 0
        for header in excel_headers:
            worksheet.write(0, col, header)
            col += 1

        probes_to_write = self.top_probes if self.top_probes is not None else self.get_top_probes(top_x = num_probes)
        # write the probes
        row, col, cntr = 1, 0, 0
        for probe_pair in probes_to_write:
            # isolate the padlock and splints
            padlock = self.padlock_leader[1:] + probe_pair[0] + self.padlock_barcode
            splint = probe_pair[1] + self.splint_connector

            # write in the padlock row
            worksheet.write(row, col, product)
            worksheet.write(row, col+1, self.gene_name + "_padlock_" + str(cntr))
            worksheet.write(row, col+2, padlock)
            worksheet.write(row, col+3, len(padlock))
            worksheet.write(row, col+4, syn_scale)
            worksheet.write(row, col+5, pur)
            worksheet.write(row, col+6, form)
            worksheet.write(row, col+7, 'Phosphate') # 5' modification for padlock

            # write the splint row
            worksheet.write(row+1, col, product)
            worksheet.write(row+1, col+1, self.gene_name + "_splint_" + str(cntr))
            worksheet.write(row+1, col+2, splint)
            worksheet.write(row+1, col+3, len(splint))
            worksheet.write(row+1, col+4, syn_scale)
            worksheet.write(row+1, col+5, pur)
            worksheet.write(row+1, col+6, form)

            # increment the rows and the counter
            row += 2
            cntr += 1
        # finished writing, close the workbook
        workbook.close()

    def sanity_check_probes(self, filename, num_probes=5, probes_to_hl=[]):
        """
        Creates a html file with the sequence in fasta format. Probe pairs are highlighted (splint and padlock). Useful for sanity checking the probes
        ---
        Input:
        filename = name and path of the docx file to write to. Required
        n_probes = number of probes to highlight. Default is 5
        probes_to_hl =  a list of probe pairs to highlight in the fasta file (NOT the probe targets in the sequence)

        Output:
        an html file with the fasta sequence and probe pairs highlighted
        """
        import webbrowser
        seq = self.sequence.upper()
        # if no probes are provided, get the top_x probes
        if len(probes_to_hl) is 0:
            # get the reverse complement of the top n probes
            top_probes = self.top_probes if self.top_probes is not None else self.get_top_probes(top_x = num_probes)
            for probe_pair in top_probes:
                probes_to_hl.append((probe_pair[0], probe_pair[1]))
        # reverse complement all of the probes
        probes_to_hl = [(rev_comp(i[0]).upper(), rev_comp(i[1]).upper()) for i in probes_to_hl]
        # split up the fasta file by the probes
        file_start = '<html>\n<head></head>\n<body><h1>{}</h1><div style="width: 50%; word-wrap: break-word;"><p>'.format(self.gene_name.upper())
        file_end = '</p><p>Padlock targets are in <span style="background-color: #FFFF00"> yellow</span>.</p><p>Splint targets are in <span style="background-color: #008000">green</span>.</p></div></body></html>'
        las_idx = 0
        for probe_pair in probes_to_hl:
            padlock = probe_pair[0]
            splint = probe_pair[1]

            # create the formatting strings
            formatted_padlock = '<span style="background-color: #FFFF00">{}</span>'.format(padlock)
            formatted_splint = '<span style="background-color: #008000">{}</span>'.format(splint)
            # replace the substring with html formatting appended
            seq = seq.replace(padlock, formatted_padlock)
            seq = seq.replace(splint, formatted_splint)

        file_contents = file_start + seq + file_end
        with open(filename, 'w') as f:
        	f.write(file_contents)
        	f.close()

        # print("Document successfully written to {}".format(filename))
        webbrowser.open_new_tab('file://' + os.path.realpath(filename))

    def check_overlap(self, p1, p2):
        """
        Checks if probes p1 and p2 overlap. True if p1 and p2 overlap and false otherwise.

        Input:
        p1 = probe sequence 1
        p2 = probe sequence 2

        Output:
        returns true if the probes overlap and false otherwise.
        """
        revcomp_seq = rev_comp(self.sequence)

        # get start and end idices for p1 and p2
        p1_startix = revcomp_seq.index(p1)
        p1_endix = p1_startix + len(p1) - 1
        p2_startix = revcomp_seq.index(p2)
        p2_endix = p2_startix + len(p2) - 1

        # debugging
        #print("p1 start {}, p1 end {}, p2 start {}, p2 end {}".format(p1_startix, p1_endix, p2_startix, p2_endix))

        # if p1 ends before p2 starts and if p2 ends before p1 starts
        if p1_endix < p2_startix or p2_endix < p2_startix:
            return(False)
        # otherwise assume overlap
        return(True)

    def check_overlaps(self, pair1, pair2):
        """
        Checks if probe pair 1 and probe pair 2 do not overlap.

        Input:
        pair1 = a probe pair in the form of a tuple (e.g. (padlock, splint,...))
        pair2 = a probe pair in the same form as pair1

        Output:
        returns true if the probe pairs overlap and false otherwise.
        """
        # check to make sure that:
        # 1. splints do not overlap
        cond1 = self.check_overlap(pair1[1], pair2[1])
        # 2. padlocks do not overlap
        cond2 = self.check_overlap(pair1[0], pair2[0])
        # 3. splint 1 does not overlap with padlock 2
        cond3 = self.check_overlap(pair1[1], pair2[0])
        # 4. padlock 1 does not overlap with splint 2
        cond4 = self.check_overlap(pair1[0], pair2[1])
        # check all conditions and add to list
        if (not cond1) and (not cond2) and (not cond3) and (not cond4):
            return(False)
        return(True)


# if __name__ == "__main__":
#     main()