# FISH Probe Designers
Created by Francesca Zaniboni and Eric Cramer
Created on 2019-02-08 14:15:38.845405
## Project Description 
A program to identify and construct sequences that can be used to design probe systems Single Molecule Fluorescence In Situ Hybridization experiments (smFISH). Specifically, the probe designer creates probes for SNAIL Proximity Ligation Assays. The SNAIL probe scheme is flexible, specific, and highly multiplexible, as demonstrated by [Wang et al (2018)](http://science.sciencemag.org/content/361/6400/eaat5691/tab-figures-data) in their STARMap method. SNAIL probes can be barcoded and read in situ using the SEDAL probe method.

From Figure 1A of the STARMap paper:
![SNAIL probe anatomy from Wang et al 2018](3_docs/img/snail-probe-example.PNG)
>Design of SNAIL probes (one component of STARmap): each primer or padlock probe has 19-25 nt (labeled by blue double-headed arrows) to hybridize with target RNA with a designed Tm of 60oC, while the complementary sequence between the primer and padlock is only 6 nt on each arm (labeled by red doubleheaded arrow) with Tm below room temperature, so that primer-padlock DNA-DNA hybridization is negligible during DNA-RNA hybridization at 40oC, but allows DNA ligation by T4 DNA ligase in the following step. Tm, melting temperature of nucleic acids. 

And from Figure 3A of the STARMap paper:
![SEDAL sequencing example from WANg et al 2018](3_docs/img/sedal-probe-example.PNG)
>SEDAL involves a T4 DNA ligase with activity strongly hindered by base mismatches, and two kinds of sequencing probes: reading probes that set the base position to be interrogated, and fluorescent decoding probes that transduce base information into colors for imaging. Unlike other sequencing-by-ligation methods which use preannealed reading probes (or equivalent), the reading probe in SEDAL is short (11 nt, with Tm near room temperature), partially degenerate (as shown here, for cycle 4, the first two base at 5’ end are N, equal amount mixture of A, T, C and G), and mixed with decoding probes and T4 DNA ligase for a one-step reaction. At room temperature, the reading probe remains in a dynamic state of annealing with and detaching from the DNA template. Only when the reading probe perfectly matches the DNA template, T4 DNA ligase (blue) ligates it to the fluorescent 8-nt decoding probe. The short reading and decoding probes are then washed away, leaving fluorescent 19-nt products stably hybridized to the DNA amplicon for imaging. For the next cycle, previous fluorescent products are stripped and the reading probe includes one more degenerate base to shift the reading frame by one base (Fig. 1E). 5’P: 5’ phosphate. 3’InvT: 3’ inverted dT base that prevents self-ligation of the reading probe. 3’OH: 3’ hydroxyl group. 

The SEDAL and SNAIL combination produces STARMap, a three-part system to detect _n_ mRNA transcripts _in situ_ with fluorescence microscopy. First tissue is cleared of lipids and unfixed proteins using the [CLARITY](http://clarityresourcecenter.org/) method. Second, a set of SNAIL splint and padlock probes are hybridized to the target RNAs. Finally, a series of _reading_ and _fluorescent_ probes are used to readout a barcode embedded in the padlock portion of the SNAIL probe (thus uniquely identifying the RNA transcript). 

## TODOs
+ change the GC content selection to a slider for G content in app
+ derive c content from 100-G content - DONE
+ add a directory selection dialog to get the output directory - DONE
+ display the names and locations of the files that are saved - DONE
+ re-format the sanity checker
+ add probe length adjustments (DONE)
+ add file dialog for direct fasta file selection
+ add multiple fasta sequence support
+ switch default arguments for spd constructor to kwargs
+ add error and variable type checking
+ add "plug and play" for probe designs
	- Ex. switching between SNAIL and PLAYR style probes
	- add options to change padlock leader and splint connector
	- barcoding options
		+ add imports for barcode to gene lookup table
+ port to executable for mac and pc
+ fix style
+ write documentation

## General Notes (to self)
For compilation with `pyinstaller`, use the conda prompt:  
`activate py34`  
`pyinstaller --onefile --windowed --distpath=D:\coding-projects\snail-probe-designer D:\coding-projects\snail-probe-designer\1_code\python\app.py`