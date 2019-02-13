# !/usr/bin/python
# -*- coding: iso-8859-1 -*-

# SNAIL Probe Designer GUI and application
# Author: Eric Cramer <eric.cramer@curie.fr>

# A GUI for the SNAIL Probe Designer
# Allows you to upload a fasta file or input a fasta sequence to generate smFISH SNAIL probes 

# TODOs:
# 1a. change the GC content selection to a slider for G content
# 1b. derive c content from 100-G content - DONE
# 2. add a directory selection dialog to get the output directory - DONE
# 3. display the names and locations of the files that are saved - DONE
# 4. fix the sanity checking bug --> keeps showing the sanity checker from the previous iteration - DONE
# 5. add options to change padlock leader, splint connector, probe length (DONE), and barcoding
# 6. add file dialog for direct fasta file selection
# 7. improve robustness --> switch default arguments for spd constructor to **kwargs, add error checking
# 8. port to executable for mac and pc
# 9. fix style
# 10. write documentation

import tkinter, os, sys
from tkinter import filedialog
import snail_probe_designer as SPD

class spd_gui(tkinter.Tk):
    """
    A GUI for the SNAIL Probe Designer

    The application allows you to upload a fasta file or input a fasta sequence to generate smFISH SNAIL probes 

    """
    def __init__(self, master=None):
        tkinter.Tk.__init__(self, master)
        self.master = master
        # initialize a snail probe designer
        # uses the default settings for gc, tm, size, and separation
        self.spd = SPD.snail_probe_designer() 
        # init the GUI
        self.initialize()

    def initialize(self):
        """
        Initializes the arrangement of widgets in the GUI and adds event handling.

        """
        # let the widget take th space of the root window
        self.grid()
        self.grid_columnconfigure(0, weight=1)

        # add a text box for the sequence
        self.sequenceVar = tkinter.StringVar()
        self.sequence = tkinter.Entry(self, textvariable=self.sequenceVar)
        self.sequence.grid(column=0, row=0, columnspan=2, sticky='EW')
        self.sequenceVar.set("Enter sequence here")

        # add a text box for the gene name
        self.geneNameVar = tkinter.StringVar()
        self.geneName = tkinter.Entry(self, textvariable=self.geneNameVar)
        self.geneName.grid(column=2, row=0, columnspan=2, sticky='EW')
        self.geneNameVar.set("Enter gene name here")

        # add a button to generate probes
        goButton = tkinter.Button(self, text="Go", command=self.goClient)
        goButton.grid(column=4, row=0)
        goButton.bind("<Return>", self.goClient)

        # add a button to choose an output directory
        self.saveFolder = ""
        saveButton = tkinter.Button(self, text="Save", command=self.saveClient)
        saveButton.grid(column=4, row=1)

        # add a button to quit
        quitButton = tkinter.Button(self, text="Quit", command=self.quitClient)
        quitButton.grid(column=4, row=7)

        # add text boxes for the melting temperature
        minMeltingTempLabel = tkinter.Label(self, text="Minimum melting temperature (C)", anchor='w', fg='black', bg='white')
        minMeltingTempLabel.grid(column=0, row=1, sticky='EW')
        self.minMeltingTempVar = tkinter.StringVar()
        self.minMeltingTemp = tkinter.Entry(self, textvariable=self.minMeltingTempVar)
        self.minMeltingTemp.grid(column=1, row=1, sticky='EW')
        self.minMeltingTempVar.set(50)
        maxMeltingTempLabel = tkinter.Label(self, text="Maximum melting temperature (C)", anchor='w', fg='black', bg='white')
        maxMeltingTempLabel.grid(column=2, row=1, sticky='EW')
        self.maxMeltingTempVar = tkinter.StringVar()
        self.maxMeltingTemp = tkinter.Entry(self, textvariable=self.maxMeltingTempVar)
        self.maxMeltingTemp.grid(column=3, row=1, sticky='EW')
        self.maxMeltingTempVar.set(65)

        # add text boxes for the probe lengths
        minLengthLabel = tkinter.Label(self, text="Minimum Length (bp)", anchor='w', fg='black', bg='white')
        minLengthLabel.grid(column=0, row=2, sticky='EW')
        self.minLengthVar = tkinter.StringVar()
        self.minLength = tkinter.Entry(self, textvariable=self.minLengthVar)
        self.minLength.grid(column=1, row=2, sticky='EW')
        self.minLengthVar.set(18)
        maxLengthLabel = tkinter.Label(self, text="Maximum Length (bp)", anchor='w', fg='black', bg='white')
        maxLengthLabel.grid(column=2, row=2, sticky='EW')
        self.maxLengthVar = tkinter.StringVar()
        self.maxLength = tkinter.Entry(self, textvariable=self.maxLengthVar)
        self.maxLength.grid(column=3, row=2, sticky='EW')
        self.maxLengthVar.set(24)

        # add text boxes for the G content
        gContentLabel = tkinter.Label(self, text="G Content (%)", anchor='w', fg='black', bg='white')
        gContentLabel.grid(column=0, row=3, sticky='EW')
        self.gContentVar = tkinter.StringVar()
        self.gContent = tkinter.Entry(self, textvariable=self.gContentVar)
        self.gContent.grid(column=1, row=3, sticky='EW')
        self.gContentVar.set(30)

        # add text boxes for the max separation
        separationLabel = tkinter.Label(self, text="Maximum Separation (bp)", anchor='w', fg='black', bg='white')
        separationLabel.grid(column=2, row=3, sticky='EW')
        self.separationVar = tkinter.StringVar()
        self.separation = tkinter.Entry(self, textvariable=self.separationVar)
        self.separation.grid(column=3, row=3, sticky='EW')
        self.separationVar.set(4)

        # add text box for number of probes to return
        numProbesLabel = tkinter.Label(self, text="Number of Probes", anchor='w', fg='black', bg='white')
        numProbesLabel.grid(column=0, row=4, sticky='EW')
        self.numProbesVar = tkinter.StringVar()
        self.numProbes = tkinter.Entry(self, textvariable=self.numProbesVar)
        self.numProbes.grid(column=1, row=4, sticky='EW')
        self.numProbesVar.set(5)

        # add a label for instructions
        self.instructionsLabelVar = tkinter.StringVar()
        instructionsLabel = tkinter.Label(self, 
        	textvariable=self.instructionsLabelVar, 
        	anchor='w', 
        	fg='black', 
        	bg='white', wraplength=500, justify=tkinter.LEFT)
        instructionsLabel.grid(column=0, row=5, columnspan=4, sticky='EW')
        instructions = "\nWelcome to the SNAIL Probe Designer\nTo use the program:\n1. Input the sequence you want to design probes for into the first box\n2. Input the name of the sequence into the second box\n3. Click 'Go'\n\nThe indicated number of probe targets will be displayed in his window in the order of best hybridization score. The same fully structured probes can be written to csv and Eurogentec formatted Excel files by clicking 'Save'. Additionally, a 'sanity check' of the sequence with the probe targets highlighted will be written to a separate html file in the same location. File names and locations will be displayed beneath the probe results.\n"
        self.instructionsLabelVar.set(instructions)

        # add a label for output
        self.outputLabelVar = tkinter.StringVar()
        outputLabel = tkinter.Label(self, textvariable=self.outputLabelVar, anchor='w', fg='blue', bg='white', wraplength=500, justify=tkinter.LEFT)
        outputLabel.grid(column=0, row=6, columnspan=4, sticky='EW')

        # add a label to show what files were saved
        self.fileOutputVar = tkinter.StringVar()
        fileOutputLabel = tkinter.Label(self, textvariable=self.fileOutputVar, anchor='w', fg='green', bg='white', wraplength=500, justify=tkinter.LEFT)
        fileOutputLabel.grid(column=0, row=7, columnspan=4, sticky='EW')

        # add a label for the copyright
        self.copyrightVar = tkinter.StringVar()
        copyrightLabel = tkinter.Label(self, textvariable=self.copyrightVar, anchor='w', fg='black', bg='white')
        copyrightLabel.grid(column=0, row=8, columnspan=5, sticky='EW')
        self.copyrightVar.set(u"\u00A9 2019 Bava Lab, Institut Curie")

        # finalize tkinter setup
        self.update()
        self.geometry(self.geometry())

    def quitClient(self):
        sys.exit()

    def saveClient(self):
        # get the directory and set file names
        dirname = tkinter.filedialog.askdirectory()
        self.saveFolder = dirname
        csv_filename = self.saveFolder+"/{}_probes.csv".format(self.geneNameVar.get().lower())
        excel_filename = self.saveFolder+"/{}s_probes_eurogentec.xlsx".format(self.geneNameVar.get().lower())
        html_filename = self.saveFolder+"/{}_sanity_check.html".format(self.geneNameVar.get().lower())
        out_str = "Files written:\n"+csv_filename+"\n"+excel_filename+"\n"+html_filename

        # get the number of probes to save
        num_probes = int(self.numProbesVar.get())

        # write out the files
        self.spd.write_probes_to_csv(csv_filename, num_probes=num_probes)
        self.spd.write_probes_to_eurogentec(excel_filename, num_probes=num_probes)
        self.spd.sanity_check_probes(html_filename, num_probes=num_probes)

        # show the ouput files
        self.fileOutputVar.set(out_str)

    def goClient(self):
    	# get the parameters from the GUI 
        # set the gc, tm, size, and separation
        self.spd.probe_tm = (float(self.minMeltingTempVar.get()), float(self.maxMeltingTempVar.get()))
        self.spd.probe_gc = (float(self.gContentVar.get()), 100-float(self.gContentVar.get())) 
        self.spd.probe_sep = int(self.separationVar.get())
        self.spd.probe_size = (int(self.minLengthVar.get()), int(self.maxLengthVar.get()))
        self.spd.prime(self.sequenceVar.get(), self.geneNameVar.get())
        self.spd.get_kmers()
        self.spd.score_kmers()

        print("{} potential probe pairs found.".format(len(self.spd.scored_probe_pairs)))

        # print out the top probe pairs, one per line
        top_probes = self.spd.get_top_probes(top_x=int(self.numProbesVar.get()))
        out_str = "No probe pairs meeting criteria exactly found."
        if len(top_probes) > 0:
            out_str = "The best {} probe pairs are:\n".format(len(top_probes))
            for i in range(len(top_probes)):
                out_str += str(top_probes[i]) + "\n"
        self.outputLabelVar.set(out_str)

    def OnPressEnter(self, event):
        """
        Event handler for when the Enter key is pressed.

        """
        self.labelVariable.set(self.entryVariable.get() + " You pressed Enter!")
        self.entry.focus_set()
        self.entry.selection_range(0, tkinter.END)

if __name__ == "__main__":
    app = spd_gui(None)
    app.title("SNAIL Probe Designer v0.0.1")
    app.mainloop()