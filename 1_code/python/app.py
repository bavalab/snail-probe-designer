# !/usr/bin/python
# -*- coding: iso-8859-1 -*-

# SNAIL Probe Designer GUI and application
# Author: Eric Cramer <eric.cramer@curie.fr>

# A GUI for the SNAIL Probe Designer
# Allows you to upload a fasta file or input a fasta sequence to generate smFISH SNAIL probes 

import tkinter, os
import snail_probe_designer as SPD

class spd_gui(tkinter.Tk):
    """
    A GUI for the SNAIL Probe Designer

    The application allows you to upload a fasta file or input a fasta sequence to generate smFISH SNAIL probes 

    """
    def __init__(self, master=None):
        tkinter.Tk.__init__(self, master)
        self.master = master
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
        self.sequence.grid(column=0, row=0, sticky='EW')
        self.sequenceVar.set("Enter sequence here")

        # add a text box for the gene name
        self.geneNameVar = tkinter.StringVar()
        self.geneName = tkinter.Entry(self, textvariable=self.geneNameVar)
        self.geneName.grid(column=1, row=0, sticky='EW')
        self.geneNameVar.set("Enter gene name here")

        # add a button to generate probes
        goButton = tkinter.Button(self, text="Go", command=self.goClient)
        goButton.grid(column=2, row=0)

        # add a button to quit
        quitButton = tkinter.Button(self, text="Quit", command=self.quitClient)
        quitButton.grid(column=3, row=0)

        # add text boxes for the melting temperature
        self.minMeltingTempVar = tkinter.StringVar()
        self.minMeltingTemp = tkinter.Entry(self, textvariable=self.minMeltingTempVar)
        self.minMeltingTemp.grid(column=0, row=1, sticky='EW')
        self.minMeltingTempVar.set("Minimum melting temperature (C)")
        self.maxMeltingTempVar = tkinter.StringVar()
        self.maxMeltingTemp = tkinter.Entry(self, textvariable=self.maxMeltingTempVar)
        self.maxMeltingTemp.grid(column=1, row=1, sticky='EW')
        self.maxMeltingTempVar.set("Maximum melting temperature (C)")

        # add text boxes for the GC content
        self.gContentVar = tkinter.StringVar()
        self.gContent = tkinter.Entry(self, textvariable=self.gContentVar)
        self.gContent.grid(column=2, row=1, sticky='EW')
        self.gContentVar.set("G Content")
        self.cContentVar = tkinter.StringVar()
        self.cContent = tkinter.Entry(self, textvariable=self.cContentVar)
        self.cContent.grid(column=3, row=1, sticky='EW')
        self.cContentVar.set("C Content")

        # add a label for instructions
        self.instructionsLabelVar = tkinter.StringVar()
        instructionsLabel = tkinter.Label(self, 
        	textvariable=self.instructionsLabelVar, 
        	anchor='w', 
        	fg='black', 
        	bg='white', wraplength=300, justify=tkinter.LEFT)
        instructionsLabel.grid(column=0, row=2, columnspan=4, sticky='EW')
        instructions = "\nWelcome to the SNAIL Probe Designer\nTo use the program:\n1. Input the sequence you want to design probes for into the first box\n2. Input the name of the sequence into the second box\n3. Click 'Go'\n\nThe most highly rated probe targets will be displayed in this window, and the full probes will be written to a csv and Eurogentec formatted excel file. Finally, a 'sanity check' of the sequence with the probe targets highlighted will be written to a separate html file."
        self.instructionsLabelVar.set(instructions)

        # add a label for output
        self.outputLabelVar = tkinter.StringVar()
        output_label = tkinter.Label(self, textvariable=self.outputLabelVar, anchor='w', fg='blue', bg='white')
        output_label.grid(column=0, row=3, columnspan=4, sticky='EW')


    def quitClient(self):
        exit()

    # def getParams(self):
    #     seq = self.sequenceVar.get()
    #     geneName = self.geneNameVar.get()
    #     minT = self.minMeltingTempVar.get()
    #     maxT = self.maxMeltingTempVar.get()
    #     gCont = self.gContentVar.get()
    #     cCont = self.cContentVar.get()

    #     # check to make sure the input is numeric    	


    def goClient(self):
    	# get the parameters from the GUI
        seq = self.sequenceVar.get()
        geneName = self.geneNameVar.get()
        minT = float(self.minMeltingTempVar.get())
        maxT = float(self.maxMeltingTempVar.get())
        tm = (minT, maxT)
        gCont = float(self.gContentVar.get())
        cCont = float(self.cContentVar.get())
        gc = (gCont, cCont)

        # uses the default settings for gc, tm, size, and separation
        spd=SPD.snail_probe_designer(tm=tm, gc=gc) 

        spd.prime(seq, geneName)
        spd.get_kmers()
        spd.score_kmers()
        spd.write_probes_to_csv("../../2_output/temp/{}_probes.csv".format(geneName.lower()), num_probes=3)
        spd.write_probes_to_eurogentec("../../2_output/temp/{}_probes_eurogentec.xlsx".format(geneName.lower()))
        print("{} potential probe pairs found.".format(len(spd.scored_probe_pairs)))

        # print out the top probe pairs, one per line
        top_probes = spd.get_top_probes()
        out_str = "The best {} probe pairs are:\n".format(len(top_probes))
        for i in range(len(top_probes)):
            out_str += str(top_probes[i]) + "\n"
        self.outputLabelVar.set(out_str)
        spd.sanity_check_probes('../../2_output/temp/sanity_check.html')

    def OnPressEnter(self, event):
        """
        Event handler for when the Enter key is pressed.

        """
        self.labelVariable.set(self.entryVariable.get() + " You pressed Enter!")
        self.entry.focus_set()
        self.entry.selection_range(0, tkinter.END)

if __name__ == "__main__":
    app = spd_gui(None)
    app.title("SNAIL Probe Designer")
    app.mainloop()