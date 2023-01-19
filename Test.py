#!/usr/bin/python
# Written in Python 3.7 in 2022 by A.L.O. Gaenssle
# Counts phylum and gene names

import math
import statistics

# Own modules
import Import_Export as IE

# Create a fasta file of all sequences (>ID [Species]\nSequence)
def CreateFasta_wholeSequence(ProteinList, OutputFile, Ask=True):
	Fasta = []
	for Protein in ProteinList:
		SubString = ">" + Protein[0] + " [" + Protein[2] + "]\n" + Protein[5] + "\n"
		Fasta.append(SubString)
	IE.ExportList(Fasta, File.rsplit(".", 1)[0] + ".fasta", Ask=Ask)


Folder = "DUF1735"
Name = "DUF1735"
# TypeList = ["Start" ,"AA"]
Ask = False
DBList = ["KEGG"]

for DB in DBList:
	File = Folder + "/Output/" + Name + "_" + DB + "_Domains_Cutoff-0.0001.txt"
	ProteinList = IE.ImportNestedList(File)
	CreateFasta(ProteinList, File, Ask=Ask)
