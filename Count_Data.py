#!/usr/bin/python
# Written in Python 3.7 in 2022 by A.L.O. Gaenssle
# Counts phylum and gene names

# Own modules
import Import_Export as IE

# Count the occurence of the same taxonomy in the file
def CountTaxonomy(Table, Header, Folder, Name, DB):
	Header = Header.strip().split("\t")
	Index = Header.index("Taxonomy")
	Count = {}
	for Gene in Table:
		if Gene[Index] in Count:
			Count[Gene[Index]] += 1
		else:
			Count[Gene[Index]] = 1
	sortedList_Count = sorted(Count.items(), key=lambda x:x[1], reverse=True)
	sortedDict_Count = dict(sortedList_Count)
	OutputFile = Folder + "/Output/" + Name + "_" + DB + "_CountTaxonomy.txt"
	IE.ExportDictionary(sortedDict_Count, OutputFile, Header="Taxonomy\tCount\n")

# Count the occurence of the same domain motif (architecture) in the file
def CountMotif(Table, Folder, Name, DB):
	Count = {}
	for Gene in Table:
		if Gene[-1] in Count:
			Count[Gene[-1]] += 1
		else:
			Count[Gene[-1]] = 1
	sortedList_Count = sorted(Count.items(), key=lambda x:x[1], reverse=True)
	sortedDict_Count = dict(sortedList_Count)
	OutputFile = Folder + "/Output/" + Name + "_" + DB + "_CountMotif.txt"
	IE.ExportDictionary(sortedDict_Count, OutputFile, Header="Motif\tCount\n")
