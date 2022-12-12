#!/usr/bin/python
# Written in Python 3.7 in 2022 by A.L.O. Gaenssle
# Counts phylum and gene names

import math

# Own modules
import Import_Export as IE

# Count the occurence of the same taxonomy in the file
def CountTaxonomy(Table, Header, Folder, Name, DB, Ask=True):
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
	IE.ExportDictionary(sortedDict_Count, OutputFile, Header="Taxonomy\tCount\n", Ask=Ask)

# Count the occurence of the same domain motif (architecture) in the file
def CountMotif(Table, Folder, Name, DB, Ask=True):
	Count = {}
	for Gene in Table:
		if Gene[-1] in Count:
			Count[Gene[-1]] += 1
		else:
			Count[Gene[-1]] = 1
	sortedList_Count = sorted(Count.items(), key=lambda x:x[1], reverse=True)
	sortedDict_Count = dict(sortedList_Count)
	OutputFile = Folder + "/Output/" + Name + "_" + DB + "_CountMotif.txt"
	IE.ExportDictionary(sortedDict_Count, OutputFile, Header="Motif\tCount\n", Ask=Ask)

# Group the numbers and count their occurence (for Start position and sequnce length)
def CountNumber(Table, Header, Folder, Name, DB, Type, Ask=True, Interval=10):
	Index = Header.split("\t").index(Type)
	Count = {}
	for Line in Table:
	    Start = Line[Index]
	    Floor = math.floor(int(Start)/Interval)*Interval
	    Group = str("{:03d}".format(Floor)) + "-" + str("{:03d}".format(Floor+Interval-1))
	    if Group in Count:
	        Count[Group] += 1
	    else:
	        Count[Group] = 1
	sortedList_Count = sorted(Count.items(), key=lambda x:x[0])
	sortedDict_Count = dict(sortedList_Count)
	if Type == "AA":
		Type = "Length"
	OutputFile = Folder + "/Output/" + Name + "_" + DB + "_Count" + Type + ".txt"
	IE.ExportDictionary(sortedDict_Count, OutputFile, Header=Type+"\tCount\n", Ask=Ask)
