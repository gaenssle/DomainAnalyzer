#!/usr/bin/python
# Written in Python 3.7 in 2022 by A.L.O. Gaenssle
# Counts phylum and gene names

import math
import statistics

# Own modules
import Import_Export as IE

# Count the occurence of the same taxonomy in the file
def CountTaxonomy(Table, Header, Folder, Name, DB, Ask=True):
	Index = Header.split("\t").index("Taxonomy") - 1
	Count = {}
	for Gene in Table:
		if Table[Gene][Index] in Count:
			Count[Table[Gene][Index]] += 1
		else:
			Count[Table[Gene][Index]] = 1
	sortedList_Count = sorted(Count.items(), key=lambda x:x[1], reverse=True)
	sortedDict_Count = dict(sortedList_Count)
	OutputFile = Folder + "/Output/" + Name + "_" + DB + "_CountTaxonomy.txt"
	IE.ExportDictionary(sortedDict_Count, OutputFile, Header="Taxonomy\tCount\n", Ask=Ask)

# Count the occurence of the same domain motif (architecture) in the file
def CountMotif(Table, Folder, Name, DB, Ask=True):
	Count = {}
	for Gene in Table:
		if Table[Gene][-1] in Count:
			Count[Table[Gene][-1]] += 1
		else:
			Count[Table[Gene][-1]] = 1
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

# Group domains according to their start position (1st, 2nd, 3rd, outlier)
def GroupDomains(Domains, Motifs, Header, Folder, Name, DB, Ask=True):
	IndexStart = int(Header.split("\t").index("Start"))
	IndexLength = int(Header.split("\t").index("AA"))
	Position = {"First":[], "Second":[], "Third":[], "Outlier":[]}
	for Line in Domains:
		Start = int(Line[IndexStart])
		# print(IndexStart, IndexLength, Start)
		# print(Line[0], Start, Line[IndexLength])
		# print(Line)
		if Start >= 100:
			try:
				ID, Suffix = Line[0].rsplit("_",1)
				if len(Suffix) == 1:
					if int(Suffix) == 2:
						Position["Second"].append([Line[0], Start, Line[IndexLength], Motifs[ID][-1]])
					elif int(Suffix) == 3:
						Position["Third"].append([Line[0], Start, Line[IndexLength], Motifs[ID][-1]])
				else:
					Position["Outlier"].append([Line[0], Start, Line[IndexLength], Motifs[Line[0]][-1]])
			except ValueError:
				Position["Outlier"].append([Line[0], Start, Line[IndexLength], Motifs[Line[0]][-1]])
		else:
			try:
				ID, Suffix = Line[0].rsplit("_",1)
				if len(Suffix) == 1:
					Position["Outlier"].append([Line[0], Start, Line[IndexLength], Motifs[ID][-1]])
				else:
					Position["First"].append([Line[0], Start, Line[IndexLength], Motifs[Line[0]][-1]])
			except ValueError:
				Position["First"].append([Line[0], Start, Line[IndexLength], Motifs[Line[0]][-1]])
	OutputFile = Folder + "/Output/" + Name + "_" + DB + ".txt"
	IE.ExportDoubleNestedDictionary(Position, OutputFile, Header="Position\tID\tStart\tLength\tMotif\n", Add= "_grouped", Ask=Ask)
	return(Position)

# Determine statistics of the grouped domains
def DomainStatistics(Position, Folder, Name, DB, Type, Ask=True):
	Statistics = {}
	if Type == "Start":
		Index = 1
	else:
		Type = "Length"
		Index = 2
	for Group in Position:
		Statistics[Group] = []
		for Entry in Position[Group]:
				Statistics[Group].append(int(Entry[Index]))
	Summary = {}
	SumHeader = "Count\tPercent\tMean\tStdev\tMin\tMax\n"
	print(DB, Type)
	print(SumHeader.replace("\n", ""))
	AllCount = sum(len(x) for x in Statistics.values())
	for Group in Statistics:
		Count = len(Statistics[Group])
		Percent = round(Count/AllCount*100,2)
		Mean = round(statistics.mean(Statistics[Group]))
		Stdev = round(statistics.stdev(Statistics[Group]))
		Min = min(Statistics[Group])
		Max = max(Statistics[Group])
		Summary[Group] = [Count, Percent, Mean, Stdev, Min, Max]
		print(*Summary[Group], sep="\t")
	OutputFile = Folder + "/Output/" + Name + "_" + DB + "_" + Type + ".txt"
	IE.ExportNestedDictionary(Summary, OutputFile, Header="Position\t" + SumHeader, Add= "_Statistics", Ask=Ask)
