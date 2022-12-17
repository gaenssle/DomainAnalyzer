#!/usr/bin/python
# Written in Python 3.7 in 2022 by A.L.O. Gaenssle
# Counts phylum and gene names

import math
import statistics

# Own modules
import Import_Export as IE

def FindOutliers(Domains, Motifs, Header, Folder, Name, DB, Ask=True):
	IndexStart = Header.split("\t").index("Start")
	IndexLength = Header.split("\t").index("AA")
	Position = {"First":[], "Second":[], "Third":[], "Outlier":[]}
	for Line in Domains:
		Start = int(Line[IndexStart])
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


Folder = "DUF1735"
Name = "DUF1735"
TypeList = ["Start" ,"AA"]
Ask = False
DBList = ["KEGG", "UniProt"]

for DB in DBList:
	File = Folder + "/Output/" + Name + "_" + DB
	Motifs = IE.ImportNestedDictionary(File + "_Domains_Motifs.txt")
	Domains, Header = IE.ImportNestedList(File + "_Domains_onlyDUF1735.txt", getHeader=True)
	Position = FindOutliers(Domains, Motifs, Header, Folder, Name, DB, Ask=Ask)
	for Type in TypeList:
		DomainStatistics(Position, Folder, Name, DB, Type, Ask=Ask)
