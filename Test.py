#!/usr/bin/python
# Written in Python 3.7 in 2022 by A.L.O. Gaenssle
# Counts phylum and gene names

import math
import statistics

# Own modules
import Import_Export as IE

Folder = "DUF1735"
Name = "DUF1735"
DB = "UniProt"
File = Folder + "/Output/" + Name + "_" + DB

Motifs = IE.ImportNestedDictionary(File + "_Domains_Motifs.txt")
Domains, Header = IE.ImportNestedList(File + "_Domains_onlyDUF1735.txt", getHeader=True)
Type = "Start"
Ask = False

def FindOutliers(Domains, Motifs, header, Folder, Name, DB, Ask=True):
	Index = Header.split("\t").index(Type)
	Position = {"First":[], "Second":[], "Third":[], "Outlier":[]}
	Statistics = {"First":[], "Second":[], "Third":[], "Outlier":[]}
	for Line in Domains:
		Start = int(Line[Index])
		if Start >= 100:
			try:
				ID, Suffix = Line[0].rsplit("_",1)
				if len(Suffix) == 1:
					if int(Suffix) == 2:
						Position["Second"].append([Line[0], Start, Motifs[ID][-1]])
						Statistics["Second"].append(Start)
					elif int(Suffix) == 3:
						Position["Third"].append([Line[0], Start, Motifs[ID][-1]])
						Statistics["Third"].append(Start)
				else:
					Position["Outlier"].append([Line[0], Start, Motifs[Line[0]][-1]])
					Statistics["Outlier"].append(Start)
			except ValueError:
				Position["Outlier"].append([Line[0], Start, Motifs[Line[0]][-1]])
				Statistics["Outlier"].append(Start)
		else:
			try:
				ID, Suffix = Line[0].rsplit("_",1)
				if len(Suffix) == 1:
					Position["Outlier"].append([Line[0], Start, Motifs[ID][-1]])
					Statistics["Outlier"].append(Start)
				else:
					Position["First"].append([Line[0], Start, Motifs[Line[0]][-1]])
					Statistics["First"].append(Start)
			except ValueError:
				Position["First"].append([Line[0], Start, Motifs[Line[0]][-1]])
				Statistics["First"].append(Start)
	Summary = {}
	for Pos in Statistics:
		Count = len(Statistics[Pos])
		Mean = round(statistics.mean(Statistics[Pos]))
		Stdev = round(statistics.stdev(Statistics[Pos]))
		Min = min(Statistics[Pos])
		Max = max(Statistics[Pos])
		Summary[Pos] = [Count, Mean, Stdev, Min, Max]
		print(Summary[Pos])

	OutputFile = Folder + "/Output/" + Name + "_" + DB + ".txt"
	# IE.ExportDoubleNestedDictionary(Position, OutputFile, Header="Position\tID\tStart\tMotif\n", Add= "_DomainPosition", Ask=Ask)
	SumHeader = "Position\tCount\tMean\tStdev\tMin\tMax\n"
	IE.ExportNestedDictionary(Summary, OutputFile, Header=SumHeader, Add= "_DomainPosition_Summary", Ask=Ask)

FindOutliers(Domains, Motifs, Header, Folder, Name, DB, Ask=Ask)
