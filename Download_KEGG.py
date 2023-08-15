#!/usr/bin/python
# Written in Python 3.7 in 2022 by A.L.O. Gaenssle
# Downloads data from the KEGG database

import re
import Bio
from Bio.KEGG import REST
import urllib.request
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

# Own modules
import Import_Export as IE

##------------------------------------------------------
## DOWNLOAD FUNCTIONS
##------------------------------------------------------
# Download Info for each protein from KEGG
def GetDetailedData(Entry, ID):
	Dict = {"ID": ID,"AASeq": ""}
	inAASeq = False
	for Line in Entry:
		Line = re.sub("\s\s+" , " ", Line)
		Line = Line.strip()
		if inAASeq == True:
			if Line.startswith("NTSEQ"):
				break
			else:
				Dict["AASeq"] += Line.strip()
		if inAASeq == False:
			if Line.startswith("ORGANISM") or Line.startswith("VIRUS"):
				Line = Line.split(" ",1)[1].strip()
				Dict["Organism"] = Line.split(" ",1)[1]
			elif "UniProt" in Line:
				Dict["UniProt"] = Line.split(" ",1)[1]
			elif Line.startswith("AASEQ"):
				Dict["AALength"] = Line.split(" ",1)[1]
				inAASeq = True
	# ProteinData = [ID, UniProt, Organism, AALength, AASeq]
	return(Dict)

# Download protein entries from KEGG -> in chunks of 10 gene IDs --> KEGG-get
def DownloadProteinEntries(Chunked_List):
	Data = []
	Entry = []
	print("Download protein info, set of", Chunked_List[0], ". . .")
	Download = REST.kegg_get(Chunked_List).read()
	Download = Download.split("\n")
	Count = 0
	for Line in Download:
		if Line.startswith("///"):
			ProteinData = GetDetailedData(Entry, Chunked_List[Count])
			Data.append(ProteinData)
			Entry = []
			Count += 1
		else:
			Entry.append(Line)
	return(Data)

# Download all genome taxonomy from KEGG --> KEGG-list
def DownloadOrganismsTemp(Name="organism"):
	print("Download organism taxonomy. . .")
	Entry = REST.kegg_list(Name).read()
	List = Entry.split("\n")
	Table = {}
	while("" in List) :
		List.remove("")
	for ID in List:
		SubList = ID.split("\t")
		Taxonomy = "-".join(SubList[3].split(";",3)[1:3]).split(" - ",1)[0]
		Table[SubList[1]] = Taxonomy
	return(Table)

# Download domain motifs (architecture) for each gene ID from KEGG
def DownloadMotif(ID):
	Data = []
	NewDomains = []
	Organism = ""
	SubList = []
	url = "https://www.kegg.jp/ssdb-bin/ssdb_motif?kid=" + ID
	with urllib.request.urlopen(url) as File:
		for Line in File:
			Line = Line.decode("utf-8").strip().replace(" : ", "")
			Line = re.sub('<[^>]*>', '|', Line)
			Line = Line.split("|")
			Line = [i for i in Line if i]
			if Line != []:
				Data.append(Line)
		for Index in range(len(Data)):
			Line = Data[Index]
			if Line[0].startswith("Organism"):
				try:
					Organism = Line[1]
				except IndexError:
					Organism = "N.D."
			if Line[0].startswith("pf:"):
				Line[0] = Line[0].split(":",1)[1]
				if Line[3] == "&nbsp;":
					Line[3] = "N.D."
					SubList = [ID, Organism] + Line
				else:
					SubList = [ID, Organism] + Line + Data[Index+1]
				SubList[3] = "{:03d}".format(int(SubList[3]))
				if SubList not in NewDomains:
					NewDomains.append(SubList)
	NewDomains.sort(key = lambda x: x[3])
	print(ID, "downloaded")
	return(NewDomains)
