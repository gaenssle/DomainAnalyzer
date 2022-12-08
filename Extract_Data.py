#!/usr/bin/python
# Written in Python 3.7 in 2022 by A.L.O. Gaenssle
# Extract data from downloaded files

import copy

# Own modules
import Import_Export as IE

# Get the indices in the nested list for the domain info
def GetIndices(Header):
	Header = Header.strip().split("\t")
	Indices = [Header.index("Sequence")]
	Indices.append(Header.index("Domain"))
	Indices.append(Header.index("Start"))
	Indices.append(Header.index("End"))
	return(Indices)

# Determine the number of domains for each gene
def GetDomainCount(Gene, Check, Offset):
	Count = int((len(Gene) - Check)/Offset)
	return(Count)

# Extract the type and sequence of domains for each gene
def ExtractMotifs(GeneTable, Header, Name, Offset, OutputFile, Ask=True):
	Motifs = []
	Indices = GetIndices(Header)
	for Gene in GeneTable:
		if Name in Gene[Indices[1]:]:
			Array = [Gene[Indices[1]]] # Create Motif
			Count = GetDomainCount(Gene, Indices[1], Offset)
			for Amount in range(1, Count): # Cycle through domains in same gene
				Array.append(Gene[Indices[1]+Offset*Amount])
			Motifs.append(Gene[:Indices[0]] + ["-".join(Array)])
	HeaderMotif = "\t".join(Header.split("\t",Indices[0])[0:Indices[0]]) + "\tMotif\n"
	IE.ExportNestedList(Motifs, OutputFile, HeaderMotif, Add="_Motifs", Ask=Ask)
	return(Motifs)

# Convert the data from 1 gene/line (>= 1 domain) to 1 domain/line (>= 1 lines/gene)
def ExtractDetails(GeneTable, Header, Name, Offset, OutputFile, Ask=True):
	Details = []
	DetailsOnly = []
	Indices = GetIndices(Header)
	for Gene in GeneTable:
		if Name in Gene[Indices[1]:]: # Ignore the genes that don't contain a [Name] domain
			Count = GetDomainCount(Gene, Indices[1], Offset)
			ID = Gene[0]
			for Amount in range(Count): # Cycle through domains in same gene
				if Amount > 0:
					ID = Gene[0] + "_" + str(Amount+1)   # make unique IDs for domains
				Shift = Offset*Amount
				Sequence = Gene[Indices[0]]
				# print(Gene[0], Gene[Indices[2]+Shift], Gene[Indices[3]+Shift])
				try:
					Fragment = Sequence[int(Gene[Indices[2]+Shift])-1:int(Gene[Indices[3]+Shift])]
				except ValueError:
					print("Data for ID", ID, "could not read be in correctly and will be skipped")
					Fragment = "XX"
				SubList = [ID, Gene[1], Gene[Indices[1]+Shift]] + Gene[2:4] + [len(Fragment), Fragment]
				if Gene[Indices[1]+Shift] == Name:
					DetailsOnly.append(SubList)
				Details.append(SubList)
	HeaderPart = Header.split("\t",4)[0:4]
	HeaderDetails = HeaderPart[0:2] + ["Domains"] + HeaderPart[2:4] + ["AA", "Sequence"]
	IE.ExportNestedList(Details, OutputFile, "\t".join(HeaderDetails) + "\n", Add="_Details", Ask=Ask)
	IE.ExportNestedList(DetailsOnly, OutputFile, "\t".join(HeaderDetails) + "\n", Add="_only" + Name, Ask=Ask)
	return(DetailsOnly)

# Create a fasta file for all domains with the given name (>ID [Taxonomy]\nSequence)
def CreateFasta(DetailsOnly, OutputFile, Ask=True):
	Fasta = []
	for Gene in DetailsOnly:
		# Option A: "> kEGG-ID [Phylum-Genus]"
		# SubString = ">" + Gene[0] + " [" + Gene[-3] + "]\n" + Gene[-1] + "\n"
		# Option B: "> short-KEGG-ID [Species]"
		if "_" in Gene[0]:
			IDshort = Gene[0].split(":",1)[0] + ":" + Gene[0].split("_",1)[1]
		else:
			IDshort = Gene[0]
		SubString = ">" + IDshort + " [" + Gene[-4] + "]\n" + Gene[-1] + "\n"
		Fasta.append(SubString)
	IE.ExportList(Fasta, OutputFile.rsplit(".", 1)[0] + ".fasta", Ask=Ask)

# Add domain motifs to gene details from KEGG
def AddMotifKEGG(Genes, Motifs, Header, Name, Cutoff, OutputFile):
	Header = Header.rstrip() + "\tDomain\tStart\tEnd\tE-Value\n"
	Abbreviations = {}
	AllDomains = copy.deepcopy(Genes)
	GoodDomains = copy.deepcopy(Genes)
	for Domain in Motifs:
		try:
			SubList = Domain[2:5] + [Domain[-2]]
			AllDomains[Domain[0]].extend(SubList)
			if float(SubList[-1]) < Cutoff:
				GoodDomains[Domain[0]].extend(SubList)
			if Domain[2] not in Abbreviations:
				Abbreviations[Domain[2]] = Domain[5]
		except ValueError:
			pass
	GoodDomainList = []
	for Gene in GoodDomains:
		if Name in GoodDomains[Gene][5:]:
			GoodDomainList.append([Gene] + GoodDomains[Gene])
	sortedList_Abbr = sorted(Abbreviations.items(), key=lambda x:x[0])
	sortedDict_Abbr = dict(sortedList_Abbr)
	IE.ExportNestedList(GoodDomainList, OutputFile, Header, Add="_Cutoff-" + str(Cutoff), Ask=False)
	IE.ExportNestedDictionary(AllDomains, OutputFile, Header, Add="_all", Ask=False)
	IE.ExportDictionary(sortedDict_Abbr, OutputFile, Header="Abbr.\tName\n", Add="_Abbreviations", Ask=False)
	return(GoodDomainList, Header)
