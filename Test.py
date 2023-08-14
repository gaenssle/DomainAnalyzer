#!/usr/bin/python
# Written in Python 3.7 in 2022 by A.L.O. Gaenssle

import os
from multiprocessing import Pool
# import re
# import math
# import urllib.request
# import ssl
# ssl._create_default_https_context = ssl._create_unverified_context
import pandas as pd

# Own modules
import Import_Export as IE
import Download_GenomeJP as Genome
import Download_KEGG as KEGG
# import Extract_Data as Extract
# import Input_Questions as Input
# import Count_Data as Count


## ---------------------------------------------------------------------------------------------------
## DEFAULT VALUES
## ---------------------------------------------------------------------------------------------------

OffsetUniProt = 3 # Offset between domain columns (Name, Start, End) for UniProt
OffsetKEGG = 4 # Offset between domain columns (Name, Start, End) for KEGG
Cutoff = 0.0001 # E-value cutoff for KEGG to get only probable domains


## ----------------------------------------------------------------------------------------------------
## HELPFER FUNCTIONS
## ----------------------------------------------------------------------------------------------------

# Speedup download by using multiprcessing, WARNING: Has to be installed first!
def MultiProcessing(IDList, Function):
	Import = []
	print("Download data for", len(IDList), "Items. . .")
	print("TEST TEST TEST TEST")
	if __name__ == '__main__':
		print("IF LOOOOOOOOOOOP")
		with Pool(10) as pool:
			print("POOOOOOOLLLLLLL")
			Import = pool.map(Function, IDList)
		pool.close()
		pool.join()
	return(Import)

# Divide the list of Gene IDs into a nested list of clusters (for multiprcessing)
def GetChunk(List, ClusterSize=250):
	ClusteredList = [List[x:x+ClusterSize] for x in range(0, len(List), ClusterSize)]
	return(ClusteredList)


## ----------------------------------------------------------------------------------------------------
## MAIN FUNCTIONS
## ----------------------------------------------------------------------------------------------------

# Download all gene IDs associated with the supplied domain name from KEGG, UniProt and PDB
def DownloadList(Domain, OutputFile, DBList, Ask=True):
	print("\nNow downloading gene IDs containing domain", Domain, ". . .\n(Speed depends on internet connection)\n")
	urlGeneList = "https://www.genome.jp/dbget-bin/get_linkdb?-t+"
	urlDomain = "+pf:"
	urlPage = "+-p+"
	urlKEGG = "genes"
	for DB in DBList:
		AddToName = "_" + DB
		DB = DB.replace("KEGG", urlKEGG)
		List, Pages = Genome.DownloadGeneList(urlGeneList + DB.lower() + urlDomain + Domain, getAmount=True)
		for Page in range(2, Pages+1):
			urlNew = urlGeneList + DB.lower() + urlPage + str(Page) + urlDomain + Domain
			List.extend(Genome.DownloadGeneList(urlNew))
		print("List contains", len(List), "items\n")
		if List != []:
			if DB == "UniProt":
				GeneTable = Genome.CleanUniProt(List)
			elif DB == "genes":
				GeneTable = Genome.CleanKEGG(List)
			elif DB == "PDB":
				GeneTable = Genome.CleanPDB(List)
			if len(GeneTable.index) != 0:
				IE.ExportDataFrame(GeneTable, OutputFile, Add=AddToName, Ask=Ask)

# Download details for all given IDs from UniProt, including taxonomy, sequence and domains
def DownloadEntryUniProt(IDList, OutputFile, OutputFolder, OutputFragments, Multiprocess="y", ClusterSize=250, Ask=False):
	ClusteredList = [IDList[x:x+ClusterSize] for x in range(0, len(IDList), ClusterSize)]
	print("Download protein data for", len(IDList), ". . .")
	for ClusterID in range(len(ClusteredList)):
		print("Download cluster", ClusterID+1, "of", len(ClusteredList))
		FragmentFile = OutputFragments + "/" + OutputFile
		print(FragmentFile)
		if os.path.exists(FragmentFile):
			print("File already exists, skip to next cluster\n")
		else:
			if Multiprocess == "y":
				ListOfDicts = MultiProcessing(ClusteredList[ClusterID],  Genome.DownloadEntryUniProt)
				# for Item in Import:
				# 	Dict =  Genome.DownloadEntryUniProt(ID)
				# 	ListOfDicts.append(Dict)
			else:
				ListOfDicts = []
				for ID in ClusteredList[ClusterID]:
					Dict =  Genome.DownloadEntryUniProt(ID)
					ListOfDicts.append(Dict)
			ProteinTable = pd.DataFrame(ListOfDicts)
			# FirstColumns = ["Protein ID", "Organism", "Taxonomy"]
			# ColumnOrder = FirstColumns + (ProteinTable.columns.drop(FirstColumns).tolist())
			# ProteinTable = ProteinTable[ColumnOrder]


					# Table[Item[0]] = Item[1:]
			print("Done!\n->", len(ProteinTable), "of", len(ClusteredList[ClusterID]), "found")
			IE.ExportDataFrame(ProteinTable, FragmentFile, Add="_" + str(ClusterID+1), Ask=Ask)
			# IE.ExportNestedDictionary(Table, FragmentFile, Header)
	# IE.CombineFiles(OutputFragments, OutputFolder, OutputFile, Header)

# # Download details for all given IDs from KEGG, including taxonomy and sequence
# def DownloadEntryKEGG(IDList, OutputFile, OutputFolder, OutputFragments, Multiprocess="y"):
# 	Header = "KEGG\tUniProt\tOrganism\tTaxonomy\tAA\tSequence\n"
# 	Organisms = KEGG.DownloadOrganismsTemp()
# 	ClusteredList = GetChunk(IDList, ClusterSize=500)
# 	print("Download motif for", len(IDList), ". . .")
# 	for ClusterID in range(len(ClusteredList)):
# 		print("Download cluster", ClusterID+1, "of", len(ClusteredList))
# 		FragmentFile = OutputFragments + "/" + OutputFile.rsplit(".", 1)[0] + "_" + str(ClusterID+1) + ".txt"
# 		ClusterChunk = GetChunk(ClusteredList[ClusterID], ClusterSize=10)
# 		if os.path.exists(FragmentFile):
# 			print("File already exists, skip to next cluster\n")
# 		else:
# 			Table = {}
# 			if Multiprocess == "y":
# 				Import = MultiProcessing(ClusterChunk, KEGG.DownloadProteinEntries)
# 				for Set in Import:
# 					for Item in Set:
# 						Table[Item[0]] = Item[1:]
# 			else:
# 				for Chunk in ClusterChunk:
# 					Set =  KEGG.DownloadProteinEntries(Chunk)
# 					for Item in Set:
# 						Table[Item[0]] = Item[1:]
# 			print("Table downloaded with", len(Table), "entries")
# 			for Gene in Table:
# 				try:
# 					Table[Gene].insert(2, Organisms[Gene.split(":",1)[0]])
# 				except KeyError:
# 					Table[Gene].insert(2, "Virus or Unknown")
# 			IE.ExportNestedDictionary(Table, FragmentFile, Header)
# 	IE.CombineFiles(OutputFragments, OutputFolder, OutputFile, Header)


## ----------------------------------------------------------------------------------------------------
## SCRIPT
## ----------------------------------------------------------------------------------------------------

# Skip manual imput section
SpeedUp = "n"
Folder = "DUF5727"
# Folder = "DUF1735"
Ask = False
Name = "DUF5727"
# Name = "DUF1735"
# DBList = ["KEGG"]
DBList = ["UniProt"]
# DBList = ["UniProt", "KEGG", "PDB"]
# DBList = ["PDB"]
Action = "d"

Cycle = 0
print("CYCLE NUMBER" ,Cycle)
Cycle += 1
IE.CreateFolder(Folder + "/Input")
IE.CreateFolder(Folder + "/Output")

# Download Sequence IDs from UniProt, KEGG and/or PDB
if "i" in Action:
	DownloadList(Name, Folder + "/Input/" + Name + ".txt", DBList, Ask=Ask)

# download protein data from UniProt and/or KEGG
if "d" in Action:
	for DB in DBList:
		# InputFile = Folder + "/Input/" + Name + "_" + DB + "_List.txt"
		InputFile = Folder + "/Input/" + Name + "_" + DB + ".txt"
		if not os.path.exists(InputFile):
			DownloadList(Name, Folder + "/Input/" + Name + ".txt", [DB])
		try:
			# FileName, Header=0, UseCols=[], Sep="\t", Stamp=False
			DataFrame = IE.ImportDataFrame(InputFile, UseCols = ["Gene ID"])
			IDList = DataFrame["Gene ID"][:50].tolist()
			# print(IDList)
			OutputFile = Name + "_" + DB + "_Protein.txt"
			OutputFragments = IE.CreateFolder(Folder + "/Output/" + Name + "_" + DB + "_ProteinFragments")
			if DB == "UniProt":
				DownloadEntryUniProt(IDList, OutputFile, Folder + "/Output/", OutputFragments, Multiprocess=SpeedUp)
			elif DB == "KEGG":
				DownloadEntryKEGG(IDList, OutputFile, Folder + "/Output/", OutputFragments, Multiprocess=SpeedUp)
		except FileNotFoundError:
			print(DB, "does not contain any items for domain", Name)