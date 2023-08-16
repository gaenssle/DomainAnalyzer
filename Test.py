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

# OffsetUniProt = 3 # Offset between domain columns (Name, Start, End) for UniProt
# OffsetKEGG = 4 # Offset between domain columns (Name, Start, End) for KEGG
Cutoff = 0.0001 # E-value cutoff for KEGG to get only probable domains (default=0.0001)
FileType = ".csv" # File types of all exported files (default="csv")
Sep = ";" # Separator between the columns (default=";")
Ask = False # Ask if files should be replaced (default=True)
ClusterSize = 20 # Define the number of entires in which the download is saved (default=250)

## ----------------------------------------------------------------------------------------------------
## HELPFER FUNCTIONS
## ----------------------------------------------------------------------------------------------------

# Speedup download by using multiprcessing, WARNING: Has to be installed first!
def MultiProcessing(IDList, Function):
	Import = []
	print("Download data for", len(IDList), "Items. . .")
	if __name__ == '__main__':
		with Pool(10) as pool:
			Import = pool.map(Function, IDList)
		pool.close()
		pool.join()
	return(Import)

## ----------------------------------------------------------------------------------------------------
## MAIN FUNCTIONS
## ----------------------------------------------------------------------------------------------------

# Download all gene IDs associated with the supplied domain name from KEGG, UniProt and PDB
def DownloadList(Name, OutputFile, DBList, FileType, Sep, Ask):
	print("\nNow downloading gene IDs containing domain", Name, ". . .\n(Speed depends on internet connection)\n")
	urlGeneList = "https://www.genome.jp/dbget-bin/get_linkdb?-t+"
	urlName = "+pf:"
	urlPage = "+-p+"
	urlKEGG = "genes"
	for DB in DBList:
		AddToName = "_" + DB
		DB = DB.replace("KEGG", urlKEGG)
		urlInitial = urlGeneList + DB.lower() + urlName + Name
		List, Pages = Genome.DownloadGeneList(urlInitial, getAmount=True)
		for Page in range(2, Pages+1):
			urlNew = urlGeneList + DB.lower() + urlPage + str(Page) + urlName + Name
			List.extend(Genome.DownloadGeneList(urlNew))
		if List != []:
			if DB == "UniProt":
				GeneTable = Genome.CleanUniProt(List)
			elif DB == "genes":
				GeneTable = Genome.CleanKEGG(List)
			elif DB == "PDB":
				GeneTable = Genome.CleanPDB(List)
			if len(GeneTable.index) != 0:
				print("Data found for", len(GeneTable.index), "entries\n")
				IE.ExportDataFrame(GeneTable, OutputFile, Add=AddToName, FileType=FileType, Sep=Sep, Ask=Ask)


# Download details for all given IDs from UniProt, including taxonomy, sequence and names
def DownloadEntryUniProt(IDList, FilePath, FileType, Sep, Multiprocess, ClusterSize, Ask):
	print("Download protein data for", len(IDList), ". . .")
	ClusteredList = [IDList[x:x+ClusterSize] for x in range(0, len(IDList), ClusterSize)]
	for ClusterID in range(len(ClusteredList)):
		print("Download cluster", ClusterID+1, "of", len(ClusteredList))
		FragmentFile = FilePath + "_" + str(ClusterID+1)
		print(FragmentFile)
		if os.path.exists(FragmentFile + FileType):
			print("File already exists, skip to next cluster\n")
		else:
			if Multiprocess == "y":
				ListOfDicts = MultiProcessing(ClusteredList[ClusterID],  Genome.DownloadEntryUniProt)
			else:
				ListOfDicts = []
				for ID in ClusteredList[ClusterID]:
					Dict =  Genome.DownloadEntryUniProt(ID)
					ListOfDicts.append(Dict)
			ProteinTable = pd.DataFrame(ListOfDicts)
			print("Done!\n->", len(ProteinTable), "of", len(ClusteredList[ClusterID]), "found")
			IE.ExportDataFrame(ProteinTable, FragmentFile, FileType=FileType, Sep=Sep, Ask=Ask)
	DataFrame = IE.CombineFiles(os.path.split(FragmentFile)[0], Sep)
	return(DataFrame)

# Download details for all given IDs from KEGG, including taxonomy and sequence
def DownloadEntryKEGG(IDList, FilePath, FileType, Sep, Multiprocess, ClusterSize, Ask):
	Organisms = KEGG.DownloadOrganismsTemp()
	print("Download protein data for", len(IDList), ". . .")
	ClusteredList = [IDList[x:x+ClusterSize] for x in range(0, len(IDList), ClusterSize)]
	for ClusterID in range(len(ClusteredList)):
		print("Download cluster", ClusterID+1, "of", len(ClusteredList))
		FragmentFile = FilePath + "_" + str(ClusterID+1)
		print(FragmentFile)
		ClusterChunk = [ClusteredList[ClusterID][x:x+10] for x in range(0, len(ClusteredList[ClusterID]), 10)]
		if os.path.exists(FragmentFile + FileType):
			print("File already exists, skip to next cluster\n")
		else:
			if Multiprocess == "y":
				ClusteredListOfDicts = MultiProcessing(ClusterChunk, KEGG.DownloadProteinEntries)
				ListOfDicts = [Entry for Cluster in ClusteredListOfDicts for Entry in Cluster]
			else:
				ListOfDicts = []
				for Chunk in ClusterChunk:
					Set =  KEGG.DownloadProteinEntries(Chunk)
					ListOfDicts.extend(Set)
			ProteinTable = pd.DataFrame(ListOfDicts)
			ProteinTable = pd.merge(ProteinTable, Organisms, on=["orgID"],  how="left")
			IE.ExportDataFrame(ProteinTable, FragmentFile, FileType=FileType, Sep=Sep, Ask=Ask)
			print("Done!\n->", len(ProteinTable), "of", len(ClusteredList[ClusterID]), "found")
	DataFrame = IE.CombineFiles(os.path.split(FragmentFile)[0], Sep)
	return(DataFrame)


## ----------------------------------------------------------------------------------------------------
## SCRIPT
## ----------------------------------------------------------------------------------------------------

# Skip manual imput section
# SpeedUp = "y"
Multiprocess = "y"
Folder = "DUF5727"
# Folder = "DUF1735"
Name = "DUF5727"
# Name = "DUF1735"
DBList = ["KEGG"]
# DBList = ["UniProt"]
# DBList = ["UniProt", "KEGG", "PDB"]
# DBList = ["PDB"]
Action = "d"


IE.CreateFolder(Folder + "/Input")
IE.CreateFolder(Folder + "/Output")

# Download Sequence IDs from UniProt, KEGG and/or PDB
if "i" in Action:
	OutputFile = os.path.join(Folder, "Input", Name)
	DownloadList(Name, OutputFile, DBList, FileType, Sep, Ask)

# download protein data from UniProt and/or KEGG
if "d" in Action:
	for DB in DBList:
		InputFile = os.path.join(Folder, "Input", Name + "_" + DB + FileType)
		if not os.path.exists(InputFile):
			OutputFile = os.path.join(Folder, "Input", Name)
			DownloadList(Name, OutputFile, [DB], FileType, Sep, Ask)
		try:
			DataFrame = pd.read_csv(InputFile, sep=Sep)
			IDList = DataFrame["ID"][:50].tolist()
			OutputPath = os.path.join(Folder, "Output", Name + "_" + DB + "_Protein")
			FragmentFolder = IE.CreateFolder(OutputPath + "Fragments")
			FragmentFile = os.path.join(FragmentFolder, Name + "_" + DB + "_Protein")
			if DB == "UniProt":
				Detailed = DownloadEntryUniProt(IDList, FragmentFile, FileType, Sep, Multiprocess, ClusterSize, Ask)
			elif DB == "KEGG":
				Detailed = DownloadEntryKEGG(IDList, FragmentFile, FileType, Sep, Multiprocess, ClusterSize, Ask)
		except FileNotFoundError:
			print(DB, "does not contain any items for domain", Name)
		DataFrame = pd.merge(DataFrame, Detailed, on=["ID"],  how="outer")
		IE.ExportDataFrame(DataFrame, OutputPath, Ask=Ask)
