#!/usr/bin/python
# Written in Python 3.10 in 2023 by A.L.O. Gaenssle

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


## ------------------------------------------------------------------------------------------------
## DEFAULT VALUES ---------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------

# OffsetUniProt = 3 # Offset between domain columns (Name, Start, End) for UniProt
# OffsetKEGG = 4 # Offset between domain columns (Name, Start, End) for KEGG
Cutoff = 0.0001 # E-value cutoff for KEGG to get only probable domains (default=0.0001)
FileType = ".csv" # File types of all exported files (default="csv")
Sep = ";" # Separator between the columns (default=";")
Ask = False # Ask if files should be replaced (default=True)
ClusterSize = 20 # Define the number of entires in which the download is saved (default=250)

## ------------------------------------------------------------------------------------------------
## HELPFER FUNCTIONS ------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------
## ================================================================================================
## Speedup download by using multiprcessing, WARNING: Has to be installed first!
def MultiProcessing(IDList, Function):
	Import = []
	print("Download data for", len(IDList), "Items. . .")
	if __name__ == '__main__':
		with Pool(10) as pool:
			Import = pool.map(Function, IDList)
		pool.close()
		pool.join()
	return(Import)

## ------------------------------------------------------------------------------------------------
## MAIN FUNCTIONS ---------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------
## ================================================================================================
## Download all gene IDs associated with the supplied domain name from KEGG, UniProt and PDB
def DownloadList(Name, OutputFile, DB, FileType, Sep, Ask):
	print("\nNow downloading gene IDs containing domain", Name, ". . .\n(Speed depends on internet connection)\n")

	# Create the required url with the following fragments
	urlGeneList = "https://www.genome.jp/dbget-bin/get_linkdb?-t+"
	urlName = "+pf:"
	urlPage = "+-p+"
	urlKEGG = "genes"
	AddToName = "_" + DB
	DB = DB.replace("KEGG", urlKEGG)
	urlInitial = urlGeneList + DB.lower() + urlName + Name

	# Download all genes from the first page and then cycle through all subsequent pages
	List, Pages = Genome.DownloadGeneList(urlInitial, getAmount=True)
	for Page in range(2, Pages+1):
		urlNew = urlGeneList + DB.lower() + urlPage + str(Page) + urlName + Name
		List.extend(Genome.DownloadGeneList(urlNew))

	# Extract all information from the entries and convert into pandas dataframe
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

## ================================================================================================
## Download details for all given IDs from UniProt, including taxonomy, sequence and names
def DownloadEntryUniProt(IDList, FilePath, FileType, Sep, Multiprocess, ClusterSize, Ask):
	print("Download protein data for", len(IDList), ". . .")

	# Create clusters of sequences to generate smaller files (in case the download crashes)
	ClusteredList = [IDList[x:x+ClusterSize] for x in range(0, len(IDList), ClusterSize)]
	for ClusterID in range(len(ClusteredList)):
		print("Download cluster", ClusterID+1, "of", len(ClusteredList))
		FragmentFile = FilePath + "_" + str(ClusterID+1)
		print(FragmentFile)

		# Ignore all files that have already been downloaded
		if os.path.exists(FragmentFile + FileType):
			print("File already exists, skip to next cluster\n")

		# Download all files that have not yet been saved
		else:
			if Multiprocess == "y":
				ListOfDicts = MultiProcessing(ClusteredList[ClusterID],  Genome.DownloadEntryUniProt)
			else:
				ListOfDicts = []
				for ID in ClusteredList[ClusterID]:
					Dict =  Genome.DownloadEntryUniProt(ID)
					ListOfDicts.append(Dict)

			# Convert to pandas dataframe and store each fragment in the fragment folder
			ProteinTable = pd.DataFrame(ListOfDicts)
			print("Done!\n->", len(ProteinTable), "of", len(ClusteredList[ClusterID]), "found")
			IE.ExportDataFrame(ProteinTable, FragmentFile, FileType=FileType, Sep=Sep, Ask=Ask)

	# After all entries have been downloaded, combine all fragments into one dataframe
	DataFrame = IE.CombineFiles(os.path.split(FragmentFile)[0], Sep)
	return(DataFrame)

## ================================================================================================
## Download details for all given IDs from KEGG, including taxonomy and sequence
def DownloadEntryKEGG(IDList, FilePath, FileType, Sep, Multiprocess, ClusterSize, Ask):
	Organisms = []
	print("Download protein data for", len(IDList), ". . .")

	# Create clusters of sequences to generate smaller files (in case the download crashes)
	ClusteredList = [IDList[x:x+ClusterSize] for x in range(0, len(IDList), ClusterSize)]
	for ClusterID in range(len(ClusteredList)):
		print("Download cluster", ClusterID+1, "of", len(ClusteredList))
		FragmentFile = FilePath + "_" + str(ClusterID+1)
		print(FragmentFile)

		# Create chunks of clusters since the data of 10 proteins can be downloaded from KEGG at once
		ClusterChunk = [ClusteredList[ClusterID][x:x+10] for x in range(0, len(ClusteredList[ClusterID]), 10)]

		# Ignore all files that have already been downloaded
		if os.path.exists(FragmentFile + FileType):
			print("File already exists, skip to next cluster\n")
		else:
			# Multiprocessing saves the entries as a list of list of dictionaries
			if Multiprocess == "y":
				ClusteredListOfDicts = MultiProcessing(ClusterChunk, KEGG.DownloadProteinEntries)
				ListOfDicts = [Entry for Cluster in ClusteredListOfDicts for Entry in Cluster]
			else:
				ListOfDicts = []
				for Chunk in ClusterChunk:
					Set =  KEGG.DownloadProteinEntries(Chunk)
					ListOfDicts.extend(Set)

			# Only download the list of organisms on KEGG if needed and add to dataframe
			if len(Organisms) == 0:
				Organisms = KEGG.DownloadOrganismsTemp()
			ProteinTable = pd.DataFrame(ListOfDicts)
			ProteinTable = pd.merge(ProteinTable, Organisms, on=["orgID"],  how="left")
			IE.ExportDataFrame(ProteinTable, FragmentFile, FileType=FileType, Sep=Sep, Ask=Ask)
			print("Done!\n->", len(ProteinTable), "of", len(ClusteredList[ClusterID]), "found")

	# After all entries have been downloaded, combine all fragments into one dataframe
	DataFrame = IE.CombineFiles(os.path.split(FragmentFile)[0], Sep)
	return(DataFrame)

## ================================================================================================
## Download all domain motifs for all given IDs from KEGG
def DownloadMotifKEGG(IDList, FilePath, FileType, Sep, Multiprocess, ClusterSize, Ask):
	print("Download motif for", len(IDList), ". . .")

	# Create clusters of sequences to generate smaller files (in case the download crashes)
	ClusteredList = [IDList[x:x+ClusterSize] for x in range(0, len(IDList), ClusterSize)]
	for ClusterID in range(len(ClusteredList)):
		print("Download cluster", ClusterID+1, "of", len(ClusteredList))
		FragmentFile = FilePath + "_" + str(ClusterID+1)
		print(FragmentFile)

		# Ignore all files that have already been downloaded
		if os.path.exists(FragmentFile):
			print("File already exists, skip to next cluster")
		else:
			ListOfLists = []
			if Multiprocess == "y":
				ClusteredListOfLists = MultiProcessing(ClusteredList[ClusterID],  KEGG.DownloadMotif)
				ListOfLists = [Entry for Cluster in ClusteredListOfLists for Entry in Cluster]
			else:
				for ID in ClusteredList[ClusterID]:
					ListOfLists.extend(KEGG.DownloadMotif(ID))
			ColNames= ["ID", "Index","Domain", "Start", "End", "Definition", "E-Value", "Score"]
			MotifTable = pd.DataFrame(ListOfLists, columns=ColNames)
			IE.ExportDataFrame(MotifTable, FragmentFile, FileType=FileType, Sep=Sep, Ask=Ask)
			print("Done!\n->", len(MotifTable), "of", len(ClusteredList[ClusterID]), "found")

	# After all entries have been downloaded, combine all fragments into one dataframe
	DataFrame = IE.CombineFiles(os.path.split(FragmentFile)[0], Sep)
	return(DataFrame)
	# 		print("Table downloaded with", len(Table), "entries")
	# 		IE.ExportNestedList(Table, FragmentFile, Header)
	# IE.CombineFiles(OutputFragments, OutputFolder, OutputFile, Header)

## ------------------------------------------------------------------------------------------------
## SCRIPT -----------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------

# Skip manual imput section
Multiprocess = "y"
Folder = "DUF5727"
# Folder = "DUF1735"
Name = "DUF5727"
# Name = "DUF1735"
# DBList = ["KEGG"]
# DBList = ["UniProt"]
DBList = ["UniProt", "KEGG"]
# DBList = ["UniProt", "KEGG", "PDB"]
# DBList = ["PDB"]
Action = "m"


IE.CreateFolder(Folder + "/Input")
IE.CreateFolder(Folder + "/Output")

for DB in DBList:
	FileName = Name + "_" + DB
	# Download Sequence IDs from UniProt, KEGG and/or PDB
	if any(s in ["i", "d", "m"] for s in Action):
		InputFile = os.path.join(Folder, "Input", FileName + FileType)
		if not os.path.exists(InputFile):
			OutputFile = os.path.join(Folder, "Input", Name)
			DownloadList(Name, OutputFile, DB, FileType, Sep, Ask)

	# download protein data from UniProt and/or KEGG
	if "d" in Action:
		try:
			DataFrame = pd.read_csv(InputFile, sep=Sep)
			IDList = DataFrame["ID"][:50].tolist()
			OutputPath = os.path.join(Folder, "Output", FileName + "_Protein")
			FragmentFolder = IE.CreateFolder(OutputPath + "Fragments")
			FragmentFile = os.path.join(FragmentFolder, FileName + "_Protein")
			if DB == "UniProt":
				Detailed = DownloadEntryUniProt(IDList, FragmentFile, FileType, Sep, Multiprocess, ClusterSize, Ask)
			elif DB == "KEGG":
				Detailed = DownloadEntryKEGG(IDList, FragmentFile, FileType, Sep, Multiprocess, ClusterSize, Ask)
		except FileNotFoundError:
			print(DB, "does not contain any items for domain", Name)
		DataFrame = pd.merge(DataFrame, Detailed, on=["ID"],  how="outer")
		IE.ExportDataFrame(DataFrame, OutputPath, Ask=Ask)

	# Download motif data from KEGG
	if "m" in Action and DB == "KEGG":
		DataFrame = pd.read_csv(InputFile, sep=Sep)
		IDList = DataFrame["ID"][:50].tolist()
		OutputPath = os.path.join(Folder, "Output", FileName + "_Motif")
		FragmentFolder = IE.CreateFolder(OutputPath + "Fragments")
		FragmentFile = os.path.join(FragmentFolder, FileName + "_Motif")
		# print(IDList)
		# OutputFile = Name + "_KEGG_Motif.txt"
		# print(OutputFile)
		# OutputFragments = IE.CreateFolder(Folder + "/Output/" + Name + "_KEGG_MotifFragments")
		DataFrame = DownloadMotifKEGG(IDList, FragmentFile, FileType, Sep, Multiprocess, ClusterSize, Ask)
		IE.ExportDataFrame(DataFrame, OutputPath, Ask=Ask)
