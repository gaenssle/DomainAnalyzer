#!/usr/bin/python
# Written in Python 3.7 in 2022 by A.L.O. Gaenssle
# Main file for domain analysis

import os
from multiprocessing import Pool

# Own modules
import Import_Export as IE
import Download_GenomeJP as Genome
import Download_KEGG as KEGG
import Extract_Data as Extract
import Input_Questions as Input
import Count_Data as Count

## ---------------------------------------------------------
## DEFAULT VALUES
## ---------------------------------------------------------

OffsetUniProt = 3 # Offset between domain columns (Name, Start, End) for UniProt
OffsetKEGG = 4 # Offset between domain columns (Name, Start, End) for KEGG
Cutoff = 0.0001 # E-value cutoff for KEGG to get only probable domains


## ----------------------------------------------------------
## HELPFER FUNCTIONS
## ----------------------------------------------------------

# Speedup download by using multiprcessing, WARNING: Has to be installed first!
def MultiProcessing(IDList, Function):
	print("Download data for", len(IDList), "Items. . .")
	if __name__ == '__main__':
		with Pool(10) as pool:
			Import = pool.map(Function, IDList)
		pool.close()
		pool.join()
	return(Import)

# Divide the list of Gene IDs into a nested list of clusters (for multiprcessing)
def GetChunk(List, ClusterSize=250):
	ClusteredList = [List[x:x+ClusterSize] for x in range(0, len(List), ClusterSize)]
	return(ClusteredList)

## ----------------------------------------------------------
## MAIN FUNCTIONS
## ----------------------------------------------------------

# Download all gene IDs associated with the supplied domain name from KEGG, UniProt and PDB
def DownloadList(Domain, OutputFile, DBList=["UniProt", "KEGG", "PDB"]):
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
		if DB == "UniProt":
			Header = "ID\tPfam\tECO\tEMBL\tName\n"
			GeneList, GeneTable = Genome.CleanUniProt(List)
		elif DB == "genes":
			Header = "ID\t#EC\tKO-ID\tName\n"
			GeneList, GeneTable = Genome.CleanKEGG(List)
		elif DB == "PDB":
			Header = "ID\tName\n"
			GeneList, GeneTable = Genome.CleanPDB(List)
		if GeneList != []:
			IE.ExportNestedList(GeneTable, OutputFile, Header, Add=AddToName + "_Table")
			IE.ExportList(GeneList, OutputFile, Add=AddToName + "_List")

# Download details for all given IDs from UniProt, including taxonomy, sequence and domains
def DownloadEntryUniProt(IDList, OutputFile, OutputFolder, OutputFragments, Multiprocess="y"):
	Header = "ID\tKEGG\tOrganism\tTaxonomy\tAA\tSequence\tDomain\tStart\tEnd\n"
	ClusteredList = GetChunk(IDList)
	print("Download motif for", len(IDList), ". . .")
	for ClusterID in range(len(ClusteredList)):
		print("Download cluster", ClusterID+1, "of", len(ClusteredList))
		FragmentFile = OutputFragments + "/" + OutputFile.rsplit(".", 1)[0] + "_" + str(ClusterID+1) + ".txt"
		print(FragmentFile)
		if os.path.exists(FragmentFile):
			print("File already exists, skip to next cluster\n")
		else:
			Table = {}
			if Multiprocess == "y":
				Import = MultiProcessing(ClusteredList[ClusterID],  Genome.DownloadEntryUniProt)
				for Item in Import:
					Table[Item[0]] = Item[1:]
			else:
				Table = {}
				for ID in ClusteredList[ClusterID]:
					Item =  Genome.DownloadEntryUniProt(ID)
					Table[Item[0]] = Item[1:]
			print("Done!\n->", len(Table), "of", len(ClusteredList[ClusterID]), "found")
			IE.ExportNestedDictionary(Table, FragmentFile, Header)
	IE.CombineFiles(OutputFragments, OutputFolder, OutputFile, Header)

# Download details for all given IDs from KEGG, including taxonomy and sequence
def DownloadEntryKEGG(IDList, OutputFile, OutputFolder, OutputFragments, Multiprocess="y"):
	Header = "KEGG\tUniProt\tOrganism\tTaxonomy\tAA\tSequence\n"
	Organisms = KEGG.DownloadOrganismsTemp()
	ClusteredList = GetChunk(IDList, ClusterSize=500)
	print("Download motif for", len(IDList), ". . .")
	for ClusterID in range(len(ClusteredList)):
		print("Download cluster", ClusterID+1, "of", len(ClusteredList))
		FragmentFile = OutputFragments + "/" + OutputFile.rsplit(".", 1)[0] + "_" + str(ClusterID+1) + ".txt"
		ClusterChunk = GetChunk(ClusteredList[ClusterID], ClusterSize=10)
		if os.path.exists(FragmentFile):
			print("File already exists, skip to next cluster\n")
		else:
			Table = {}
			if Multiprocess == "y":
				Import = MultiProcessing(ClusterChunk, KEGG.DownloadProteinEntries)
				for Set in Import:
					for Item in Set:
						Table[Item[0]] = Item[1:]
			else:
				for Chunk in ClusterChunk:
					Set =  KEGG.DownloadProteinEntries(Chunk)
					for Item in Set:
						Table[Item[0]] = Item[1:]
			print("Table downloaded with", len(Table), "entries")
			for Gene in Table:
				try:
					Table[Gene].insert(2, Organisms[Gene.split(":",1)[0]])
				except KeyError:
					Table[Gene].insert(2, "Virus or Unknown")
			IE.ExportNestedDictionary(Table, FragmentFile, Header)
	IE.CombineFiles(OutputFragments, OutputFolder, OutputFile, Header)

# Download all domain motifs for all given IDs from KEGG
def DownloadMotifKEGG(IDList, OutputFile, OutputFolder, OutputFragments, Multiprocess="y"):
	Header = "ID\tOrganism\tDomain\tStart\tEnd\tName\tE-Value\tScore\n"
	ClusterSize = 250
	ClusteredList = GetChunk(IDList)
	print("Download motif for", len(IDList), ". . .")
	for ClusterID in range(len(ClusteredList)):
		print("Download cluster", ClusterID+1, "of", len(ClusteredList))
		FragmentFile = OutputFragments + "/" + OutputFile.rsplit(".", 1)[0] + "_" + str(ClusterID+1) + ".txt"
		print(FragmentFile)
		if os.path.exists(FragmentFile):
			print("File already exists, skip to next cluster")
		else:
			Table = []
			if Multiprocess == "y":
				Import = MultiProcessing(ClusteredList[ClusterID],  KEGG.DownloadMotif)
				for Protein in Import:
					Table.extend(Protein)
			else:
				for ID in ClusteredList[ClusterID]:
					Table.extend(KEGG.DownloadMotif(ID))
			print("Table downloaded with", len(Table), "entries")
			IE.ExportNestedList(Table, FragmentFile, Header)
	IE.CombineFiles(OutputFragments, OutputFolder, OutputFile, Header)

## ----------------------------------------------------------
## ACTUAL SCRIPT
## ----------------------------------------------------------

# Questionaire to get User Input
Input.PrintHeader()
DirectoryName, ScriptName = os.path.split(os.path.abspath(__file__))
os.chdir(DirectoryName)
Folder, Ask = Input.GetFolder(DirectoryName)
Name = Input.GetName(Folder)
DBList, Action = Input.GetAction()
SpeedUp = Input.UseMultiprocess(Action)

# Download Sequence IDs from UniProt, KEGG and/or PDB
if "i" in Action:
	DownloadList(Name, Folder + "/Input/" + Name + ".txt", DBList)

# download protein data from UniProt and/or KEGG
if "d" in Action:
	for DB in DBList:
		InputFile = IDList = Folder + "/Input/" + Name + "_" + DB + "_List.txt"
		if not os.path.exists(InputFile):
			DownloadList(Name, Folder + "/Input/" + Name + ".txt", [DB])
		try:
			IDList = IE.ImportList(InputFile)
			OutputFile = Name + "_" + DB + "_Protein.txt"
			OutputFragments = IE.CreateFolder(Folder + "/Output/" + Name + "_" + DB + "_ProteinFragments")
			if DB == "UniProt":
				DownloadEntryUniProt(IDList, OutputFile, Folder + "/Output/", OutputFragments, Multiprocess=SpeedUp)
			elif DB == "KEGG":
				DownloadEntryKEGG(IDList, OutputFile, Folder + "/Output/", OutputFragments, Multiprocess=SpeedUp)
		except FileNotFoundError:
			print(DB, "does not contain any items for domain", Name)


# Download motif data from KEGG
if "m" in Action:
	if "KEGG" in DBList:
		InputFile = IDList = Folder + "/Input/" + Name + "_KEGG_List.txt"
		if not os.path.exists(InputFile):
			DownloadList(Name, Folder + "/Input/" + Name + ".txt", ["KEGG"])
		IDList = IE.ImportList(InputFile)
		OutputFile = Name + "_KEGG_Motif.txt"
		OutputFragments = IE.CreateFolder(Folder + "/Output/" + Name + "_KEGG_MotifFragments")
		DownloadMotifKEGG(IDList, OutputFile, Folder + "/Output/", OutputFragments, Multiprocess=SpeedUp)

# Extract data from UniProt and/or KEGG
if "e" in Action:
	for DB in DBList:
		InputFile = Folder + "/Output/" + Name + "_" + DB + "_Protein_all.txt"
		OutputFile = Folder + "/Output/" + Name + "_" + DB + "_Domains.txt"
		if DB == "UniProt":
			GeneTable, Header = IE.ImportNestedList(InputFile, getHeader=True)
			Motifs = Extract.ExtractMotifs(GeneTable, Header, Name, OffsetUniProt, OutputFile, Ask=Ask)
			DetailsOnly = Extract.ExtractDetails(GeneTable, Header, Name, OffsetUniProt, OutputFile, Ask=Ask)
			Extract.CreateFasta(DetailsOnly, OutputFile, Ask=Ask)
		elif DB == "KEGG":
			MotifFile = Folder + "/Output/" + Name + "_" + DB + "_Motif_all.txt"
			GeneTable, Header = IE.ImportNestedDictionary(InputFile, getHeader=True)
			MotifTable = IE.ImportNestedList(MotifFile)
			GoodDomains, Header = Extract.AddMotifKEGG(GeneTable, MotifTable, Header, Name, Cutoff, OutputFile, Ask=Ask)
			Motifs = Extract.ExtractMotifs(GoodDomains, Header, Name, OffsetKEGG, OutputFile, Ask=Ask)
			DetailsOnly = Extract.ExtractDetails(GoodDomains, Header, Name, OffsetKEGG, OutputFile, Ask=Ask)
			Extract.CreateFasta(DetailsOnly, OutputFile, Ask=Ask)

# Count occurence of same taxonomy and domain motifs
if "c" in Action:
	for DB in DBList:
		try:
			InputFile = Folder + "/Output/" + Name + "_" + DB + "_Domains_Motifs.txt"
			GeneTable, Header = IE.ImportNestedList(InputFile, getHeader=True)
			Count.CountMotif(GeneTable, Folder, Name, DB, Ask=Ask)
			Count.CountTaxonomy(GeneTable, Header, Folder, Name, DB, Ask=Ask)
		except FileNotFoundError:
			pass

Input.PrintFooter()
