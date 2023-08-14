#!/usr/bin/python
# Written in Python 3.7 in 2022 by A.L.O. Gaenssle

# import os
# import re
# import math
# import urllib.request
# import ssl
# ssl._create_default_https_context = ssl._create_unverified_context
# import pandas as pd

# Own modules
import Import_Export as IE
import Download_GenomeJP as Genome



# Download all gene IDs associated with the supplied domain name from KEGG, UniProt and PDB
def DownloadList(Domain, OutputFile, DBList, Ask):
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
				IE.ExportDataFrame(GeneTable, OutputFile, Add=AddToName + "_Table", Ask=Ask)
				IE.ExportDataFrame(GeneTable, OutputFile, Columns=["Gene ID"], Header=False, Add=AddToName + "_List", Ask=Ask)




# Skip manual imput section
Folder = "DUF5727"
# Folder = "DUF1735"
Ask = False
Name = "DUF5727"
# Name = "DUF1735"
# DBList = ["KEGG"]
DBList = ["UniProt", "KEGG", "PDB"]
# DBList = ["PDB"]
Action = "i"

IE.CreateFolder(Folder + "/Input")
IE.CreateFolder(Folder + "/Output")

# Download Sequence IDs from UniProt, KEGG and/or PDB
if "i" in Action:
	DownloadList(Name, Folder + "/Input/" + Name + ".txt", DBList, Ask=Ask)

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