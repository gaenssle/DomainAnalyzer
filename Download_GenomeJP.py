#!/usr/bin/python
# Written in Python 3.7 in 2022 by A.L.O. Gaenssle
# Downloads gene IDs from Genome.jp

import re
import math
import urllib.request
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

##------------------------------------------------------
## DOWNLOAD FUNCTIONS
##------------------------------------------------------

# Download gene list from Genome.jp
def DownloadGeneList(url, getAmount=False):
	Amount = 0
	with urllib.request.urlopen(url) as File:
		Add = False
		List = []
		for Line in File:
			Line = Line.decode("utf-8").strip()
			Line = re.sub('<[^>]*>', '', Line)
			if Line != "":
				if Add == False:
					if Line.startswith("Hits:"):
						Amount = Line.split(" ",2)[1]
						Amount = math.ceil(int(Amount)/1000)
					elif Line.startswith("----------"):
						Add = True
				else:
					if Line.startswith("DBGET integrated"):
						break
					else:
						List.append(Line.strip())
	print(url, "downloaded")
	print(len(List), "gene IDs added")
	if getAmount:
		return(List, Amount)
	else:
		return(List)


# Download UniProt protein entry
def DownloadEntryUniProt(ID):
	url = "https://www.genome.jp/entry/up:" + ID
	with urllib.request.urlopen(url) as File:
		InSequence = False
		InDomain = False
		KEGG = ""
		Organism = ""
		TaxString = ""
		Taxonomy = ""
		Length = ""
		Domains = []
		Sequence = ""
		for Line in File:
			Line = Line.decode("utf-8").strip()
			Line = re.sub('<[^>]*>', '', Line)
			if Line != "":
				try:
					if Line.startswith("ID"):
						Length = Line.rsplit(";",1)[1].rsplit(" ",1)[0].strip()
					elif Line.startswith("OS"):
						Organism = Line.split(" ",1)[1].replace(".", "").strip()
					elif Line.startswith("OC"):
						TaxString += Line.split(" ",1)[1].replace(".", "").strip()
					elif Line.startswith("DR"):
						if "KEGG" in Line:
							KEGG = Line.split(";")[1].strip()
					elif Line.startswith("FT"):
						if "DOMAIN" in Line:
							InDomain = True
							String = Line.rsplit(" ",1)[1]
							Domains.extend(String.split(".."))
						elif InDomain and "/note=" in Line and (len(Domains)+1) % 3 == 0:
							Domains.insert(-2, Line.split("\"")[1])
					elif InSequence:
						if Line.startswith("//"):
							break
						else:
							Sequence += Line.replace(" ", "")
					elif Line.startswith("SQ"):
						InSequence = True
				except IndexError:
					pass
	try:
		Kingdom, Phylum = TaxString.split("; ",3)[:2]
		Taxonomy = Kingdom + "-" + Phylum
	except ValueError:
		Taxonomy = "-".join(TaxString.split("; ",))
	print("ID", ID, "downloaded")
	Details = [ID, KEGG, Organism, Taxonomy, Length, Sequence] + Domains
	return(Details)


##------------------------------------------------------
## CLEANUP DATA FUNCTIONS
##------------------------------------------------------

# Convert downloaded UniProt text to table
def CleanUniProt(Data):
	GeneTable = []
	GeneList = []
	for Line in Data:
		try:
			GeneID, String = Line.split(" ",1)
			String = String.strip().split("Full=",1)[1]
			Name, IDString = String.split("{",1)
			IDList = IDString.split("};",1)[0].split("|")
			ECO = ""
			Pfam = ""
			EMBL = ""
			for ID in IDList:
				if len(ID) > 10:
					if ID.startswith("ECO:"):
						ECO = ID.split(":",1)[1]
					elif ID.startswith("EMBL:"):
						EMBL = ID.split(":",1)[1]
					elif ID.startswith("Pfam:"):
						Pfam = ID.split(":",1)[1]
			GeneTable.append([GeneID, Pfam, ECO, EMBL, Name])
			GeneList.append(GeneID)
		except ValueError:
			pass
	return(GeneList, GeneTable)

# Convert downloaded KEGG gene text to table
def CleanKEGG(Data):
	GeneTable = []
	GeneList = []
	for Line in Data:
		try:
			GeneID, String = Line.split(" ",1)
			String = String.strip()
			if "no KO assigned" in String:
				if "Bank" in String:
					Name = String.split("Bank) ")[1]
				elif "RefSeq" in String:
					Name = String.split("RefSeq)")[1]
				else:
					Name = String.split(")")[1]
				KOID = ""
			else:
				if "|" in String:
					String = String.split("|")[0]
					KOID, Name = String.split(" ", 1)
			if "[EC" in Name:
				Name, EC = Name.split(" [EC:")
				EC = EC.replace("]", "")
			else:
				EC = ""
		except ValueError:
			pass
		GeneTable.append([GeneID, EC, KOID, Name])
		GeneList.append(GeneID)
	return(GeneList, GeneTable)

# Convert downloaded PDB gene text to table
def CleanPDB(Data):
	GeneTable = []
	GeneList = []
	for Line in Data:
		try:
			GeneID, Name = Line.split(" ",1)
			Name = Name.strip()
		except ValueError:
			pass
		GeneTable.append([GeneID, Name])
		GeneList.append(GeneID)
	return(GeneList, GeneTable)
