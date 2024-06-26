#!/usr/bin/python
# Written in Python 3.8 in 2023 by A.L.O. Gaenssle

import os
import re
from multiprocessing import Pool
import pandas as pd
import argparse

# Own modules
import Import_Export as IE
import Download_GenomeJP as Genome
import Download_KEGG as KEGG


## ------------------------------------------------------------------------------------------------
## INPUT ARGUMENTS --------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="DOMAIN ANALYZER"
	"\nThis program downloads sequences from databases (e.g. KEGG or UniProt) via Genome.jp"
    "\nThe input is an id from e.g PFAM or Prosite"
    "\nAll gene IDs associated to this ID from the desired databases are downloaded"
    "\nIt then downloads relevant details for each gene ID, e.g. organism and domain architecture"
    "\nThe data can be counted regarding e.g. taxonomy, domain architecture and sequence length")
parser.add_argument("name", 
	help="name of the target that should be downloaded (e.g. domain ID or kegg orthology ID)")
parser.add_argument("-m", "--multiprocess", 
	help="turn on mutltiprocessing (only for Linux)",
	action="store_true")
parser.add_argument("-ask", "--askoverwrite", 
	help="ask before overwriting files",
	action="store_true")
parser.add_argument("-db", "--dblist", 
	help="list databases to be searched, separated by ';' (default: %(default)s)", 
	default="UniProt;KEGG;PDB;swissprot")
parser.add_argument("-a", "--action", 
	help="add actions to be conducted: "
	"a=all, i=entry IDs, d=protein data, m=KEGG motif, e=extract (default: %(default)s)", 
	default="a")
parser.add_argument("-st", "--searchtype", 
	help="type of the searched id; may be pf=domain, ko=KEGG gene ID, gn=KEGG genome ID, ec=enzyme class (default: %(default)s)", 
	default="pf")
parser.add_argument("-c", "--cutoff", 
	help="min E-Value of Pfam domains (default: %(default)s)", 
	default=0.0001, 
	type=float)
parser.add_argument("-sam", "--samplesize", 
	help="max number of downloaded entries (default: %(default)s)", 
	default=0, 
	type=int)
parser.add_argument("-f", "--folder", 
	help="name of the parent folder (default: same as 'name')")
parser.add_argument("-cs", "--clustersize", 
	help="entries/frament files (default: %(default)s)", 
	default=100, 
	type=int)
parser.add_argument("-ft", "--filetype", 
	help="type of the generated files (default: %(default)s)", 
	default=".csv")
parser.add_argument("-sep", "--separator", 
	help="separator between columns in the output files (default: %(default)s)", 
	default=";")

args = parser.parse_args()

# Set folder name to searched ID if not set
if args.folder == None:
    args.folder = args.name

# Check if the given action is valid and replace with the list if == 'a'
while all(ch in "aidmec" for ch in args.action) == False:
    args.action = input("\nWhich action do you want to conduct?"
        "\n- a\tconduct all actions\n- i\tdownload sequence IDs\n- d\tdownload data"
        "\n- m\tdownload motif (for KEGG)\n- e\textract data\n- c\tcount data"
        "\nPlease enter any or multiple of letters (e.g 'a' or 'dme' [without ''])\n")
if "a" in args.action:
    args.action = "idmec"

# Check if given list of db is valid
while True:
	try:
		args.dblist = args.dblist.replace("-","").lower().split(";")
		break
	except IndexError:
		args.dblist = input("\nCannot read given list of databases"
	        "\nPlease enter the names of the databases separated by ';'"
	        "\nWithout spaces but can be in lower case or caps\n")

# Supported databases have been implemented and can be downloaded
SupportedDBs = ["uniprot", "swissprot", "kegg"]

# Check for tab separator
if args.separator in ["\\t", "tab", "'\\t'", "{tab}"]:
	args.separator = "\t"

# Set distance (amino acids) below which the domains are merged
# (They are then shown as one domain, separated by '|')
MergeDistance = 10


## ------------------------------------------------------------------------------------------------
## HELPFER FUNCTIONS ------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------
## ================================================================================================
## Speedup download by using multiprcessing, WARNING: Has to be installed first!
def MultiProcessing(IDList, Function):
	Import = []
	print(f"Download data for {len(IDList)} Items. . .")
	if __name__ == "__main__":
		with Pool(10) as pool:
			Import = pool.map(Function, IDList)
		pool.close()
		pool.join()
	return(Import)


## ================================================================================================
## Find alternative name for PFAM domain if no hits are found in the downloaded list
def FindAlternativeName(ProteinData, GivenName):
	print("\n\nWARNING\nThe domain name has not been found! Please choose an alternative name")

	TempData = ProteinData.copy()

	# Pivot dataframe to 1 domain/row to get all domain names
	TempData = pd.wide_to_long(TempData, stubnames=["Name"], 
		i=["ID"], j="Domain", sep="-", suffix=r"[\w\d]+")
	TempData.dropna(subset = ["Name"], inplace=True)

	# Count all obtained domain names
	TempData = TempData.groupby(["Name"]).size() \
		.sort_values(ascending=False) \
		.reset_index(name="Entries")

	# Print and request index of alternative domain name
	print(TempData)
	Answer = input("\nFound domain names in the obtained data." +
		f"\nPlease enter the index of the correct alternative name (left number, 0 to {str(len(TempData)-1)})" +
		"\nTo quit, write 'quit'\n")
	while True:
		if Answer == "quit":
			quit()
		elif Answer.isdigit() and int(Answer) in range(len(TempData)):
			AlternativeName = TempData["Name"][int(Answer)]
			print(f"Selected alternative name is: {AlternativeName}\n\n")
			break
		else:
			Answer = input(f"\nPlease enter either 'quit' or an index between 0 and {str(len(TempData)-1)}\n")
	return(AlternativeName)


## ------------------------------------------------------------------------------------------------
## MAIN FUNCTIONS ---------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------
## ================================================================================================
## Download all gene IDs associated with the supplied domain name from KEGG, UniProt and PDB
def DownloadList(Name, OutputFile, DB, SearchType, FileType, Sep, Ask):
	print(f"Download gene IDs containing domain {Name} . . .", 
		"\n(Speed depends on internet connection)\n")

	# Create the required url with the following fragments
	urlGeneList = "https://www.genome.jp/dbget-bin/get_linkdb?-t+"
	urlName = "+" + SearchType + ":"
	# SearchType = "+ps:"
	urlPage = "+-p+"
	urlKEGG = "genes"
	AddToName = "_" + DB
	DB = DB.replace("kegg", urlKEGG).lower()
	urlInitial = urlGeneList + DB + urlName + Name

	# Download all genes from the first page and then cycle through all subsequent pages
	IDList, Pages = Genome.DownloadGeneList(urlInitial, getAmount=True)
	for Page in range(2, Pages+1):
		urlNew = urlGeneList + DB + urlPage + str(Page) + urlName + Name
		IDList.extend(Genome.DownloadGeneList(urlNew))

	# Extract all information from the entries and convert into pandas dataframe
	if IDList:
		if DB == "uniprot" or DB == "swissprot":
			GeneTable = Genome.CleanUniProt(IDList)
		elif DB == "genes":
			GeneTable = Genome.CleanKEGG(IDList)
		else:
			GeneTable = Genome.CleanPDB(IDList)
		if len(GeneTable.index) != 0:
			print(f"Data found for {len(GeneTable.index)} entries\n")
			IE.ExportDataFrame(GeneTable, OutputFile, 
				Add=AddToName, FileType=FileType, Sep=Sep, Ask=Ask)
			if DB == "genes":
				IE.ExportDataFrame(GeneTable, OutputFile, 
					Add=AddToName + "_GeneIDs", Columns= ["ID"], FileType=FileType, Sep=Sep, Ask=Ask, Header=False)
	return(bool(IDList))


## ================================================================================================
## Download details for all given IDs from UniProt, including taxonomy, sequence and names
def DownloadEntryUniProt(IDList, FilePath, FileType, Sep, Multiprocess, ClusterSize, Ask):
	print(f"Download protein data for {len(IDList)} . . .")

	# Create clusters of sequences to generate smaller files (in case the download crashes)
	ClusteredList = [IDList[x:x+ClusterSize] for x in range(0, len(IDList), ClusterSize)]
	for ClusterID in range(len(ClusteredList)):
		print(f"Download cluster {ClusterID+1} of {len(ClusteredList)}")
		FragmentFile = FilePath + "_" + str(ClusterID+1)
		print(FragmentFile)

		# Ignore all files that have already been downloaded
		if os.path.exists(FragmentFile + FileType):
			print("File already exists, skip to next cluster\n")

		# Download all files that have not yet been saved
		else:
			if Multiprocess:
				ListOfDicts = MultiProcessing(ClusteredList[ClusterID], Genome.DownloadEntryUniProt)
			else:
				ListOfDicts = []
				for ID in ClusteredList[ClusterID]:
					Dict =  Genome.DownloadEntryUniProt(ID)
					ListOfDicts.append(Dict)

			# Convert to pandas dataframe and store each fragment in the fragment folder
			ProteinTable = pd.DataFrame(ListOfDicts)
			print(f"Done!\n-> {len(ProteinTable)} of {len(ClusteredList[ClusterID])} found")
			IE.ExportDataFrame(ProteinTable, FragmentFile, FileType=FileType, Sep=Sep, Ask=Ask)

	# After all entries have been downloaded, combine all fragments into one dataframe
	DataFrame = IE.CombineFiles(os.path.split(FragmentFile)[0], Sep, FileType)
	return(DataFrame)


## ================================================================================================
## Download details for all given IDs from KEGG, including taxonomy and sequence
def DownloadEntryKEGG(IDList, FilePath, FileType, Sep, Multiprocess, ClusterSize, Ask):
	Organisms = None
	print(f"Download protein data for {len(IDList)} . . .")

	# Create clusters of sequences to generate smaller files (in case the download crashes)
	ClusteredList = [IDList[x:x+ClusterSize] for x in range(0, len(IDList), ClusterSize)]
	for ClusterID in range(len(ClusteredList)):
		print(f"Download cluster {ClusterID+1} of {len(ClusteredList)}")
		FragmentFile = FilePath + "_" + str(ClusterID+1)
		print(FragmentFile)

		# Create chunks of clusters since data of 10 proteins can be downloaded from KEGG at once
		ClusterChunk = [ClusteredList[ClusterID][x:x+10] for x in range(0, len(ClusteredList[ClusterID]), 10)]

		# Ignore all files that have already been downloaded
		if os.path.exists(FragmentFile + FileType):
			print("File already exists, skip to next cluster\n")
		else:
			# Multiprocessing saves the entries as a list of list of dictionaries
			if Multiprocess:
				ClusteredListOfDicts = MultiProcessing(ClusterChunk, KEGG.DownloadProteinEntries)
				ListOfDicts = [Entry for Cluster in ClusteredListOfDicts for Entry in Cluster]
			else:
				ListOfDicts = []
				for Chunk in ClusterChunk:
					Set =  KEGG.DownloadProteinEntries(Chunk)
					ListOfDicts.extend(Set)

			# Only download the list of organisms on KEGG if needed and add to dataframe
			if Organisms is None:
				Organisms = KEGG.DownloadOrganismsTemp()
			ProteinTable = pd.DataFrame(ListOfDicts)
			ProteinTable = pd.merge(ProteinTable, Organisms, on=["orgID"],  how="left")
			IE.ExportDataFrame(ProteinTable, FragmentFile, FileType=FileType, Sep=Sep, Ask=Ask)
			print(f"Done!\n-> {len(ProteinTable)} of {len(ClusteredList[ClusterID])} found")

	# After all entries have been downloaded, combine all fragments into one dataframe
	DataFrame = IE.CombineFiles(os.path.split(FragmentFile)[0], Sep, FileType)
	return(DataFrame)


## ================================================================================================
## Download all domain motifs for all given IDs from KEGG
def DownloadMotifKEGG(IDList, FilePath, CutOff, FileType, Sep, Multiprocess, ClusterSize, Ask):
	print(f"Download motif for {len(IDList)} . . .")

	# Create clusters of sequences to generate smaller files (in case the download crashes)
	ClusteredList = [IDList[x:x+ClusterSize] for x in range(0, len(IDList), ClusterSize)]
	for ClusterID in range(len(ClusteredList)):
		print(f"Download cluster {ClusterID+1} of {len(ClusteredList)}")
		FragmentFile = FilePath + "_" + str(ClusterID+1)
		print(FragmentFile)

		# Ignore all files that have already been downloaded
		if os.path.exists(FragmentFile + FileType):
			print("File already exists, skip to next cluster")
		else:
			ListOfLists = []
			if Multiprocess:
				ClusteredListOfLists = MultiProcessing(ClusteredList[ClusterID],  KEGG.DownloadMotif)
				ListOfLists = [Entry for Cluster in ClusteredListOfLists for Entry in Cluster]
			else:
				for ID in ClusteredList[ClusterID]:
					ListOfLists.extend(KEGG.DownloadMotif(ID))
			ColNames= ["ID", "Index","Name", "Start", "End", "Definition", "E-Value", "Score"]
			MotifTable = pd.DataFrame(ListOfLists, columns=ColNames)
			IE.ExportDataFrame(MotifTable, FragmentFile, FileType=FileType, Sep=Sep, Ask=Ask)
			print(f"Done!\n-> {len(MotifTable)} domains for {len(ClusteredList[ClusterID])} entries found")

	# After all entries have been downloaded, combine all fragments into one dataframe
	DataFrame = IE.CombineFiles(os.path.split(FragmentFile)[0], Sep, FileType)

	# Remove all rows with E-Values that are empty ("-") or below the cutoff and condense the drataframe to 1 row/protein
	DataFrame = DataFrame.sort_values(by=["ID", "Start"])	# Sort domains by start position to ensure correct order
	DataFrame["E-Value"] = pd.to_numeric(DataFrame["E-Value"], errors="coerce") # Remove empty-E-Values and convert to floats
	ConcatDataFrame = DataFrame.copy()	# Make a copy to avoid errors
	ConcatDataFrame = ConcatDataFrame[ConcatDataFrame["E-Value"] < CutOff] 	# Remove all rows with E-Values below CutOff
	ConcatDataFrame["Index"] = ConcatDataFrame.groupby(["ID"]).cumcount()+1	# Reset the domain counting for the remaining domains

	# Find and merge overlapping domains (starting/ending within the MergeDistance from each other)
	ConcatDataFrame = KEGG.MergeOverlapping(ConcatDataFrame, MergeDistance)

	# Move info for all domains for the same protein to one row (pivot dataframe and rename domain columns)
	ConcatDataFrame = ConcatDataFrame.pivot(index="ID", columns="Index", 
		values=["Name", "Start", "End", "E-Value"])
	ConcatDataFrame = ConcatDataFrame.sort_index(axis=1, level=1, sort_remaining=False) # Ensure correct order of columns
	ConcatDataFrame.columns = [f"{x}-D{y}" for x, y in ConcatDataFrame.columns]	# Give the new columns names starting with D[n]-
	return(DataFrame, ConcatDataFrame)


## ================================================================================================
## Remove all entries missing the target domain, including overlapping domains
def FilterDomains(ProteinData, DomainNameCols, DomainName):
	AlternativeName = ""
	Hit = False

	# Get input or alternative name and remove all entries that miss the domain
	for Column in DomainNameCols:
		if ProteinData[Column].str.contains(DomainName).any():
			Hit = True
	if Hit == True:
		SearchPhrase = DomainName
	else:
		AlternativeName = FindAlternativeName(ProteinData, DomainName)
		SearchPhrase = AlternativeName
	ProteinData["Hit"] = ProteinData[DomainNameCols].apply(lambda x: 
		x.str.contains(SearchPhrase,case=False)).any(axis=1).astype(int)

	ProteinData.drop(ProteinData.loc[ProteinData["Hit"]==0].index, inplace=True)
	ProteinData.dropna(axis=1, how="all", inplace=True)
	ProteinData.drop(["Hit"], axis=1, inplace=True)

	# Create a dataframe that only contains the target domains (pivot to 1 domain/row)
	Domains = ProteinData.copy()
	DomainCols = ["Name", "Start", "End", "E-Value"]
	Domains = pd.wide_to_long(Domains, stubnames=DomainCols, 
		i=["ID"], j="Domain", sep="-", suffix=r"[\w\d]+")
	Domains.dropna(subset = ["Name"], inplace=True)
	if AlternativeName == "":
		Domains = Domains[Domains["Name"].str.contains(DomainName)]
	else:
		Domains = Domains[Domains["Name"].str.contains(AlternativeName)]

	# Extract the amino acid sequence for each domain
	Domains[["Start", "End"]] = Domains[["Start", "End"]].astype(int)
	Domains["Sequence"] = Domains.apply(lambda x: 
		x["Sequence"][x["Start"]:x["End"]+1], axis=1)
	Domains["Length"] = Domains["Sequence"].apply(len)

	return(ProteinData, Domains)


## ------------------------------------------------------------------------------------------------
## SCRIPT -----------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------
# Print program header
print('\n{:=<70}'.format(''))
print('{:=^70}'.format('  DOMAIN ANALYZER  '))
print('{:=<70}'.format(''))
print('{: ^70}\n\n'.format('2024, by A.L.O. Gaenssle'))


IE.CreateFolder(os.path.join(args.folder, "Input"))
IE.CreateFolder(os.path.join(args.folder, "Output"))



for DB in args.dblist:
	FileName = args.name + "_" + DB

	# Download Sequence IDs from e.g. UniProt, KEGG and/or PDB
	if any(s in ["i", "d", "m"] for s in args.action):
		InputFile = os.path.join(args.folder, "Input", FileName + args.filetype)
		if not os.path.exists(InputFile):
			print(f"\nNo input file found for {DB.upper()}")
			OutputFile = os.path.join(args.folder, "Input", args.name)
			Hits = DownloadList(args.name, OutputFile, DB, 
				args.searchtype, args.filetype, args.separator, args.askoverwrite)
		if os.path.exists(InputFile):
			print(f"\nInput file for {DB.upper()} already exists, loading file . . .\n({InputFile})")
			DataFrame = pd.read_csv(InputFile, sep=args.separator)
			IDList = DataFrame["ID"].tolist()
			Hits = bool(IDList)
			if args.samplesize != 0 and args.samplesize <= len(IDList):
				IDList = IDList[:args.samplesize]
		if Hits == False or DB not in SupportedDBs:
			continue


	# download protein data from UniProt and/or KEGG
	if "d" in args.action:
		OutputPath = os.path.join(args.folder, "Output", FileName + "_Protein")
		FragmentFolder = IE.CreateFolder(OutputPath + "Fragments")
		FragmentFile = os.path.join(FragmentFolder, FileName + "_Protein")
		if DB == "kegg":
			Detailed = DownloadEntryKEGG(IDList, FragmentFile, args.filetype, 
				args.separator, args.multiprocess, args.clustersize, args.askoverwrite)
		else:
			Detailed = DownloadEntryUniProt(IDList, FragmentFile, args.filetype, 
				args.separator, args.multiprocess, args.clustersize, args.askoverwrite)
		DataFrame = pd.merge(DataFrame, Detailed, on=["ID"],  how="right")
		IE.ExportDataFrame(DataFrame, OutputPath, 
			FileType=args.filetype, Sep=args.separator, Ask=args.askoverwrite)


	# Download motif data from KEGG
	if "m" in args.action and DB == "kegg":
		OutputPath = os.path.join(args.folder, "Output", FileName + "_Motif")
		FragmentFolder = IE.CreateFolder(OutputPath + "Fragments")
		FragmentFile = os.path.join(FragmentFolder, FileName + "_Motif")
		Motifs, ConcatMotifs = DownloadMotifKEGG(IDList, FragmentFile, args.cutoff, 
			args.filetype, args.separator, args.multiprocess, args.clustersize, args.askoverwrite)
		IE.ExportDataFrame(Motifs, OutputPath, 
			FileType=args.filetype, Sep=args.separator, Ask=args.askoverwrite)

		# Remove all entries that don't contain the target PFAM domain
		if args.searchtype == "pf":	
			DomainNameCols = [col for col in ConcatMotifs if col.startswith("Name-")]
			ConcatMotifs["Hit"] = ConcatMotifs[DomainNameCols].apply(lambda x: 
				x.str.contains(args.name,case=False)).any(axis=1).astype(int)
			ConcatMotifs.drop(ConcatMotifs.index[ConcatMotifs["Hit"] == 0], inplace = True)
			ConcatMotifs.dropna(how="all", axis=1, inplace=True)
			ConcatMotifs.drop(columns=["Hit"], inplace=True)
		try:
			ProteinFile = os.path.join(args.folder, "Output", FileName + "_Protein" + args.filetype)
			ProteinData = pd.read_csv(ProteinFile, sep=args.separator)
			DataFrame = pd.merge(ProteinData, ConcatMotifs, on=["ID"],  how="right")
		except FileNotFoundError:
			DataFrame = pd.merge(DataFrame, ConcatMotifs, on=["ID"],  how="right")
			print("Motif data is stored without protein data since the protein file is missing")
			continue
		CutoffFilePath = os.path.join(args.folder, "Output", 
			FileName + "_Protein_cutoff" + str(args.cutoff))
		IE.ExportDataFrame(DataFrame, CutoffFilePath, 
			FileType=args.filetype, Sep=args.separator, Ask=args.askoverwrite)
		CutoffIDFilePath = os.path.join(args.folder, "Input", 
			FileName + "_GeneIDs_cutoff" + str(args.cutoff))
		IE.ExportDataFrame(DataFrame, CutoffIDFilePath, Columns=["ID"],
			FileType=args.filetype, Sep=args.separator, Ask=args.askoverwrite, Header=False)


	# Extract data from UniProt and/or KEGG
	if "e" in args.action:
		FilePath = os.path.join(args.folder, "Output", FileName)
		if DB == "kegg":
			InputFile = FilePath + "_Protein_cutoff" + str(args.cutoff) + args.filetype
		else:
			InputFile = FilePath + "_Protein" + args.filetype
		try:
			ProteinData = pd.read_csv(InputFile, sep=args.separator)
			ProteinData.dropna(axis=1, how="all", inplace=True)
			DomainNameCols = [col for col in ProteinData if col.startswith("Name-")]
		except FileNotFoundError:
			print(f"No input file found for database: {DB}\n\n")
			continue

		# Remove all entries missing the target domain, including overlapping domains
		if args.searchtype == "pf":
			ProteinData, Domains = FilterDomains(ProteinData, DomainNameCols, args.name)
			IE.CreateFasta(Domains, FilePath, only=True)
			IE.ExportDataFrame(Domains, FilePath + "_Domain_Details", Index=True,
				FileType=args.filetype, Sep=args.separator, Ask=args.askoverwrite)

		# Create Fasta files of complete sequences
		IE.CreateFasta(ProteinData, FilePath)

		# Create a file with the domain architecture (domains are joined by '+')
		DomainNameCols = [col for col in ProteinData if col.startswith("Name-")]
		MotifCols = ["ID", "Organism", "Taxonomy", "Length"] + DomainNameCols
		Motifs = ProteinData[MotifCols].copy()
		Motifs["Domains"] = Motifs[DomainNameCols].apply(lambda x: "+".join(x.dropna()), axis=1)
		Motifs.drop(columns=DomainNameCols, inplace=True)
		IE.ExportDataFrame(Motifs, FilePath + "_DomainArchitecture", 
			FileType=args.filetype, Sep=args.separator, Ask=args.askoverwrite)
		
		# Count the number of entries for each taxonomy, organism and domain architecture
		ColName = "Entries"
		CountDomains = Motifs.groupby(["Domains"]).size() \
			.reset_index(name=ColName) \
			.sort_values([ColName], ascending=False)
		CountTax = Motifs.groupby(["Taxonomy"]).size() \
			.reset_index(name=ColName) \
			.sort_values([ColName], ascending=False)
		CountOrganisms = Motifs.groupby(["Taxonomy","Organism"]).size() \
			.reset_index(name=ColName) \
			.sort_values(["Taxonomy",ColName], ascending=False)
		IE.ExportDataFrame(CountDomains, FilePath + "_CountDomains", 
			FileType=args.filetype, Sep=args.separator, Ask=args.askoverwrite)
		IE.ExportDataFrame(CountTax, FilePath + "_CountTax", 
			FileType=args.filetype, Sep=args.separator, Ask=args.askoverwrite)
		IE.ExportDataFrame(CountOrganisms, FilePath + "_CountOrganisms", 
			FileType=args.filetype, Sep=args.separator, Ask=args.askoverwrite)

print('{:=^70}'.format('  End of program  '))