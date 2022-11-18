#!/usr/bin/python
# Written in Python 3.7 in 2022 by A.L.O. Gaenssle
# Ask user for what to do

import os

# Own modules
import Import_Export as IE

# Print header with info on program
def PrintHeader():
	print("\n","-"*75,"\n","-"*75,"\n")
	print("DOMAIN ANALYZER\tby A.L.O. Gaenssle, 2022")
	print("\n","-"*75,"\n","-"*75,"\n")
	print("This program help automating sequence analysis "
		"of domains")
	print("\nThe program has been written to:"
		"\n- Download sequences from KEGG and UniProt "
		"(based on domain name)"
		"\n- Count taxonomy and domain architecture"
		"\n- Extract domain sequences from protein sequences"
		"\n- Save data in the following files:\n\t- Gene IDs"
		"\n\t- Gene details\n\t- Fasta files"
		"\n\t- Genes split into domains"
		"\n\t- Files summarizing domain architecture and taxonomy")
	print("\nContact A. Lucie Gaenssle for help and adaptions "
		"(A.L.O.Gaenssle@rug.nl)")
	print("\n","-"*75)
	print("\nInformation:\n- Complete input by clicking enter"
		"\n- Navigate within and between the inputs using the arrow keys"
		"\n- Terminate the program any time by:"
		"\n\t- Closing the terminal (window)"
		"\n\t- Ctrl + C")
	print("\n","-"*75,"\n QUESTIONAIRE\n","-"*75)


# Get name of the folder the data should be saved in
def GetFolder(DirectoryName):
	Confirmation = False
	Folder = input("\nEnter your folder (project) name:"
		"\n- If folder is a subfolder new to the python script: e.g. Test"
		"\n- Otherwise enter full path: e.g. "
		"X:\Test\DUF1735"
		"\n\nYour current directory is:\n%s\n"
		% DirectoryName)
	while Confirmation == False:
		if os.path.exists(Folder):
			Process = input("\nThis folder already exists"
				"\nDo you want to:\n- a\tappend/replace data to this folder\n- n\tenter new folder name\n")
			while Process not in ("a", "n"):
				Process = input("\nPlease enter 'a' (append) or 'n' (new)\n")
			if Process == "a":
				break
			else:
				Folder = input("\nEnter your folder (project) name:\n")
		else:
			Process = input("\nThis folder does not exist yet"
					"\nDo you want to:\n- c\tcreate a new project\n- n\tenter new folder name\n")
			while Process not in ("c", "n"):
				Process = input("\nPlease enter 'c' (create) or 'n' (re-enter)\n")
			if Process == "c":
				break
		if Process == "n":
			Folder = input("\nEnter your folder (project) name:\n")
	IE.CreateFolder(Folder + "/Input")
	IE.CreateFolder(Folder + "/Output")
	# print("Created folder:", Folder, "with the Subfolders Input and Output")
	return(Folder)

# Ask from which databases should be search (options: UniProt, KEGG and PDB)
def GetDB():
	DBList = input("\nWhich databases do you want to download from?"
		"\n- a\tall (KEGG, UniProt, PDB)\n- k\tKEGG\n- u\tUniProt\n- p\tPDB (only for sequence IDs)"
		"\n-> enter e.g 'a', 'k' or 'ku'\n")
	while all(ch in "akup" for ch in DBList) == False:
		DBList = input("\nPlease enter 'a', 'k', 'u' or 'p' (alone or combined)\n")
	if "a" in DBList:
		DBList = "kup"
	DBList = DBList.replace("k", "KEGG\t").replace("p", "PDB\t").replace("u", "UniProt\t")
	DBList = DBList.split("\t")[:-1]
	return(DBList)

# Ask which type of action should be conducted (download genes, extract/count data, etc)
def GetAction():
	DBList = GetDB()
	Action = input("\nWhich data would you like to add/replace?"
		"\n- a\tconduct all actions\n- i\tdownload sequence IDs\n- d\tdownload data"
		"\n- m\tdownload motif (for KEGG)\n- e\textract data\n- c\tcount data"
		"\n-> enter e.g 'a', 'e', or 'dme'\n")
	while all(ch in "aidmec" for ch in Action) == False:
		Action = input("\nPlease enter any or multiple of the following: a, i, d, m, e, c\n")
	if "a" in Action:
		Action = "idmec"
	return(DBList, Action)

# Get the name of the domain which should be search for, e.g. DUF1735
def GetName(Folder):
	Verified = False
	try:
		Name = os.listdir(Folder+"/Input")[0].split("_",1)[0]
	except:
		Name = input("\nEnter the name/ID of the domain (e.g. DUF1735)\n")
	while Verified == False:
		Correct = input("\nThe domain name is %s -> is this correct?\n(y=yes, n=no)\n" % Name)
		if Correct == "y":
			break
		else:
			Name = input("\nEnter the name/ID of the domain (e.g. DUF1735)\n")
	return(Name)

# Ask if multiprocessing should be used
def UseMultiprocess():
	Confirmation = input("\nUse multiprocess? -> Speeds up process 8x but can crash"
		"\n-> also needs to be installed first (module: multiprocessing)\n(y=yes, n=no)\n")
	while Confirmation not in ("y", "n"):
		Confirmation = input("\nPlease enter 'y' or 'n'!\n")
	return(Confirmation)

def PrintFooter():
	print("\n","-"*75,"\n END OF PROGAM\n","-"*75)
