#!/usr/bin/python
# Written in Python 3.7 in 2022 by A.L.O. Gaenssle
# Module for importing and exporting text files

import pandas as pd
import os
import re


##------------------------------------------------------
## HELPER FUNCTIONS
##------------------------------------------------------

# Create new Folder
def CreateFolder(NewPath):
	if not os.path.exists(NewPath):
		os.makedirs(NewPath)
		print("Created folder:", NewPath)
	else:
		print("Files will be added to:", NewPath)
	return(NewPath)

# # Check file type (Species names, Genome IDs or Table of both)
# def CheckFileType(FileName):
# 	with open(FileName, 'r') as InputFile:
# 		Text = InputFile.readlines()
# 		if "\t" in Text[0]:
# 			Type = "both"
# 			print(FileName, "contains Species names and Genome IDs")
# 		elif re.match(r"T[0-9]{5}", Text[0]):
# 			Type = "ID"
# 			print(FileName, "contains Genome IDs")
# 		elif len(Text[0]) < 6:
# 			Type = "ID"
# 			print(FileName, "contains Genome IDs")
# 		else:
# 			Type = "Name"
# 			print(FileName, "contains Species names")
# 	return(Type)

# Check if the file already exists and if it should be replaced
def CheckFileExists(FileName, Ask):
	if Ask:
		Replace = "n"
	else:
		Replace = "y"
	while Replace == "n":
		if not os.path.exists(FileName):
			break
		else:
			Replace = input("\nFile " + FileName + " already exits -> should it be replaced?"
				"\n(y=yes, n=no)\n")
			while Replace not in ("y", "n"):
				Replace = input("\nPlease enter 'y' or 'n'!\n")
		if Replace == "n":
			FileName = input("\nEnter a new filename\n")
	return(FileName)

# # Copy all data from Fragment files (250-500 genes/file) into one large file
# def CombineFiles(OutputFragments, OutputFolder, OutputFile, Header):
# 	Data = []
# 	FragmentList = os.listdir(OutputFragments)
# 	for Fragment in FragmentList:
# 		print(Fragment)
# 		Data.extend(ImportNestedList(OutputFragments + "/" + Fragment))
# 	ExportNestedList(sorted(Data), OutputFolder + OutputFile, Header, Add="_all")


##------------------------------------------------------
## IMPORT FILE FUNCTION
##------------------------------------------------------

# # Import files as list (1D) [Line1, Line2, Line3]
# def ImportList(FileName, Stamp=False):
# 	with open(FileName, 'r') as InputFile:
# 		if Stamp == True:
# 			print("Import File:", FileName)
# 		List = []
# 		for Line in InputFile:
# 			if Line != "":
# 				List.append(Line.strip())
# 		return(List)


##------------------------------------------------------
## EXPORT FILE FUNCTION
##------------------------------------------------------

# Export DataFrame
def ExportDataFrame(DataFrame, FileName, Add="", Columns="", Ask=True, Header=True):
	if Add != "":
		Name, FileType = FileName.rsplit(".", 1)
		FileName = Name + Add + "." + FileType
	FileName = CheckFileExists(FileName, Ask)
	if Columns == "":
		Columns = list(DataFrame)
	DataFrame.to_csv(FileName, sep="\t", columns = Columns, index=False, header=Header)
	print("File saved as:", FileName, "\n")

