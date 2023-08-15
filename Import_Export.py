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

# Copy all data from Fragment files (250-500 genes/file) into one large file
def CombineFiles(Folder, Sep):
	FileList = os.listdir(Folder)
	DataList = []
	for File in FileList:
		FilePath = os.path.join(Folder, File)
		DataFrame = pd.read_csv(FilePath, sep=Sep)
		DataList.append(DataFrame)
	DataFrame = pd.concat(DataList, axis=0, ignore_index=True)
	DataFrame = DataFrame.sort_values(DataFrame.columns[0])
	return(DataFrame)


##------------------------------------------------------
## IMPORT FILE FUNCTION
##------------------------------------------------------

# # Import pandas dataframe
# def ImportDataFrame(FileName, UseCols=[], FileType=".csv", Sep=";", Stamp=False):
# 	FileName = FileName + FileType
# 	if Stamp == True:
# 		print("Import File:", FileName)
# 	if UseCols == []:
# 		DataFrame = pd.read_csv(FileName, sep=Sep)
# 	else:
# 		DataFrame = pd.read_csv(FileName, sep=Sep, usecols=UseCols)
# 	return(DataFrame)


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

# Export pandas dataframe
def ExportDataFrame(DataFrame, FileName, Add="", Columns="", FileType=".csv", Sep=";", Ask=True, Header=True):
	FileName = FileName + Add + FileType
	FileName = CheckFileExists(FileName, Ask)
	if Columns == "":
		Columns = list(DataFrame)
	DataFrame.to_csv(FileName, sep=Sep, columns = Columns, index=False, header=Header)
	print("File saved as:", FileName, "\n")
