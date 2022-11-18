#!/usr/bin/python
# Written in Python 3.7 in 2022 by A.L.O. Gaenssle
# Module for importing and exporting text files

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

# Check file type (Species names, Genome IDs or Table of both)
def CheckFileType(FileName):
	with open(FileName, 'r') as InputFile:
		Text = InputFile.readlines()
		if "\t" in Text[0]:
			Type = "both"
			print(FileName, "contains Species names and Genome IDs")
		elif re.match(r"T[0-9]{5}", Text[0]):
			Type = "ID"
			print(FileName, "contains Genome IDs")
		elif len(Text[0]) < 6:
			Type = "ID"
			print(FileName, "contains Genome IDs")
		else:
			Type = "Name"
			print(FileName, "contains Species names")
	return(Type)

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
def CombineFiles(OutputFragments, OutputFolder, OutputFile, Header):
	Data = []
	FragmentList = os.listdir(OutputFragments)
	for Fragment in FragmentList:
		print(Fragment)
		Data.extend(ImportNestedList(OutputFragments + "/" + Fragment))
	ExportNestedList(sorted(Data), OutputFolder + OutputFile, Header, Add="_all")


##------------------------------------------------------
## IMPORT FILE FUNCTIONS
##------------------------------------------------------

# Import files as list (1D) [Line1, Line2, Line3]
def ImportList(FileName, Stamp=False):
	with open(FileName, 'r') as InputFile:
		if Stamp == True:
			print("Import File:", FileName)
		List = []
		for Line in InputFile:
			if Line != "":
				List.append(Line.strip())
		return(List)

# Import file as nested list (2D) [[L1-C1, L1-C2], [L2-C1, L2-C2]]
def ImportNestedList(FileName, Stamp=False, getHeader=False):
	with open(FileName, 'r') as InputFile:
		if Stamp == True:
			print("Import File:", FileName)
		Header = next(InputFile)
		List = []
		for Line in InputFile:
			Line = Line.rstrip()
			if Line != "":
				List.append(Line.split("\t"))
		if getHeader:
			return(List, Header)
		else:
			return(List)

# Import file as dictionary {L1-C1:L1-C2, L2-C1:L2-C2}
def ImportDictionary(FileName, Stamp=False, getHeader=False):
	with open(FileName, 'r') as InputFile:
		if Stamp == True:
			print("Import File:", FileName)
		Header = next(InputFile)
		Dictionary = {}
		for Line in InputFile:
			Line = Line.strip()
			if Line != "" and "\t" in Line:
				Key, Value = Line.split("\t")
				Dictionary[Key] = Value
		if getHeader:
			return(Dictionary, Header)
		else:
			return(Dictionary)

# Import file as dictionary {L1C1:[L1C2, L1C3],L2C1:[L2C2, L2C3]}
def ImportNestedDictionary(FileName, Stamp=False, getHeader=False):
	with open(FileName, 'r') as InputFile:
		if Stamp == True:
			print("Import File:", FileName)
		Header = next(InputFile)
		Dictionary = {}
		for Line in InputFile:
			Line = Line.strip()
			if Line != "" and "\t" in Line:
				List = Line.split("\t")
				Dictionary[List[0]] = List[1:]
		if getHeader:
			return(Dictionary, Header)
		else:
			return(Dictionary)


# Import file as dictionary with same names in C1 {N1:[[L1],[L2]], N2:[[Ln],[Lm]]}
def ImportDoubleNestedDictionary(FileName, Stamp=False, getHeader=False):
	with open(FileName, 'r') as InputFile:
		if Stamp == True:
			print("Import File:", FileName)
		Header = next(InputFile)
		Dictionary = {}
		for Line in InputFile:
			Line = Line.strip()
			if Line != "" and "\t" in Line:
				List = Line.split("\t")
				if List[0] not in Dictionary:
					Dictionary[List[0]] = [List[1:]]
				else:
					Dictionary[List[0]].append(List[1:])
		if getHeader:
			return(Dictionary, Header)
		else:
			return(Dictionary)


##------------------------------------------------------
## EXPORT FILE FUNCTIONS
##------------------------------------------------------

# Export list (1D) as text file [Line1, Line2, Line3]
def ExportList(Data, FileName, Add="", Ask=True):
	if Add != "":
		Part1, Part2 = FileName.rsplit(".", 1)
		FileName = Part1 + Add + "." + Part2
	FileName = CheckFileExists(FileName, Ask)
	OutputFile = open(FileName, "w")
	for Item in Data:
		String = Item + "\n"
		OutputFile.write(String)
	print("File saved as:", FileName, "\n")

# Export nested list (2D) as text file [[L1-W1, L1-W2], [L2-W1, L2-W2]]
def ExportNestedList(Data, FileName, Header, Add="", Ask=True):
	if Add != "":
		Part1, Part2 = FileName.rsplit(".", 1)
		FileName = Part1 + Add + "." + Part2
	FileName = CheckFileExists(FileName, Ask)
	OutputFile = open(FileName, "w")
	OutputFile.write(Header)
	for Item in Data:
		Item = [str(i) for i in Item]
		String = "\t".join(Item) + "\n"
		OutputFile.write(String)
	print("File saved as:", FileName, "\n")

# Export double nested list (3D) as text file [[[111,112],[121,122]],[[211,212],[221,222]]]
def ExportDoubleNestedList(Data, FileName, Header, Add="", Ask=True):
	if Add != "":
		Part1, Part2 = FileName.rsplit(".", 1)
		FileName = Part1 + Add + "." + Part2
	FileName = CheckFileExists(FileName, Ask)
	OutputFile = open(FileName, "w")
	OutputFile.write(Header)
	for Item in Data:
		for SubItem in Item:
			SubItem = [str(i) for i in SubItem]
			String = str(Data.index(Item)+1) + "\t" + "\t".join(SubItem) + "\n"
			OutputFile.write(String)
	print("File saved as:", FileName, "\n")

# Export dictionary as text file {A:1, B:2, C:3}
def ExportDictionary(Data, FileName, Header, Add="", Ask=True):
	if Add != "":
		Part1, Part2 = FileName.rsplit(".", 1)
		FileName = Part1 + Add + "." + Part2
	FileName = CheckFileExists(FileName, Ask)
	OutputFile = open(FileName, "w")
	OutputFile.write(Header)
	for Name, ID in Data.items():
		String = str(Name) + "\t" + str(ID) + "\n"
		OutputFile.write(String)
	print("File saved as:", FileName, "\n")

# Export nested dictionary as text file {A:[1,2,3], B:[1,2,3]}
def ExportNestedDictionary(Data, FileName, Header, Add="", Ask=True):
	if Add != "":
		Part1, Part2 = FileName.rsplit(".", 1)
		FileName = Part1 + Add + "." + Part2
	FileName = CheckFileExists(FileName, Ask)
	OutputFile = open(FileName, "w")
	OutputFile.write(Header)
	for Key in Data:
		Data[Key] = [str(i) for i in Data[Key]]
		String = str(Key) + "\t" + "\t".join(Data[Key]) + "\n"
		OutputFile.write(String)
	print("File saved as:", FileName, "\n")

# Export double nested dictionary as text file {A:[[11,12],[21,22]],B:[[11,12],[21,22]]}
def ExportDoubleNestedDictionary(Data, FileName, Header, Add="", Ask=True):
	if Add != "":
		Part1, Part2 = FileName.rsplit(".", 1)
		FileName = Part1 + Add + "." + Part2
	FileName = CheckFileExists(FileName, Ask)
	OutputFile = open(FileName, "w")
	OutputFile.write(Header)
	for Key in Data:
		for SubItem in Data[Key]:
			SubItem = [str(i) for i in SubItem]
			String = str(Key) + "\t" + "\t".join(SubItem) + "\n"
			OutputFile.write(String)
	print("File saved as:", FileName, "\n")
