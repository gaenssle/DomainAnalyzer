import argparse


parser = argparse.ArgumentParser(description="DOMAIN ANALYZER\nThis program downloads sequences from various databases (e.g. KEGG or UniProt) via Genome.jp"
    " The input is an id from e.g PFAM or Prosite")
parser.add_argument("name", help="name of the domain")
parser.add_argument("-m", "--multiprocess", help="turn on mutltiprocessing", action="store_true")
parser.add_argument("-ask", "--askoverwrite", help="ask before overwriting files", action="store_true")
parser.add_argument("-db", "--dblist", help="list databases to be searched, separated by ',' (default: %(default)s)", default="UniProt,KEGG;PDB,swissprot")
parser.add_argument("-a", "--action", help="add actions to be conducted (default: %(default)s)", default="a")
parser.add_argument("-st", "--searchtype", help="type of the searched id (default: %(default)s)", default="pf")
parser.add_argument("-c", "--cutoff", help="min E-Value of Pfam domains (default: %(default)s)", default=0.0001, type=float)
parser.add_argument("-f", "--folder", help="name of the parent folder (default: same as 'name')")
parser.add_argument("-cs", "--clustersize", help="entries/frament files (default: %(default)s)", default=20, type=int)
parser.add_argument("-ft", "--filetype", help="type of the produced files (default: %(default)s)", default=".csv")
parser.add_argument("-sep", "--separator", help="separator between columns in the output files (default: %(default)s)", default=";")

args = parser.parse_args()
if args.folder == None:
    args.folder = args.name

while all(ch in "aidmec" for ch in args.action) == False:
    args.action = input("\nWhich action do you want to conduct?"
        "\n- a\tconduct all actions\n- i\tdownload sequence IDs\n- d\tdownload data"
        "\n- m\tdownload motif (for KEGG)\n- e\textract data\n- c\tcount data"
        "\nPlease enter any or multiple of letters (e.g 'a' or 'dme' [without ''])\n")
if "a" in args.action:
    args.action = "idmec"


print(args)



# Gene/Protein IDs are downloaded from genome.jp with the using DBGET
# Each url consits of two variables
# 1) The database it searches
# --> Can be from KEGG (genes, mgenes)
# --> Or UniProt/SWISS-PROT (uniprot, swissprot)
# 2) The type of id it searches for
# --> KEGG ID (gn, ko)
# --> UniProt/SWISS-PROT (up, sp)
# --> Enzyme related (ec [# EC], pf [PFAM], ps [Prosite], rs [RefSeq])


# Multiprocess = "y"
# Folder = "DUF5727"
# Folder = "Test"
# Folder = "DUF1735"
# Name = "PS00812"
# Name = "DUF1735"
# DBList = ["KEGG"]
# DBList = ["UniProt"]
# DBList = ["UniProt", "KEGG"]
# DBList = ["UniProt", "KEGG", "PDB", "swissprot"]
# DBList = ["PDB"]
# Action = "dm"
# SearchType = "ps"
# Cutoff = 0.0001 # E-value cutoff for KEGG to get only probable domains (default=0.0001)
# FileType = ".csv" # File types of all exported files (default="csv")
# Sep = ";" # Separator between the columns (default=";")
# Ask = False # Ask if files should be replaced (default=True)
# ClusterSize = 20 # Define the number of entires in which the download is saved (default=250)
