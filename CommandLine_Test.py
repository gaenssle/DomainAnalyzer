import argparse


parser = argparse.ArgumentParser(description="calculate X to the power of Y")
group = parser.add_mutually_exclusive_group()
group.add_argument("-v", "--verbose", action="store_true")
group.add_argument("-q", "--quiet", action="store_true")
parser.add_argument("x", type=int, help="the base")
parser.add_argument("y", type=int, help="the exponent")
args = parser.parse_args()
answer = args.x**args.y

if args.quiet:
    print(answer)
elif args.verbose:
    print("{} to the power {} equals {}".format(args.x, args.y, answer))
else:
    print("{}^{} == {}".format(args.x, args.y, answer))

parser = argparse.ArgumentParser(description="download")
parser.add_argument("name", help="name of the domain")
parser.add_argument("-f", "--folder", help="name of the parent folder")

# Gene/Protein IDs are downloaded from genome.jp with the using DBGET
# Each url consits of two variables
# 1) The database it searches
# --> Can be from KEGG (genes, mgenes)
# --> Or UniProt/SWISS-PROT (uniprot, swissprot)
# 2) The type of id it searches for
# --> KEGG ID (gn, ko)
# --> UniProt/SWISS-PROT (up, sp)
# --> Enzyme related (ec [# EC], pf [PFAM], ps [Prosite], rs [RefSeq])
