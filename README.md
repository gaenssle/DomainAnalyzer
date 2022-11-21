# Domain Analyzer
created 2022 by gaenssle
written in Python 3.7

for questions, write to A.L.O.Gaenssle@rug.nl

This program downloads all protein data available on KEGG and UniProt via Genome.jp

-> Input is a domain name (PFAM or UniProt ID)
The program can:
- Count taxonomy
- Count domain architecture
- Extract domain sequences from protein sequences
- Save data in the following files:
  - Gene IDs
  - Gene details (organism, architecture, sequence, etc)
  - Summary domain architecture
  - Summary of taxonomic distribution
  - Fasta files (only containing the target domain sequence)
