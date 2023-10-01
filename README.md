# Domain Analyzer
created 2023 by gaenssle
written in Python 3.8

for questions, contact A.L.O.Gaenssle@rug.nl

The program downloads all protein data available on databases such as KEGG or UniProt via Genome.jp

-> Input is a any ID, e.g. domain name (PFAM or UniProt ID)
The program can:
- Download all sequences associated with the input ID for each database via Genome.jp, including:
  * Organism and taxonomy
  * Sequence
  * Domain architecture
  * Available IDs from other databases   
- Count taxonomy
- Count domain architecture
- Extract domain sequences from protein sequences
- Save data in the following files:
  - Gene IDs
  - Gene details (organism, architecture, sequence, etc)
  - Summary domain architecture
  - Summary of taxonomic distribution
  - Fasta files (containing entire sequence)
  - Fasta files (only containing the target domain sequence)
