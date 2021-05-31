# georgySewer

This program find mutations for Sewer samples - from pileup tables.

How to run - 

make sure the following files are in the working directory

a) region table (regions.csv in this github)

b) list of uk mutations (b117muts.csv in this github)

c) Fasta for the reference sequence (added to this github)

in the CLI ->

1) sudo -s

2) conda activate nextstrain (maybe installing "xlsxwriter", "biopython" packages is required)

3) python main.py (Path)

i.e. python main.py /data3/sewer/

Program Flow - 

Iterating over all folders in files in a given path and inserting them to list IF they are also in a pre-proccessed text file that filtering them by month

We iterate file by file in this list and row by row each mutation pileup file, and searching for mutations ( there is frequency of 5>% of nucleotide that different from the ref nucleotide)

We are parsing the information from the mutation pileup file, make a heuristic translation against the ref sequence and build all mutations table.

This table is the first output - "All_mutations.csv" file

after that the program read this output and make calculations - Count, % of total, Average of the frequencies.

and sort them by protein, than Average, than percent of total.
