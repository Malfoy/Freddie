#!/usr/bin/env python3

#kmer sizes to test below 63
low_kmer_sizes = ["21","31","41","51","61"]

#kmer sizes to test above 63
high_kmer_sizes = ["71","81","91","101"]

#Query file to evaluate
contig_file="contigs.fa" 

#File of file containing parents fasta files
parents_file_of_file="fof.txt"

#Read N times all the files but uses N time less memory
number_iteration="1"

#Add the results to this text file
out_file="results.txt"

# Freddie_path=""
Freddie_path="./"

for w in low_kmer_sizes:
    command=Freddie_path+"Freddie_MS "+parents_file_of_file+" "+contig_file+" "+w+" "+number_iteration+" >> "+out_file+" ;"
    print(command)

for w in high_kmer_sizes:
    command=Freddie_path+"Jason_MS "+parents_file_of_file+" "+contig_file+" "+w+" "+number_iteration+" >> "+out_file+" ;"
    print(command)
