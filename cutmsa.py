#!/usr/bin/env python3

from position import *

msa_file="subsampled_alignment.fasta"
query_sequence="ATGTACTCATTCGTTTCGGAAGAGACAGGTACGTTAATAGTTAATAGCGTACTTCTTTTTCTTGCTTTCGTGGTATTCTTGCTAGTTACACTAGCCATCCTTACTGCGCTTCGATTGTGTGCGTACTGCTGCAATATTGTTAACGTGAGTCTTGTAAAACCTTCTTTTTACGTTTACTCTCGTGTTAAAAATCTGAATTCTTCTAGAGTTCCTGATCTTCTGGTCTAA"

def cutMsa(file_msa,positions):
	salida = open(file_msa.split(".")[:-1][0]+".cut.fasta","w")
	for record in SeqIO.parse(file_msa, "fasta"):
		secuencia=record.seq			
		i=0
		cut_secuencia=''
		check=False
		for i in range(len(secuencia)):
			if i == positions[0]:
				check=True
			if check==True and i<positions[1]:
				cut_secuencia=cut_secuencia+secuencia[i]
		salida.write(record.id+"\n")
		salida.write(cut_secuencia+"\n")
	return 0

cutMsa(msa_file,getPosition(query_sequence,dic_msa(msa_file)))
