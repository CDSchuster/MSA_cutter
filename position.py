#!/usr/bin/env python3

from Bio import SeqIO

def dic_msa(file_msa):
	dic={}
	for record in SeqIO.parse(file_msa, 'fasta'):
		dic[record.description]=record.seq
	return dic

def getPosition(query_seq, dic_msa):

	# Importante: Toma una referencia de todo el msa para buscar el inicio y final de la query
	secuencia=dic_msa[list(dic_msa.keys())[0]]
	query_seq_index=0
	query_seq_length=len(query_seq)
	largo_seqs=len(secuencia)
	error=-1

	for i in range(largo_seqs):
		if query_seq_index < query_seq_length and secuencia[i] != query_seq[query_seq_index]:
			error=error+1
		if  query_seq_index < query_seq_length and secuencia[i] == query_seq[query_seq_index]:
			query_seq_index=query_seq_index+1

		if error>3 or query_seq_index == error:
			query_seq_index=0
			error=-1
		if query_seq_length == query_seq_index:
				inicio= i - query_seq_length - error 
				final= i + 1
				query_seq_length=0

	print([inicio, final])

	return [inicio, final]
