#!/usr/bin/env python3

from Bio import SeqIO

# Revierte un string de secuencia
def revert(seq):
	seq=seq[::-1]
	return seq

# Le un archivo fasta y revierte las secuencias que deben ser revertidas
# Como sabe cual revertir: en el nombre al final va a haber una letra F o R. R debe ser revertida.
# Devuelve el mismo formato fasta
def revertFile(file):
	records=[]
	for record in SeqIO.parse(file, 'fasta'):
		reverse_check=record.description.split(" ")[-1]
		if reverse_check == "R":
			sequence=revert(record.seq)
		else:
			sequence=record.seq
		record.seq=sequence
		records.append(record)
	file_out=file.split(".")[:-1][0]+".revert.fasta"

	SeqIO.write(records, file_out, "fasta")
	print("Reverted "+file+" to "+file_out)
	return 0