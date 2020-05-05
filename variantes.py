#! python3
from Bio import SeqIO

def countLetters(lista, divisor, no_letter, proportion_fiter):
	uniques=set(lista)
	out_list=[]
	for element in uniques:
		count = lista.count(element)
		proportion=count*100/divisor
		if proportion>proportion_fiter and element != no_letter:
			out_list.append(element + ": "+ "{0:.3f}".format(proportion)) 
	return out_list

def variantes(file, type_file_prot, proportion_fiter):
	
	# Nombre del archivo de salida
	salida_file=open("variantes_file.txt","w")

	if type_file_prot:
		no_letter="X"
	else:
		no_letter="N"

	sequences={}
	cantidad_seq=0
	for record in SeqIO.parse(file, 'fasta'):
		sequences[record.description]=record.seq
		largo_seqs=len(record.seq)
		cantidad_seq+=1

	column_letters=[]
	column_dif=0
	for i in range(largo_seqs):
		column_letters=[]
		for seq in sequences:
			secuencia=sequences[seq]
			column_letters.append(secuencia[i])
		if ((no_letter in column_letters and len(set(column_letters))>2) 
			or (no_letter not in column_letters and len(set(column_letters))>1)):

			countLetters_list=countLetters(column_letters,cantidad_seq,no_letter, proportion_fiter)

			if len(countLetters_list)>1:
				salida_file.write(">Posicion: " + str(i+1)+"\n")
				countLetters_list='\n'.join(countLetters_list)
				salida_file.write(countLetters_list+"\n")
				column_dif=column_dif+1
			
	salida_file.write("Cantidad de columnas con variantes: "+ str(column_dif)+"\n")
	salida_file.write("cantidad de secuencias: "+str(cantidad_seq)+"\n")
	salida_file.write("Largo de las secuencias: "+str(largo_seqs))

	salida_file.close()

	return 0

# La funcion recibe, un archivo de entrada(msa comun o ya cortado) 
# Un termino booleano para ver si es de proteinas(True) o nucleotidos(False)
# Filtro de proporcion de variantes. Solo se veran las variantes mayores a ese valor
# Escribe un archivo de salida con los resultados
variantes("/home/agustin/workspace/covid-19/msa_colum_extract/protein/S_prot.fasta",True,0.01)