from Bio import SeqIO
import argparse

def getQuerySeq(query):
    """Toma la secuencia de un query y la devuelve como string"""
    q_seq=str(list(SeqIO.parse(query, "fasta"))[0].seq)
    return q_seq

def MSA_to_dict(file_msa):
    """Toma un MSA y lo convierte en diccionario"""
    dic={}
    for record in SeqIO.parse(file_msa, 'fasta'):
        dic[record.description]=record.seq
    return dic

def getPosition(query_seq, dic_msa):
    """recibe seq y dicc de MSA
       Importante: Toma una referencia de todo el msa
       para buscar el inicio y final de la query"""
    ref=dic_msa[list(dic_msa.keys())[0]] #referencia es primer seq
    query_seq_index=0
    query_seq_length=len(query_seq)
    largo_seqs=len(ref)
    error=-1
    for i in range(largo_seqs):
        print(query_seq_index, error, query_seq[query_seq_index], ref[i], file=aa)
        if query_seq_index<query_seq_length and ref[i]!=query_seq[query_seq_index]:
            error=error+1
        if  query_seq_index<query_seq_length and ref[i]==query_seq[query_seq_index]:
            query_seq_index=query_seq_index+1
        if error>3 or query_seq_index==error:
            query_seq_index=0
            error=-1
        if query_seq_length==query_seq_index:
            inicio=i-query_seq_length-error
            final=i+1
            query_seq_length=0
    return [inicio, final]

def cutMSA(MSA, lims, outFile):
    """Recibe el MSA original, el archivo salida y las posiciones limite
    y recorta el MSA a esas posiciones. El resultado lo escribe en la salida"""
    fasta=list(SeqIO.parse(MSA, "fasta"))
    outf=open(outFile, "w")
    for rec in fasta:
        print(">"+rec.description+"\n"+rec.seq[lims[0]:lims[1]], file=outf)
    outf.close()
    return 0

def parse_arguments():
    """Parsea los argumentos de entrada del script"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-msa", nargs='+')
    parser.add_argument("-q", nargs='+')
    parser.add_argument("-o", nargs='+')
    return parser.parse_args()

def main():
    """Funcion donde se ejcutan todas las demas funcionas en serie"""
    args=parse_arguments()
    dict_MSA=MSA_to_dict(args.msa[0])
    query_seq=getQuerySeq(args.q[0])
    pos=getPosition(query_seq, dict_MSA)
    cutMSA(args.msa[0], pos, args.o[0])
    return 0

if __name__=='__main__':
    main()
