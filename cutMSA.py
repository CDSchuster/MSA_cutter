from Bio import SeqIO
import argparse

def getQuerySeq(query, rev):
    """Toma la secuencia de un query y la devuelve como string"""
    q_seq=str(list(SeqIO.parse(query, "fasta"))[0].seq)
    if rev: q_seq=q_seq[::-1]
    return q_seq

def getPosition(query_seq, ref):
    """recibe un query y una refernecia en el MSA, las alinea y
       devuelve posiciones de inicio y fin del alineamiento
       del query contra la referencia"""
    query_seq_index=0
    query_seq_length=len(query_seq)
    largo_seqs=len(ref)
    error, i = 0, 0
    while i<largo_seqs and query_seq_index<query_seq_length:
        if ref[query_seq_index+i]!=query_seq[query_seq_index]:
            error+=1
            query_seq_index=query_seq_index+1
        elif ref[query_seq_index+i]==query_seq[query_seq_index]:
            query_seq_index=query_seq_index+1
        if error>10 or query_seq_index==error:
            query_seq_index=0
            error=0
            i+=1
    final=i+query_seq_index
    return [i, final]

def cutMSA(MSA, query, outFile, rev):
    """Recibe el MSA, el archivo salida y el query, ejecuta las demas funciones,
    obtiene las posiciones de alineamiento y recorta el MSA a esas posiciones.
    El resultado lo escribe en la salida"""
    query_seq=getQuerySeq(query, rev)
    fasta=list(SeqIO.parse(MSA, "fasta"))
    ref=fasta[0].seq
    lims=getPosition(query_seq, ref)
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
    parser.add_argument("-rev", nargs='+', type=bool)
    parser.add_argument("-o", nargs='+')
    return parser.parse_args()

def main():
    """Funcion donde se ejcutan todas las demas funcionas en serie"""
    args=parse_arguments()
    cutMSA(args.msa[0], args.q[0], args.o[0], args.rev[0])
    return 0

if __name__=='__main__':
    main()
