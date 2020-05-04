from Bio import SeqIO
import argparse

def main(MSA, inf, sup, outFile):
    fasta=list(SeqIO.parse(MSA, "fasta"))
    outf=open(outFile, "w")
    for rec in fasta:
        print(">"+rec.description+"\n"+rec.seq[inf-1:sup], file=outf)
    outf.close()

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-msa", nargs='+')
    parser.add_argument("-inf", nargs='?', type=int)
    parser.add_argument("-sup", nargs='?', type=int)
    parser.add_argument("-o", nargs='?')
    args=parser.parse_args()
    main(args.msa[0], args.inf, args.sup, args.o)
