from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import Gapped, generic_dna, IUPAC
import pandas as pd

def loadFasta(fasFile):
    fasta=list(SeqIO.parse(fasFile, "fasta"))
    return fasta

def translate_seq(seq):
    seq=seq.upper()
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    protein =""
    nts=["A", "T", "C", "G"]
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            if codon[0] in nts and codon[1] in nts and codon[2] in nts:
                protein+= table[codon]
            else: protein+="X"
    return protein

def translate_MSA(fasta):
    records=[]
    for rec in fasta:
        name=rec.description.replace(" ", "")
        aa=Seq(translate_seq(str(rec.seq)), IUPAC.protein)
        newRec=SeqRecord(aa, id=name, description="")
        records.append(newRec)
    return records

def writeFasta(records, outFile):
    SeqIO.write(records, outFile, "fasta")
    return 0

def main(inFile, outFile):
    fasta=loadFasta(inFile)
    records=translate_MSA(fasta)
    writeFasta(records, outFile)
    return 0

if __name__=='__main__':
    main(inFile, outFile)
