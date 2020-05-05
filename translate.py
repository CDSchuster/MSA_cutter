from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

class translator(object):

    def __init__(self, inFile, outFile):
        fasta=self.__loadFasta__(inFile)
        records=self.__translate_MSA__(fasta)
        self.__writeFasta__(records, outFile)

    def __loadFasta__(self, fasFile):
        fasta=list(SeqIO.parse(fasFile, "fasta"))
        return fasta

    def __translate_seq__(self, seq):
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

    def __translate_MSA__(self, fasta):
        records=[]
        for rec in fasta:
            name=rec.description.replace(" ", "")
            aa=Seq(self.__translate_seq__(str(rec.seq)), IUPAC.protein)
            newRec=SeqRecord(aa, id=name, description="")
            records.append(newRec)
        return records

    def __writeFasta__(self, records, outFile):
        SeqIO.write(records, outFile, "fasta")
        return 0
