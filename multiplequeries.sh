for file in ./*query*; do python3 cutMSA.py -msa subsampled_alignment.fasta -q $file -o $file.out; done
