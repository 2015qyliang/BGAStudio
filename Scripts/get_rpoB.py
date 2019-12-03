from Bio import SeqIO
import os
import sys

workPath = sys.argv[1]
output = sys.argv[2]

os.chdir(workPath)

genomeFiles = [ fn.split('.f')[0] for fn in os.listdir('01_GenomeFna') ]
rpoBSEQ = []

for fn in genomeFiles:
	faaFile = '02_ProteinsFaa/' +  fn + '.aa.fasta'
	UniprotResult = '05_ResultUniprot/' + fn + '.aa.txt'

	rpoBseqHeader = [ line.split('\t')[0] for line in open(UniprotResult).readlines() if "GN=rpoB	" in line ][0]
	faaSeq = SeqIO.parse(faaFile, 'fasta')
	for seq in faaSeq:
		if str(seq.id) == rpoBseqHeader:
			rpoBseq = '>' + fn + '\n' + str(seq.seq)[:-1] + '\n'
			rpoBSEQ.append(rpoBseq)
			print(' -- rpoB -- ', fn, '--- ' )
			break

open(output + '/rpoB.aa.fasta', 'w').writelines(rpoBSEQ)
