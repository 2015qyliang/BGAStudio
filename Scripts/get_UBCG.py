from Bio import SeqIO
import os
import sys

workPath = sys.argv[1]
output = sys.argv[2]

os.chdir(workPath)

UBCGgeneSet = [ line.strip() for line in open('Database/UBCGgeneSet.txt').readlines() ]
genomeFiles = [ fn.split('.f')[0] for fn in os.listdir('01_GenomeFna') ]
UBCGgeneSEQ = []

for fn in genomeFiles:
	fnaFile = '03_GenesFna/' +  fn + '.na.fasta'
	UniprotResult = '05_ResultUniprot/' + fn + '.aa.txt'

	ubcgHeaders = []
	for ubcg in UBCGgeneSet:
		for line in open(UniprotResult).readlines():
			if ubcg + '	' in line:
				ubcgHeaders.append(line.split('\t')[0])
				break

	ubcgSeq = '>' + fn + '\n'
	fnaSeqDict =  { str(seq.id):str(seq.seq) for seq in SeqIO.parse(fnaFile, 'fasta') }
	for header in ubcgHeaders:
		ubcgSeq += fnaSeqDict[header][:-1] + '----'

	UBCGgeneSEQ.append(ubcgSeq[:-4] + '\n')
	print(' -- UBCG -- ', fn, '--- ' )

open(output + '/UBCG.fasta', 'w').writelines(UBCGgeneSEQ)
