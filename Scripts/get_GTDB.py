from Bio import SeqIO
import os
import sys

workPath = sys.argv[1]
output = sys.argv[2]

os.chdir(workPath)

GTDBgeneSet = [ line.split('\t')[0] for line in open('Database/GTDBmarkerFunc.txt').readlines() ]
genomeFiles = [ fn.split('.f')[0] for fn in os.listdir('01_GenomeFna') ]
GTDBgeneSEQ = []

for fn in genomeFiles:
	fnaFile = '03_GenesFna/' +  fn + '.na.fasta'
	GTDBResult = '09_ResultGTDB/' + fn + '.aa.txt'

	gtdbHeaders = []
	for gtdb in GTDBgeneSet:
		for line in open(GTDBResult).readlines():
			if gtdb in line:
				gtdbHeaders.append(line.split('\t')[0])
				break

	gtdbSeq = '>' + fn + '\n'
	fnaSeqDict =  { str(seq.id):str(seq.seq) for seq in SeqIO.parse(fnaFile, 'fasta') }
	for header in gtdbHeaders:
		gtdbSeq += fnaSeqDict[header][:-1] + '----'

	GTDBgeneSEQ.append(gtdbSeq[:-4] + '\n')
	print(' -- GTDB -- ', fn, '--- ' )

open(output + '/GTDB.fasta', 'w').writelines(GTDBgeneSEQ)
