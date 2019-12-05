from Bio import SeqIO
import os
import sys

def extractGeneseq(start, end, num):
	GTDBgeneSEQ = []
	for fn in genomeFiles:
		fnaFile = '03_GenesFna/' +  fn + '.na.fasta'
		GTDBResult = '09_ResultGTDB/' + fn + '.aa.txt'

		gtdbHeaders = []
		for gtdb in GTDBgeneSet[start:end]:
			for line in open(GTDBResult).readlines():
				if gtdb in line:
					gtdbHeaders.append(line.split('\t')[0])
					break

		gtdbSeq = '>' + fn + '\n'
		fnaSeqDict =  { str(seq.id):str(seq.seq) for seq in SeqIO.parse(fnaFile, 'fasta') }
		for header in gtdbHeaders:
			gtdbSeq += fnaSeqDict[header][:-1] + '----'

		GTDBgeneSEQ.append(gtdbSeq[:-4] + '\n')
		print(' - GTDB -', str(num), '--', fn)
	open(output + '/GTDB_' + str(num) + '.fasta', 'w').writelines(GTDBgeneSEQ)

workPath = sys.argv[1]
output = sys.argv[2]

os.chdir(workPath)

GTDBgeneSet = [ line.split('\t')[0] for line in open('Database/GTDBmarkerFunc.txt').readlines() ]
genomeFiles = [ fn.split('.f')[0] for fn in os.listdir('01_GenomeFna') ]

num = 1
for i in range(0, len(GTDBgeneSet), 2):
	extractGeneseq(i, i + 2, num)
	num += 1
