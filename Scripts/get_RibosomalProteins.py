from Bio import SeqIO
import os
import sys

workPath = sys.argv[1]
output = sys.argv[2]

os.chdir(workPath)

ribosomalGeneSet = [ line.strip() for line in open('Database/RibosomalProteinSet.txt').readlines() ]
genomeFiles = [ fn.split('.f')[0] for fn in os.listdir('01_GenomeFna') ]
ribosomalProteinsSEQ = []

for fn in genomeFiles:
	faaFile = '02_ProteinsFaa/' +  fn + '.aa.fasta'
	UniprotResult = '05_ResultUniprot/' + fn + '.aa.txt'

	ribosomalHeaders = []
	for ribosomal in ribosomalGeneSet:
		for line in open(UniprotResult).readlines():
			if ribosomal in line:
				ribosomalHeaders.append(line.split('\t')[0])
				break

	ribosomalSeq = '>' + fn + '\n'
	faaSeqDict =  { str(seq.id):str(seq.seq) for seq in SeqIO.parse(faaFile, 'fasta') }
	for header in ribosomalHeaders:
		ribosomalSeq += faaSeqDict[header][:-1] + '----'

	ribosomalProteinsSEQ.append(ribosomalSeq[:-4] + '\n')
	print(' -- Ribosomal -- ', fn, '--- ' )

open(output + '/ribosomalProteins.aa.fasta', 'w').writelines(ribosomalProteinsSEQ)
