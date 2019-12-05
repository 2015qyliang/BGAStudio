from Bio import SeqIO
import os
import sys

def getMaxLength(alignFile):
	seqFile = SeqIO.parse(alignFile, 'fasta')
	lenList = []
	for seq in seqFile:
		lenList.append(len(seq.seq))
	return max(lenList)

def reWriteAlign(alignFile, maxLen):
	newSeqs = []
	seqFile = SeqIO.parse(alignFile, 'fasta')
	for seq in seqFile:
		header = '>' + str(seq.id) + '\n'
		sequence = str(seq.seq) + '-'*(maxLen - len(seq.seq)) + '\n'
		newSeqs.append(header + sequence)
	return newSeqs

def getGenomeGeneSet(alignPath, genomeName):
	geneSeqs = ''
	alignFiles = os.listdir(alignPath)
	for fn in alignFiles:
		alignSeqs = SeqIO.parse(alignPath + fn, 'fasta')
		for seq in alignSeqs:
			if genomeName == str(seq.id):
				geneSeqs += str(seq.seq) + '----'
				break
	genomeSEQ = '>' + genomeName + '\n' + geneSeqs[:-4] + '\n'
	return genomeSEQ

workPath = sys.argv[1]
os.chdir(workPath)
genomeFiles = [ fn.split('.f')[0] for fn in os.listdir('01_GenomeFna') ]

# merge UBCG & GTDB align
UBCG = []
GTDB = []
for fn in genomeFiles:
	print('-- Merge --', fn)
	gnUBCG = getGenomeGeneSet('PhyloAnalysisResult/UBCGalignGeneSet/', fn)
	UBCG.append(gnUBCG)
	gnGTDB = getGenomeGeneSet('PhyloAnalysisResult/GTDBalignGeneSet/', fn)
	GTDB.append(gnGTDB)

open('PhyloAnalysisResult/UBCGalign.fasta', 'w').writelines(UBCG)
open('PhyloAnalysisResult/GTDBalign.fasta', 'w').writelines(GTDB)

# add tail
seqs = reWriteAlign('PhyloAnalysisResult/UBCGalign.fasta', getMaxLength('PhyloAnalysisResult/UBCGalign.fasta'))
open('PhyloAnalysisResult/UBCGalign.fasta', 'w').writelines(seqs)
seqs = reWriteAlign('PhyloAnalysisResult/GTDBalign.fasta', getMaxLength('PhyloAnalysisResult/GTDBalign.fasta'))
open('PhyloAnalysisResult/GTDBalign.fasta', 'w').writelines(seqs)
