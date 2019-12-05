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

workPath = sys.argv[1]
os.chdir(workPath)

UBCGalignFiles = [ 'PhyloAnalysisResult/UBCGalignGeneSet/' + fn for fn in os.listdir('PhyloAnalysisResult/UBCGalignGeneSet/') ]
GTDBalignFiles = [ 'PhyloAnalysisResult/GTDBalignGeneSet/' + fn for fn in os.listdir('PhyloAnalysisResult/GTDBalignGeneSet/') ]

for fn in UBCGalignFiles:
	seqs = reWriteAlign(fn, getMaxLength(fn))
	open(fn, 'w').writelines(seqs)
for fn in GTDBalignFiles:
	seqs = reWriteAlign(fn, getMaxLength(fn))
	open(fn, 'w').writelines(seqs)

