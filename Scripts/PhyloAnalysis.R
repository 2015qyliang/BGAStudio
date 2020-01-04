
genomeFiles <- list.files('02_ProteinsFaa')
idCutoff <- 35 # the value could be reset
queryCoverage <- 50 # 50 or 70

diamondProcess <- c(FALSE, # TRUE or FALSE -- Uniprot
                    FALSE # TRUE or FALSE -- GTDBmarker
                    )

###############################################################################
if (diamondProcess[1]) {

  if (!dir.exists('05_ResultUniprot')) {
    dir.create('05_ResultUniprot')
  }
  unlink("05_ResultUniprot", recursive = TRUE)
  dir.create('05_ResultUniprot')

  dbtype <- "Uniprot"
  for (fn in genomeFiles) {
    fnHeader <- strsplit(fn, '.fasta')[[1]]
    blast.comm1 <- paste0('./Lib/diamond.exe blastp -q 02_ProteinsFaa/', fn)
    blast.comm2 <- paste0(' -d Database/', dbtype)
    blast.comm3 <- paste0(' --max-target-seqs 1 -e 1e-8 --id ', idCutoff,
                          ' -o 05_ResultUniprot/', fnHeader, '.txt',
                          ' --tmpdir 04_DiamondTmpDir/',
                          ' --query-cover ', queryCoverage)
    system(paste0(blast.comm1, blast.comm2, blast.comm3))
  }
}

###############################################################################
if (diamondProcess[2]) {

  if (!dir.exists('09_ResultGTDB')) {
    dir.create('09_ResultGTDB')
  }
  unlink("09_ResultGTDB", recursive = TRUE)
  dir.create('09_ResultGTDB')

  dbtype <- "GTDBmarker"
  for (fn in genomeFiles) {
    fnHeader <- strsplit(fn, '.fasta')[[1]]
    blast.comm1 <- paste0('./Lib/diamond.exe blastp -q 02_ProteinsFaa/', fn)
    blast.comm2 <- paste0(' -d Database/', dbtype)
    blast.comm3 <- paste0(' --max-target-seqs 1 -e 1e-8 --id ', idCutoff,
                          ' -o 09_ResultGTDB/', fnHeader, '.txt',
                          ' --tmpdir 04_DiamondTmpDir/',
                          ' --query-cover ', queryCoverage)
    system(paste0(blast.comm1, blast.comm2, blast.comm3))
  }
}

###############################################################################
if (!dir.exists('PhyloAnalysisResult')) {
  dir.create('PhyloAnalysisResult')
}
unlink("PhyloAnalysisResult", recursive = TRUE)
dir.create('PhyloAnalysisResult')

###############################################################################
# summary rpoB
system( paste0('python Scripts/get_rpoB.py ', getwd(), '  PhyloAnalysisResult') )

# summary ribosomal proteins
system( paste0('python Scripts/get_RibosomalProteins.py ', getwd(), '  PhyloAnalysisResult') )

# summary UBCG
system( paste0('python Scripts/get_UBCG.py ', getwd(), '  PhyloAnalysisResult') )

# summary GTDB
system( paste0('python Scripts/get_GTDB.py ', getwd(), '  PhyloAnalysisResult') )


###############################################################################
#####################        Sequences   Alignment        #####################
muscle <-' ./Lib/muscle3.8.31_i86win32.exe -maxiters 2 -in PhyloAnalysisResult/'

system( paste0(muscle, 'ribosomalProteins.aa.fasta ',
               ' -out  PhyloAnalysisResult/ribosomalProteins.aa_align.fasta') )
system( paste0(muscle, 'rpoB.aa.fasta ',
               ' -out  PhyloAnalysisResult/rpoB.aa_align.fasta') )

# summary UBCG
dir.create('PhyloAnalysisResult/UBCGrawGeneSet')
dir.create('PhyloAnalysisResult/UBCGalignGeneSet')
UBCGfiles <- list.files('PhyloAnalysisResult/',pattern = 'UBCG_')
for (fn in UBCGfiles) {
  print('====================================================================================')
  print(paste0('>>>>>>>>>>>>>>>> sequences alignment by Muscle ===============> ', fn) )
  fileHeader <- strsplit(fn, '.fasta')[[1]][1]
  system( paste0(muscle, fn, ' -out  PhyloAnalysisResult/', fileHeader, '_align.fasta') )
  file.copy(paste0('PhyloAnalysisResult/', fn), 
            paste0('PhyloAnalysisResult/UBCGrawGeneSet/', fn))
  file.copy(paste0('PhyloAnalysisResult/', fileHeader, '_align.fasta'), 
            paste0('PhyloAnalysisResult/UBCGalignGeneSet/', fileHeader, '_align.fasta'))
  unlink(paste0('PhyloAnalysisResult/', fn), recursive = T)
  unlink(paste0('PhyloAnalysisResult/', fileHeader, '_align.fasta'), recursive = T)
}

# summary GTDB
dir.create('PhyloAnalysisResult/GTDBrawGeneSet')
dir.create('PhyloAnalysisResult/GTDBalignGeneSet')
GTDBfiles <- list.files('PhyloAnalysisResult/',pattern = 'GTDB_')
for (fn in GTDBfiles) {
  print('====================================================================================')
  print(paste0('>>>>>>>>>>>>>>>> sequences alignment by Muscle ===============> ', fn) )
  fileHeader <- strsplit(fn, '.fasta')[[1]][1]
  system( paste0(muscle, fn, ' -out  PhyloAnalysisResult/', fileHeader, '_align.fasta') )
  file.copy(paste0('PhyloAnalysisResult/', fn), 
            paste0('PhyloAnalysisResult/UBCGrawGeneSet/', fn))
  file.copy(paste0('PhyloAnalysisResult/', fileHeader, '_align.fasta'), 
            paste0('PhyloAnalysisResult/GTDBalignGeneSet/', fileHeader, '_align.fasta'))
  unlink(paste0('PhyloAnalysisResult/', fn), recursive = T)
  unlink(paste0('PhyloAnalysisResult/', fileHeader, '_align.fasta'), recursive = T)
}

###############################################################################
# add tail to obtain equal length
system( paste0('python Scripts/addTail.py ', getwd()) )


###############################################################################
#####################        Merge       Alignment        #####################
system( paste0('python Scripts/mergeAlignUBCG&GTDB.py ', getwd()) )

###############################################################################
# trim alignment  -in Align.fasta -out AlignTrimed.fasta -automated1
trim <- './Lib/trimAl/bin/trimal.exe  -automated1 -in '
system( paste0(trim, 'PhyloAnalysisResult/rpoB.aa_align.fasta -out PhyloAnalysisResult/rpoB_AlignTrimed.fasta') )
system( paste0(trim, 'PhyloAnalysisResult/ribosomalProteins.aa_align.fasta -out PhyloAnalysisResult/ribosomalP_AlignTrimed.fasta') )
system( paste0(trim, 'PhyloAnalysisResult/UBCGalign.fasta -out PhyloAnalysisResult/UBCG_AlignTrimed.fasta') )
system( paste0(trim, 'PhyloAnalysisResult/GTDBalign.fasta -out PhyloAnalysisResult/GTDB_AlignTrimed.fasta') )

###############################################################################
dir.create('PhyloAnalysisResult/rawAlignFiles')
file.copy('PhyloAnalysisResult/rpoB.aa_align.fasta', 
          'PhyloAnalysisResult/rawAlignFiles/rpoB.aa_align.fasta')
unlink('PhyloAnalysisResult/rpoB.aa_align.fasta', recursive = T)
file.copy('PhyloAnalysisResult/ribosomalProteins.aa_align.fasta', 
          'PhyloAnalysisResult/rawAlignFiles/ribosomalProteins.aa_align.fasta')
unlink('PhyloAnalysisResult/ribosomalProteins.aa_align.fasta', recursive = T)
file.copy('PhyloAnalysisResult/UBCGalign.fasta', 
          'PhyloAnalysisResult/rawAlignFiles/UBCGalign.fasta')
unlink('PhyloAnalysisResult/UBCGalign.fasta', recursive = T)
file.copy('PhyloAnalysisResult/GTDBalign.fasta', 
          'PhyloAnalysisResult/rawAlignFiles/GTDBalign.fasta')
unlink('PhyloAnalysisResult/GTDBalign.fasta', recursive = T)
file.copy('PhyloAnalysisResult/rpoB.aa.fasta', 
          'PhyloAnalysisResult/rawAlignFiles/rpoB.aa.fasta')
unlink('PhyloAnalysisResult/rpoB.aa.fasta', recursive = T)
file.copy('PhyloAnalysisResult/ribosomalProteins.aa_align.fasta', 
          'PhyloAnalysisResult/rawAlignFiles/ribosomalProteins.aa.fasta')
unlink('PhyloAnalysisResult/ribosomalProteins.aa.fasta', recursive = T)

