
genomeFiles <- list.files('02_ProteinsFaa')
idCutoff <- 35 # the value could be reset
queryCoverage <- 50 # 50 or 70

diamondProcess <- c(FALSE, # TRUE or FALSE -- Uniprot
                    FALSE # TRUE or FALSE -- GTDBmarker
                    ) 

###############################################################################

if (!dir.exists('PhyloAnalysisResult')) {
  dir.create('PhyloAnalysisResult')
}
unlink("PhyloAnalysisResult", recursive = TRUE)
dir.create('PhyloAnalysisResult')


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
# summary rpoB
system( paste0('python Scripts/get_rpoB.py ', getwd(), '  PhyloAnalysisResult') )

# summary ribosomal proteins
system( paste0('python Scripts/get_RibosomalProteins.py ', getwd(), '  PhyloAnalysisResult') )

# summary ribosomal proteins
system( paste0('python Scripts/get_UBCG.py ', getwd(), '  PhyloAnalysisResult') )

# summary GTDB
system( paste0('python Scripts/get_GTDB.py ', getwd(), '  PhyloAnalysisResult') )

