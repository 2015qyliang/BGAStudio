
genomeFiles <- list.files('02_ProteinsFaa')
idCutoff <- 35 # the value could be reset

diamondProcess <- c(
  FALSE, # 1 - Uniprot -- TRUE or FALSE
  FALSE, # 2 - OrthoMclOG5 -- TRUE or FALSE
  FALSE, # 3 - COG -- TRUE or FALSE
  TRUE # 4 - CAZy -- TRUE or FALSE
)

###############################################################################
if (!dir.exists('05_ResultUniprot')) {
  dir.create('05_ResultUniprot')
}
unlink("05_ResultUniprot", recursive = TRUE)
dir.create('05_ResultUniprot')

if (!dir.exists('06_ResultOrthoMclOG5')) {
  dir.create('06_ResultOrthoMclOG5')
}
unlink("06_ResultOrthoMclOG5", recursive = TRUE)
dir.create('06_ResultOrthoMclOG5')

if (!dir.exists('07_ResultCOG')) {
  dir.create('07_ResultCOG')
}
unlink("07_ResultCOG", recursive = TRUE)
dir.create('07_ResultCOG')

if (!dir.exists('08_ResultCAZy')) {
  dir.create('08_ResultCAZy')
}
unlink("08_ResultCAZy", recursive = TRUE)
dir.create('08_ResultCAZy')

###############################################################################
if (diamondProcess[1]) {
  dbtype <- "Uniprot"
  for (fn in genomeFiles) {
    fnHeader <- strsplit(fn, '.fasta')[[1]]
    blast.comm1 <- paste0('./Lib/diamond.exe blastp -q 02_ProteinsFaa/', fn)
    blast.comm2 <- paste0(' -d Database/', dbtype)
    blast.comm3 <- paste0(' --max-target-seqs 1 -e 1e-5 --id ', idCutoff,
                          ' -o 05_ResultUniprot/', fnHeader, '.txt',
                          ' --tmpdir 04_DiamondTmpDir/')
    system(paste0(blast.comm1, blast.comm2, blast.comm3))
  }
}

###############################################################################
if (diamondProcess[2]) {
  dbtype <- "OrthoMclOG5"
  for (fn in genomeFiles) {
    fnHeader <- strsplit(fn, '.fasta')[[1]]
    blast.comm1 <- paste0('./Lib/diamond.exe blastp -q 02_ProteinsFaa/', fn)
    blast.comm2 <- paste0(' -d Database/', dbtype)
    blast.comm3 <- paste0(' --max-target-seqs 1 -e 1e-5 --id ', idCutoff,
                          ' -o 06_ResultOrthoMclOG5/', fnHeader, '.txt',
                          ' --tmpdir 04_DiamondTmpDir/')
    system(paste0(blast.comm1, blast.comm2, blast.comm3))
  }
}

###############################################################################
if (diamondProcess[3]) {
  dbtype <- "COG2014"
  for (fn in genomeFiles) {
    fnHeader <- strsplit(fn, '.fasta')[[1]]
    blast.comm1 <- paste0('./Lib/diamond.exe blastp -q 02_ProteinsFaa/', fn)
    blast.comm2 <- paste0(' -d Database/', dbtype)
    blast.comm3 <- paste0(' --max-target-seqs 1 -e 1e-5 --id ', idCutoff,
                          ' -o 07_ResultCOG/', fnHeader, '.txt',
                          ' --tmpdir 04_DiamondTmpDir/')
    system(paste0(blast.comm1, blast.comm2, blast.comm3))
  }
}

###############################################################################
if (diamondProcess[4]) {
  dbtype <- "CAZyDB"
  for (fn in genomeFiles) {
    fnHeader <- strsplit(fn, '.fasta')[[1]]
    blast.comm1 <- paste0('./Lib/diamond.exe blastp -q 02_ProteinsFaa/', fn)
    blast.comm2 <- paste0(' -d Database/', dbtype)
    blast.comm3 <- paste0(' --max-target-seqs 1 -e 1e-5 --id ', idCutoff,
                          ' -o 08_ResultCAZy/', fnHeader, '.txt',
                          ' --tmpdir 04_DiamondTmpDir/')
    system(paste0(blast.comm1, blast.comm2, blast.comm3))
  }
}

###############################################################################
