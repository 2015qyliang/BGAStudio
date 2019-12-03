
genomeFiles <- list.files('01_GenomeFna')

###############################################################################
if (!dir.exists('02_ProteinsFaa')) {
  dir.create('02_ProteinsFaa')
}
unlink("02_ProteinsFaa", recursive = TRUE)
dir.create('02_ProteinsFaa')

if (!dir.exists('03_GenesFna')) {
  dir.create('03_GenesFna')
}
unlink("03_GenesFna", recursive = TRUE)
dir.create('03_GenesFna')

###############################################################################
for (fn in genomeFiles) {
  fnHeader <- strsplit(fn, '.fasta')[[1]]
  prodigalComd1 <- paste0('./Lib/prodigal.exe -i 01_GenomeFna/', fn, 
                          ' -d 03_GenesFna/', fnHeader, '.na.fasta')
  prodigalComd2 <- paste0(' -a 02_ProteinsFaa/', fnHeader, '.aa.fasta')
  system(paste0(prodigalComd1, prodigalComd2))
}



