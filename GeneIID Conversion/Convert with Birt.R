#install the libraries
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("clusterProfiler")

genelist <- c('ENSG00000198938.2',
              'ENSG00000198712.1',
              'ENSG00000198804.2')

#Remove all numbers after the dot
genelist <- sub("[.][0-9]*", "", genelist)

library(clusterProfiler)
library(org.Hs.eg.db)
new_genelist <- bitr(
  genelist,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db,
  drop = TRUE
)
new_genelist

# ENSEMBL             ENTREZID
# ENSG00000198938     4514
# ENSG00000198712     4513
# ENSG00000198804     4512

# use keytypes(org.Hs.eg.db)  to get
#all the kyes type

keytypes(org.Hs.eg.db)
#[1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"
#[8] "EVIDENCE"     "EVIDENCEALL"  "GENENAME"     "GENETYPE"     "GO"           "GOALL"        "IPI"
#[15] "MAP"          "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"
#[22] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"       "UNIPROT"     