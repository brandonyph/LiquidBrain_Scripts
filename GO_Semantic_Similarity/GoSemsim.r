
library(GOSemSim)
hsGO <- godata('org.Hs.eg.db', ont="MF")
goSim("GO:0004022", "GO:0005515", semData=hsGO, measure="Wang")
      
go1 = c("GO:0004022","GO:0004024","GO:0004174")
go2 = c("GO:0009055","GO:0005515")
mgoSim(go1, go2, semData=hsGO, measure="Wang", combine=NULL)
     
hsGO2 <- godata('org.Hs.eg.db', keytype = "SYMBOL", ont="MF", computeIC=FALSE)   
genes <- c("CDC45", "MCM10", "CDC20", "NMU", "MMP1")
mgeneSim(genes, semData=hsGO2, measure="Wang", combine="BMA")
     


gs1 <- c("835", "5261","241", "994", "514", "533")
gs2 <- c("578","582", "400", "409", "411")
clusterSim(gs1, gs2, semData=hsGO, measure="Wang", combine="BMA")
