
##########################################################################
library(tidyr)
library(readr)

count <- read_delim("final.full.count_table", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

contig <- read_delim("final.taxonomy", delim = "\t", 
                     escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)

contig_tax <- contig %>% separate(2,sep = ";",into= c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
contig_tax <- contig_tax[,-8]
head(contig_tax)  


OTU <- read.delim("D:/LindaMetagenomeSeq/MiSOP/final.opti_mcc.0.03.cons.taxonomy")
OTU_tax <- OTU %>% separate(3,sep = ";",into= c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
OTU_tax <- OTU_tax[,-9]

##########################################################################

library(ggplot2)
library(ComplexHeatmap)

count_high <- count[count$total >100,]
count_high2 <- count_high[,-c(1,2)]
count_high2 <- log10(count_high2+1)
count_high2  <- as.matrix(count_high2,ncol=19)
rownames(count_high2) <- count$Representative_Sequence[1:100]

Heatmap(count_high2, row_names_gp = gpar(fontsize = 5))

##########################################################################

OTU_tax_high <- OTU_tax[OTU_tax$Size>100,]

OTU_tax_high$Size <- log10(OTU_tax_high$Size+1)

ggplot(data = OTU_tax_high, aes(x = OTU,y=Size, fill=Class)) + 
  geom_bar(aes(), stat="identity", position="stack") +
  theme(axis.text.x = element_text(angle=45, vjust=.5, hjust=1,size = 5))


##########################################################################

#Rarefaction  Analysis 

rarefaction <- read_delim("final.opti_mcc.groups.ave-std.summary", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)


rarefaction_coverage <- rarefaction[rarefaction$method == "ave",]

ggplot(data=rarefaction_coverage,aes(x=group,y=coverage))+ geom_col()+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1,size = 8))

##########################################################################

pcoa <- read_delim("final.opti_mcc.thetayc.0.03.lt.ave.pcoa.axes", 
                                                      delim = "\t", escape_double = FALSE, 
                                                      trim_ws = TRUE)


ggplot(data=pcoa,aes(x=axis1,y=axis2,label=group))+geom_text(size=3)

##########################################################################
## To export sample phylo Tree

#cluster.split(fasta=final.fasta, count=final.count_table, taxonomy=final.taxonomy, cutoff=0.03)
#make.shared(list=final.opti_mcc.list, count=final.count_table)
#tree.shared(shared=final.opti_mcc.shared)
##########################################################################

library(ggtree)
library(ggplot2)
library(treeio)

mothurtre <- read.newick(file = "final.opti_mcc.thetayc.0.03.tre")

normalize <- function(x, na.rm = TRUE) {
  return((x - min(x)) / (max(x) - min(x)))
}

mothurtre$edge.length <- normalize(mothurtre$edge.length - 0.5)

label_new <- c(mothurtre$tip.label, rep("", 18))

ggtree(mothurtre)


ggtree(mothurtre) +
  geom_text(
    label = label_new,
    size = 4,
    hjust = -0.2,
    color = "firebrick"
  )

ggtree(mothurtre, ladderize = FALSE) +
  geom_text(
    label = label_new,
    size = 4,
    hjust = -0.2,
    color = "firebrick"
  )

ggtree(mothurtre,
       layout = "circular",
       alpha = 0.9,
       size = 3) +
  geom_tiplab(aes(angle = angle), color = 'blue')

label_manual <- seq(1,37)

ggtree(mothurtre,
       layout = "circular",
       alpha = 0.9,
       size = 1) +
  #geom_tiplab(aes(angle = angle), color = 'blue') +
  geom_text(
    aes(label = label_manual),
    size = 6,
    vjust = 0,
    color = "firebrick"
  )

ggtree(mothurtre,
       layout = "fan",
       alpha = 0.9,
       size = 1) +
  geom_text(
    aes(label = label_new),
    size = 3,
    vjust = 0,
    color = "red"
  ) +
  geom_cladelab(
    node = 27,
    label = "This is a clade",
    align = TRUE,
    offset = 0.2,
    geom = 'label',
    fill = 'lightblue'
  )


ggtree(mothurtre,
       layout = "fan",
       alpha = 0.9,
       size = 1) +
  geom_text(
    aes(label = label_new),
    size = 3,
    vjust = 0,
    color = "red"
  ) +
  geom_cladelab(
    node = 5,
    label = "This is a clade",
    align = TRUE,
    offset = 0.2,
    geom = 'label',
    fill = 'lightblue'
  ) +
  geom_hilight(node=30, 
               fill="red")+
  geom_hilight(node=21, 
               fill="green")
  
###########################################################################

    

