library("Biostrings")
library(vioplot)
library(stringr)
small<- readDNAStringSet("uniprot-40s+ribosomal+protein-filtered-reviewed_yes+AND+organism__Homo+sap--.fasta")
big <- readDNAStringSet("uniprot-60s+ribosomal+protein+AND+reviewed_yes+AND+organism__Homo+sapiens+--.fasta")
proteome<-readDNAStringSet("uniprot-reviewed_yes+AND+organism__Homo+sapiens+(Human)+[9606]_.fasta")
dss2df <- function(dss) data.frame(width=width(dss), seq=as.character(dss), names=names(dss))
small<-dss2df(small)
big<-dss2df(big)
proteome<-dss2df(proteome)

small<-small[grep("40S ribosomal protein",small$names),]
big<-big[grep("60S ribosomal protein",big$names),]

vioplot(
  c(str_count(as.character(proteome$seq), "K")/str_count(as.character(proteome$seq))*100),
  c(str_count(as.character(big$seq), "K")/str_count(as.character(big$seq))*100),
  c(str_count(as.character(small$seq), "K")/str_count(as.character(small$seq))*100),
  c(str_count(as.character(proteome$seq), "R")/str_count(as.character(proteome$seq))*100),
  c(str_count(as.character(big$seq), "R")/str_count(as.character(big$seq))*100),
  c(str_count(as.character(small$seq), "R")/str_count(as.character(small$seq))*100),
  ylab="proportion of residues", names=c("\"K\" proteome","\"K\" RPL","\"K\" RPS","\"R\" proteome","\"R\" RPL","\"R\" RPS"),las=2,
  col=c("blue","red","green","blue","red","green"),cex.axis=0.7)

summary(str_count(as.character(big$seq), "R")/str_count(as.character(big$seq))*100)
summary(str_count(as.character(small$seq), "R")/str_count(as.character(small$seq))*100)
summary(str_count(as.character(big$seq), "K")/str_count(as.character(big$seq))*100)
summary(str_count(as.character(small$seq), "K")/str_count(as.character(small$seq))*100)
summary(str_count(as.character(proteome$seq), "K")/str_count(as.character(proteome$seq))*100)

#rps2
vioplot(
  c(c(str_count(as.character(small$seq)[2], "K")/str_count(as.character(small$seq)[2])*100),
    c(str_count(as.character(small$seq)[grep("RPS20",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS20",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS21",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS21",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS6",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS6",small$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPLP2",big$names)], "K")/str_count(as.character(big$seq)[grep("RPLP2",big$names)])*100),
    
    
    c(str_count(as.character(small$seq)[grep("RPS3",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS3",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS9",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS9",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS26",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS26",small$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL17",big$names)], "K")/str_count(as.character(big$seq)[grep("RPL17",big$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS17L",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS17L",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS28",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS28",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPSA",small$names)], "K")/str_count(as.character(small$seq)[grep("RPSA",small$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL10",big$names)], "K")/str_count(as.character(big$seq)[grep("RPL10",big$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS27A",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS27A",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS15A",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS15A",small$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL15",big$names)], "K")/str_count(as.character(big$seq)[grep("RPL15",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL37A",big$names)], "K")/str_count(as.character(big$seq)[grep("RPL37A",big$names)])*100)),
  
  
  c(c(str_count(as.character(small$seq)[grep("RPS7",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS7",small$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL38",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL38",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL19",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL19",big$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS4X",small$names)], "K")/ str_count(as.character(small$seq)[grep("RPS4X",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS16",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS16",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS14",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS14",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS11",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS11",small$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL21",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL21",big$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS10",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS10",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS23",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS23",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS19",small$names)], "K")/ str_count(as.character(small$seq)[grep("RPS19",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS12",small$names)], "K")/ str_count(as.character(small$seq)[grep("RPS12",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS5",small$names)], "K")/ str_count(as.character(small$seq)[grep("RPS5",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS3A",small$names)], "K")/ str_count(as.character(small$seq)[grep("RPS3A",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS25",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS25",small$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL23",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL23",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL26",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL26",big$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS8",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS8",small$names)])*100)),
  
  c(
    c(str_count(as.character(big$seq)[grep("RPL28",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL28",big$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS27L",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS27L",small$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL35A",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL35A",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL32",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL32",big$names)])*100),
    
    c(str_count(as.character(big$seq)[grep("RPL29",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL29",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL22",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL22",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL31",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL31",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL30",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL30",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL8",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL8",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL12",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL12",big$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS13",small$names)], "K")/ str_count(as.character(small$seq)[grep("RPS13",small$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL27",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL27",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL13A",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL13A",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL9",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL9",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL18",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL18",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL7",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL7",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL35AL",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL35AL",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL14",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL14",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL34",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL34",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL13",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL13",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL35",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL35",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL36",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL36",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL15",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL15",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL27A",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL27A",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL24",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL24",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL3",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL3",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL11",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL11",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL6",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL6",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL7A",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL7A",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL4",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL4",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL23A",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL23A",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL10A",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL10A",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL18A",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL18A",big$names)])*100)),
  
  col=c("pink","blue","yellow"),ylab="proportion of K residues", ylim=c(0,40), names=c("stay at cytosol","back to nucleus","always at nucleus"),las=2 ,cex.axis=0.62)

pinktranscript<-
  c(c(str_count(as.character(small$seq)[2], "K")/str_count(as.character(small$seq)[2])*100),
    c(str_count(as.character(small$seq)[grep("RPS20",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS20",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS21",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS21",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS6",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS6",small$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPLP2",big$names)], "K")/str_count(as.character(big$seq)[grep("RPLP2",big$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS3",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS3",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS9",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS9",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS26",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS26",small$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL17",big$names)], "K")/str_count(as.character(big$seq)[grep("RPL17",big$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS17L",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS17L",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS28",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS28",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPSA",small$names)], "K")/str_count(as.character(small$seq)[grep("RPSA",small$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL10",big$names)], "K")/str_count(as.character(big$seq)[grep("RPL10",big$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS27A",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS27A",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS15A",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS15A",small$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL15",big$names)], "K")/str_count(as.character(big$seq)[grep("RPL15",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL37A",big$names)], "K")/str_count(as.character(big$seq)[grep("RPL37A",big$names)])*100))

bluetranscript<-
  c(c(str_count(as.character(small$seq)[grep("RPS7",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS7",small$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL38",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL38",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL19",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL19",big$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS4X",small$names)], "K")/ str_count(as.character(small$seq)[grep("RPS4X",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS16",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS16",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS14",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS14",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS11",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS11",small$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL21",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL21",big$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS10",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS10",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS23",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS23",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS19",small$names)], "K")/ str_count(as.character(small$seq)[grep("RPS19",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS12",small$names)], "K")/ str_count(as.character(small$seq)[grep("RPS12",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS5",small$names)], "K")/ str_count(as.character(small$seq)[grep("RPS5",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS3A",small$names)], "K")/ str_count(as.character(small$seq)[grep("RPS3A",small$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS25",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS25",small$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL23",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL23",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL26",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL26",big$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS8",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS8",small$names)])*100))

yellowtranscript<-
  c(c(str_count(as.character(big$seq)[grep("RPL28",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL28",big$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS27L",small$names)], "K")/str_count(as.character(small$seq)[grep("RPS27L",small$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL35A",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL35A",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL32",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL32",big$names)])*100),
    
    c(str_count(as.character(big$seq)[grep("RPL29",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL29",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL22",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL22",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL31",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL31",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL30",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL30",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL8",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL8",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL12",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL12",big$names)])*100),
    c(str_count(as.character(small$seq)[grep("RPS13",small$names)], "K")/ str_count(as.character(small$seq)[grep("RPS13",small$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL27",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL27",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL13A",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL13A",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL9",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL9",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL18",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL18",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL7",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL7",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL35AL",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL35AL",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL14",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL14",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL34",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL34",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL13",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL13",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL35",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL35",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL36",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL36",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL15",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL15",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL27A",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL27A",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL24",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL24",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL3",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL3",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL11",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL11",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL6",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL6",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL7A",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL7A",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL4",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL4",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL23A",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL23A",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL10A",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL10A",big$names)])*100),
    c(str_count(as.character(big$seq)[grep("RPL18A",big$names)], "K")/ str_count(as.character(big$seq)[grep("RPL18A",big$names)])*100))


t.test(yellowtranscript,bluetranscript)