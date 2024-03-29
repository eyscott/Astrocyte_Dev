data<-read.table('RawDataDaniela_2.txt', header = T, stringsAsFactors = F)
library(wesanderson)
library(ggplot2)
library(ggdendro)
library(ape)
library(magrittr)
library(dendextend)
library(gridExtra)
library(stringr)
#library(data.table)
#devtools::unload("data.table")

#factor all variable I will later probe
data$birth<-paste(str_split_fixed(data$TimEpoint,"_",2)[,1])
data$Layer_2<-paste("L",data$Layer,sep="_")
data$TimEpoint <- factor(x = data$TimEpoint, levels = c("E16_P21","E16_P56","P0_P21","P0_P56"))
data$Layer_2 <- factor(x = data$Layer_2, levels = c("L_1","L_2_3","L_4","L_5","L_6"))
data$Layer <- factor(x = data$Layer, levels = c("1","2_3","4","5","6"))
data$birth <- factor(x = data$birth, levels = c("E16","P0"))
data$merger<-paste(data$TimEpoint,data$Layer, sep = '_')
data$merger <- factor(x = data$merger, levels = c("E16_P21_1","E16_P21_2_3","E16_P21_4","E16_P21_5","E16_P21_6","E16_P56_1","E16_P56_2_3",
                                                 "E16_P56_4","E16_P56_5","E16_P56_6","P0_P21_1","P0_P21_2_3","P0_P21_4","P0_P21_5","P0_P21_6","P0_P56_1","P0_P56_2_3",
                                                 "P0_P56_4","P0_P56_5","P0_P56_6"))

#gather raw data versus z-caled data
data.data <- as.data.frame(data[ ,c(6:14)])
data.data.scale <- scale(as.data.frame(data[ ,c(6:14)]))

##convert scaled data into matrix with labeled rows
data.matrix<-as.matrix(data.data.scale)
rownames(data.matrix)<-data$TimEpoint

##start dendrogram construction
# Compute distances and hierarchical clustering
dd <- dist((data.matrix), method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
# Convert hclust into a dendrogram and plot
hcd <- as.dendrogram(hc)
##briefly plot
Time <- factor(x = data$TimEpoint, levels = c("E16_P21","E16_P56","P0_P21","P0_P56"))
Layer<-factor(x = data$Layer, levels = c("1","2_3","4","5","6"))
mycolT <- c("#EABE94", "#0B775E", "#35274A" ,"#F2300F")[Time]
mycolL <- c("#FAD510", "#CB2314", "#BDC881", "#354823", "#1E1E1E")[Layer]


allcols<-cbind(mycolT,mycolL)
#as phylo tree
plot(as.phylo(hc), type = "unrooted", cex = 0.6,
     no.margin = TRUE,tip.color = mycolT,label.offset=0)
#as circular plot
plot(as.phylo(hc), type = "fan",cex = 0.6,
     no.margin = TRUE,tip.color = mycolT,label.offset=0.2)

# Build dendrogram object from hclust results
dend <- data.matrix %>% # data
  dist %>% # calculate a distance matrix, 
  hclust(method = "ward.D2") %>% # Hierarchical clustering 
  as.dendrogram # Turn the object into a dendrogram.

dend %>% set("labels_col", c("#EABE94", "#0B775E", "#35274A" ,"#F2300F")) %>% # change color
  set("labels_cex", 1) %>% # Change size
  set("branches_k_color", k = 7) %>%
  plot(main = " ") 

##colour bar
# Set the plot margin: bottom, left, top & right
par(mar = c(10, 3, 3, 4) + 0.1,
    xpd = NA) # allow content to go into outer margin 
# Setup the color bar based on E16 & P0
plot(dend)
the_bars_birth <- ifelse(grepl("E16", dend_labels),"#0B775E", "#35274A")
the_bars_end <- ifelse(grepl("P56", dend_labels),"firebrick3", "beige")
the_bars <- cbind(the_bars_birth, the_bars_end)
#colored_bars(colors = the_bars, dend = dend, rowLabels = c("birth", "end"))
#colored_bars(colors = mycolL, dend = dend,rowLabels = "Layer")
colored_bars(colors = allcols, dend = dend,rowLabels = c("Time","Layer"))

#more refined example
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")
# Set the plot margin: bottom, left, top & right
par(mar = c(5, 3, 3, 4) + 0.1,
    xpd = NA) # allow content to go into outer margin 
# Setup the color bar based on E16 & P0
dend %>% set("branches_k_color", k = 7,value =sample(wes_palette(n=7, name="Zissou1Continuous"))) %>%
  plot(main = " ",leaflab = "none") # plot horizontle with "horiz = TRUE, "
colored_bars(colors = allcols, dend = dend,rowLabels = c("Time","Layer"))

dend <- data.matrix %>% dist %>% 
  hclust(method = "ward.D2") %>% as.dendrogram  %>%
  set("branches_k_color", k = 7,value =sample(wes_palette(n=7, name="Zissou1Continuous"))) %>%
  set("labels_cex", 0.2)

df2<-colored_bars(colors = allcols, dend = dend,rowLabels = c("Time","Layer"))
labels<- dend %>% labels
df2<- cbind(df2,'ID'= as.data.frame(dend %>% labels))
colnames(df2)[3]<-"ID"

ggd1 <- as.ggdend(dend)
ggplot(ggd1, theme = theme_minimal()) + ylab('Height') + xlab(" ")
library(data.table)
ggd1_df <- as.data.frame(rbindlist(ggd1, fill = TRUE))
write.table(ggd1_df,"Dend_input.txt")

p1<-ggplot(ggd1, theme = theme_minimal()) + ylab('Height') + xlab(" ")
p2<-ggplot(df2,aes(ID,y=1,fill=factor(mycolT)))+geom_tile()+
  scale_y_continuous(expand=c(0,0))+
  theme(axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        legend.position="none")
p1/p2

dend_prop<-as.data.frame(t(dendextend:::cutree.dendrogram(dend,h=15)), stringsAsFactors=T)
t_dend_prop = setNames(data.frame(t(dend_prop[,-1])), dend_prop[,1])
t_dend_prop <- tibble::rownames_to_column(t_dend_prop, "ID")
colnames(t_dend_prop)<-c("ID","cluster")
#remove .X after each label
t_dend_prop_2<-(cbind(t_dend_prop,(t(data.frame(strsplit(as.character(t_dend_prop$ID),".",fixed = T))))[ ,1]))
colnames(t_dend_prop_2)<-c("ID","cluster","trmt")
#tab <- table(t_dend_prop_2$trmt, t_dend_prop_2$cluster)
tab <- table(t_dend_prop_2$cluster,t_dend_prop_2$trmt)
tab

ptab <- as.data.frame(prop.table(tab, margin=1))

#colx <- c("#91BAB6","#6FB2C1","#E3B710","#3A9AB2","#A5C2A3","#BDC881",)
colx <-c("#EABE94", "#0B775E", "#35274A" ,"#F2300F")
ptab$Var2<-factor(x = ptab$Var2, levels = c("E16_P21","E16_P56","P0_P21","P0_P56"))
Pt <- data.table(ptab)
Pt[,midpoint:=cumsum(Freq) - Freq/2,by=Var1]
ptab_col<-Pt[ ,colx:=colx,by=Var1]



vysg <- ggplot(ptab_col, aes(x=1,y=Freq,fill=Var2)) + 
  geom_bar(stat="identity",width=2) + 
  coord_polar(theta='y')+
  theme(axis.ticks=element_blank(), axis.title=element_blank(),axis.text.y = element_blank(),axis.text.x = element_blank(), panel.grid  = element_blank())+
  #geom_text(aes(x = 2.5, y = midpoint, label = Var1, colour = I(colx)))+
  #scale_x_continuous(limits=c(-1,2.5))+
  scale_fill_manual('Timepoints',values=colx)
vysg<-vysg+facet_wrap(~ Var1, ncol = 7)
vysg
write.table(ptab_col,"Time_ident_Pie.txt")
##now for pies for layers
rownames(data.matrix)<-data$Layer_2

dend <- data.matrix %>% # data
  dist %>% # calculate a distance matrix, 
  hclust(method = "ward.D2") %>% # Hierarchical clustering 
  as.dendrogram # Turn the object into a dendrogram.

dend %>% set("labels_col", c("#9A8822", "#F5CDB4", "#F8AFA8", "#FDDDA0", "#74A089")) %>% # change color
  set("labels_cex", 1) %>% # Change size
  set("branches_k_color", k = 7) %>%
  plot(main = " ") 

dend_prop<-as.data.frame(t(dendextend:::cutree.dendrogram(dend,h=15)), stringsAsFactors=T)
t_dend_prop = setNames(data.frame(t(dend_prop[,-1])), dend_prop[,1])
t_dend_prop <- tibble::rownames_to_column(t_dend_prop, "ID")
colnames(t_dend_prop)<-c("ID","cluster")
#remove .X after each label
t_dend_prop_2<-(cbind(t_dend_prop,(t(data.frame(strsplit(as.character(t_dend_prop$ID),".",fixed = T))))[ ,1]))
colnames(t_dend_prop_2)<-c("ID","cluster","Layer")
#tab <- table(t_dend_prop_2$trmt, t_dend_prop_2$cluster)
library(data.table)
tab <- table(t_dend_prop_2$cluster,t_dend_prop_2$Layer)
tab

ptab <- as.data.frame(prop.table(tab, margin=1))

colx <- c("#91BAB6","#F8AFA8","#E3B710","#3A9AB2","#BDC881")
ptab$Var2<-factor(x = ptab$Var2, levels = c("L_1","L_2_3","L_4","L_5","L_6"))
Pt <- data.table(ptab)
Pt[,midpoint:=cumsum(Freq) - Freq/2,by=Var1]
ptab_col<-Pt[ ,colx:=colx,by=Var1]

vysg <- ggplot(ptab_col, aes(x=1,y=Freq,fill=Var2)) + 
  geom_bar(stat="identity",width=2) + 
  coord_polar(theta='y')+
  theme(axis.ticks=element_blank(), axis.title=element_blank(),axis.text.y = element_blank(),axis.text.x = element_blank(), panel.grid  = element_blank())+
  #geom_text(aes(x = 2.5, y = midpoint, label = Var1, colour = I(colx)))+
  #scale_x_continuous(limits=c(-1,2.5))+
  scale_fill_manual('Layer',values=colx)
vysg<-vysg+facet_wrap(~ Var1, ncol = 7)
vysg
write.table(ptab_col, "Layer_input_Pie.txt")
##now pies for merger
rownames(data.matrix)<-data$merger

dend <- data.matrix %>% # data
  dist %>% # calculate a distance matrix, 
  hclust(method = "ward.D2") %>% # Hierarchical clustering 
  as.dendrogram # Turn the object into a dendrogram.

dend %>% set("labels_col", wes_palette("Zissou1", 21, type = c("continuous"))) %>% # change color
  set("labels_cex", 1) %>% # Change size
  set("branches_k_color", k = 7) %>%
  plot(main = " ") 


dend_prop<-as.data.frame(t(dendextend:::cutree.dendrogram(dend,h=15)), stringsAsFactors=T)
t_dend_prop = setNames(data.frame(t(dend_prop[,-1])), dend_prop[,1])
t_dend_prop <- tibble::rownames_to_column(t_dend_prop, "ID")
colnames(t_dend_prop)<-c("ID","cluster")
#remove .X after each label
t_dend_prop_2<-(cbind(t_dend_prop,(t(data.frame(strsplit(as.character(t_dend_prop$ID),".",fixed = T))))[ ,1]))
colnames(t_dend_prop_2)<-c("ID","cluster","merger")
#tab <- table(t_dend_prop_2$trmt, t_dend_prop_2$cluster)
library(data.table)
tab <- table(t_dend_prop_2$cluster,t_dend_prop_2$merger)
tab

ptab <- as.data.frame(prop.table(tab, margin=1))

colx <- wes_palette("Zissou1", 20, type = c("continuous"))
ptab$Var2<-factor(x = ptab$Var2, levels = c("E16_P21_1","E16_P21_2_3","E16_P21_4","E16_P21_5","E16_P21_6","E16_P56_1","E16_P56_2_3",
                                            "E16_P56_4","E16_P56_5","E16_P56_6","P0_P21_1","P0_P21_2_3","P0_P21_4","P0_P21_5","P0_P21_6","P0_P56_1","P0_P56_2_3",
                                            "P0_P56_4","P0_P56_5","P0_P56_6"))
Pt <- data.table(ptab)
Pt[,midpoint:=cumsum(Freq) - Freq/2,by=Var1]
ptab_col<-Pt[ ,colx:=colx,by=Var1]

vysg <- ggplot(ptab_col, aes(x=1,y=Freq,fill=Var2)) + 
  geom_bar(stat="identity",width=2) + 
  coord_polar(theta='y')+
  theme(axis.ticks=element_blank(), axis.title=element_blank(),axis.text.y = element_blank(),axis.text.x = element_blank(), panel.grid  = element_blank())+
  #geom_text(aes(x = 2.5, y = midpoint, label = Var1, colour = I(colx)))+
  #scale_x_continuous(limits=c(-1,2.5))+
  scale_fill_manual('Merged',values=colx)
vysg<-vysg+facet_wrap(~ Var1, ncol = 7)
vysg
write.table(ptab_col, "All_ident.txt")
