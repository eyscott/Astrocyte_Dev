library(reshape2)
library(plyr)
library(dplyr)

#Get expression data and gene length
#TSO_1=PFA test, TSO_2=E16, TSO_3=E16, TSO_4=PO, TSO_5=PO
#TSO_6=Jenny, TSO_9=PO, TSO_10=PO
#only retrieving Ines samples (TSO_10,2,3,4,5,9--60 cells total)
exp <- read.table("Ines_counts.txt", header = T, stringsAsFactors = F)[ ,c(1,6,7:16,27:66,97:106)]

#readcounts
sum_counts<-data.frame(Sums=colSums(exp[,-c(1,2)]))
rownames(sum_counts) <- paste("1",rownames(sum_counts),sep="_")
sum_counts$id<-c(rep('new',10), rep('old',40),rep('new',10))

library(ggplot2)
library(scales)
library(wesanderson)
library(viridis)
# Basic violin plot
setwd('/Users/erica_1/Desktop/May_2023visit/Ines_PFA')
pdf(file='readcounts.pdf', width=6, height=4,bg="white")
ggplot(sum_counts, aes(x=id,y=Sums,fill=id)) + 
  geom_violin() + ylab("Read counts") + xlab(" ") +
  scale_fill_manual(values=(wes_palette(n=3, name="Royal1"))) +
  theme(legend.position="none") +geom_dotplot(binaxis='y', stackdir='center',
                                              position=position_dodge(1))
dev.off()

##most of the new ones (which are only P0) are better

##process data a bit
library(biomaRt)
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "mmusculus_gene_ensembl")

annot <- getBM(attributes = c('ensembl_gene_id','external_gene_name','gene_biotype'),
               values = exp$Geneid, 
               mart = ensembl)

exp_annot <- merge(exp,annot, by.x='Geneid', by.y="ensembl_gene_id")
exp_annot_unique <- distinct(exp_annot, external_gene_name, .keep_all= TRUE)#55252 
exp_unique_nopseudo <- exp_annot_unique[!grepl("pseudogene|misc", exp_annot_unique$gene_biotype),]#now 41031
rownames(exp_unique_nopseudo)<-exp_unique_nopseudo$external_gene_name

exp_annot_melt <- melt(exp_unique_nopseudo[ ,c(2:63)],id.vars= c("external_gene_name","Length"))

##normalize for gene length
exp_annot_melt$Length_kb <- (exp_annot_melt$Length)/1000
exp_annot_melt$value_L <- (exp_annot_melt$value/exp_annot_melt$Length_kb)
#get rid of data we don't need anymore
exp_annot_melt_red <- exp_annot_melt[ ,c("external_gene_name","variable","value_L")]

#normalize for library depth
exp_annot_melt_red_presum<- ddply(exp_annot_melt_red, c("variable"), summarise,
                                      Count_sum = sum(value_L))
exp_annot_melt_red_presum$RPM <- ((exp_annot_melt_red_presum$Count_sum)/1000000)
exp_annot_melt_red_RPM <- merge(exp_annot_melt_red, exp_annot_melt_red_presum, by="variable")
exp_annot_melt_red_RPM$value_L_RPM <- (exp_annot_melt_red_RPM$value_L/exp_annot_melt_red_RPM$RPM)
##get rid of some of the temp columns before making gene matrix
exp_annot_melt_red_RPM_red <- exp_annot_melt_red_RPM[ ,c("external_gene_name","variable", "value_L_RPM")]

exp_NORM <- dcast(exp_annot_melt_red_RPM_red, external_gene_name~variable,value.var='value_L_RPM',
                      fun.aggregate = mean, sep = "\t")
setwd('/Users/erica_1/Desktop/May_2023visit/Ines_PFA')
write.table(exp_NORM, "Ines_Norm_mat.txt")

#######

#convert NAs to zero
exp_NORM[is.na(exp_NORM)] <- 0

genes_per_cell <- data.frame(Matrix::colSums(exp_NORM>0))
ggplot(genes_per_cell, aes(x=rownames(genes_per_cell),y=genes_per_cell[ ,1])) + 
  geom_violin() + 
  ylab("Gene counts") + xlab(" ") +
  theme(legend.position="none") + geom_dotplot(binaxis='y', stackdir='center',
                                               position=position_dodge(1),dotsize = 0.6, color='white')

#look at gene content
exp_NORM_prep <- merge(exp_annot_unique,exp_NORM, by='external_gene_name')
exp_NORM_prepP <- exp_NORM_prep[!grepl("^IG|*_gene$|TEC|*_pseudogene$|ribozyme|pseudogene", exp_NORM_prep$gene_biotype),]#now 35747
exp_NORM_prepP_red<-exp_NORM_prepP[ ,c(1,64:124)]
exp_NORM_prepP_melt<-melt(exp_NORM_prepP_red, id.vars = c('external_gene_name','gene_biotype'))
colnames(exp_NORM_prepP_melt)<-c("Gene_name","gene_biotype","variable","value")

library(RColorBrewer)
n <-20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ggplot(exp_NORM_prepP_melt,aes(x=factor(1),weight=value,fill=gene_biotype)) + 
  geom_bar(width=1,alpha = 1) + xlab(" ") + 
  ylab(" ") + coord_polar("y") + 
  guides(fill=guide_legend(title="RNA type")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 10)) +
  theme(axis.text = element_text(colour="black", size = 8)) +
  scale_fill_manual(values = col_vector) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())

##remove cells with less than 100 genes detected
keep <- rownames(subset(genes_per_cell, Matrix..colSums.exp_NORM...0. >120))
exp_NORM_keep<-exp_NORM[,keep]
numbers<-c(11,seq(1,9,1))
colnames(exp_NORM_keep) <- c("GeneName",c(paste("P0_1", numbers, sep="."),paste("E16_1", c(11,2,3,4,5,6,8), sep="."), paste("E16_2", c(11,2,3,4,5,6,7,8), sep="."),
                       paste("P0_2", c(11,2,3,4,5,6,8), sep="."),paste("P0_3", c(11,2,3,4,5), sep="."),paste("P0_4", c(11,2,4:9), sep=".")))

##45 cells total: 15 E16, 30 P0


library(edgeR)
exp_unique_nopseudo
exp_raw_keep<-exp_unique_nopseudo[,keep]
rownames(exp_raw_keep)<-exp_raw_keep$external_gene_name
edgeR_input<-exp_raw_keep[ ,-1]
rep<-c(rep(1,10),rep(2,15),rep(1,20))
edgeR_input_norm = calcNormFactors(edgeR_input, method='TMM')
edgeR_input_disp<- estimateCommonDisp(edgeR_input, group=rep,
                                       lib.size=edgeR_input_norm, tol=1e-06,
                                       rowsum.filter=5, verbose=FALSE)

d <- DGEList(counts=edgeR_input, group=rep, lib.size=edgeR_input_norm)
##also to get gene names associated with the Gene IDs
dim(d) #[1] 41029    50

#testing using GLM
design <- model.matrix(~factor(rep))
y <- estimateDisp(d,design)
fit <- glmQLFit(y,design)
qlf_1v2 <- glmQLFTest(fit,contrast=c(-1,1))
#coef=ncol(fit$design) AKA number of factors here
qlf_topTags <- topTags(qlf_1v2, n=41031)
#get individual cpm normalized values for each sample to make a heatmap
qlf_topTags_ind <-rownames(topTags(qlf_1v2, n=41031)$table)
tagsbysample_qlf<-as.data.frame(cpm(d)[qlf_topTags_ind, order(d$samples$group)])
tagsbysample_qlf$Gene_id <-rownames(tagsbysample_qlf)
#retrieve genes with a significant pvalue
qlf_topTags_id<-subset(qlf_topTags$table, FDR < 0.05) #77
qlf_topTags_id$Gene_id <-rownames(qlf_topTags_id)
Genes_1v2<-qlf_topTags_id$Gene_id

##subset 1v2 genes
tagsbysample_qlf_sub<- subset(tagsbysample_qlf, Gene_id %in% Genes_1v2)
pheat_2<- as.data.frame(t(tagsbysample_qlf_sub[ ,c(1:45)]))
pheat_3<-tibble::rownames_to_column(pheat_2, "Gene_id")
#pheat_3$Gene_id<-gsub("X","",as.character(pheat_3$Gene_id))
rownames(pheat_3)<-pheat_3$Gene_id
pheat_4<-pheat_3[ ,-1]

##make metadata sheet for pheatmap
meta <- data.frame(cells = c(rownames(pheat_4)),
                  time = c(rep(2,30),rep(1,15)))
meta_ordered <- as.data.frame(meta[order(meta$time),])
Row_annot<-as.data.frame(meta_ordered[ ,-1])
rownames(Row_annot)<-meta_ordered$cells
colnames(Row_annot)<-'time'
time<-c('1'="#5A283E",'2'="#E8BD8C")
my_colour = list(time = time)

reorder_idx <- match(rownames(Row_annot),rownames(pheat_4)) 
pheat_5<-(pheat_4[reorder_idx])
pheat_6<-scale(pheat_4[reorder_idx])
library(pheatmap)
breaks<- c(0,seq(100,1000,100),seq(1001,5000,1000))
#breaks<- c(seq(-1,1,0.2))
pheatmap(pheat_5, annotation_row = Row_annot,
         cluster_rows = F,
         color = cividis(15),
         breaks = breaks,
         fontsize = 8,
         annotation_colors = my_colour)


##add upper and lower layer data to the comparison
layer<-c(rep('L',7),rep('U',8),rep('L',3),rep('U',2),rep('L',5),rep('L',7),
         'U','L',rep('U',3),'U','CC','U',rep('L',3),'CC','U')
date<-c(rep('E16',15),rep('P0',30))
meta_ordered_2<-data.frame(meta_ordered,
                           Birth= date,
                           Layer=layer)

Row_annot_2<-as.data.frame(meta_ordered_2[ ,-c(1,2)])
rownames(Row_annot_2)<-meta_ordered_2$cells
Birth<-c('E16'="#5A283E",'P0'="#E8BD8C")
Layer<-c('L'="#2E324C",'U'="#F2300F",'CC'="#649373")
my_colour_2 = list(Birth = Birth,
                 Layer=Layer)

pheatmap(pheat_5, annotation_row = Row_annot_2,
         cluster_rows = F,
         color = cividis(15),
         breaks = breaks,
         fontsize = 8,
         annotation_colors = my_colour_2,
         show_rownames = F,
         gaps_row = 30)
df_test<-as.data.frame(pheat_5)
row.names(df_test)<-paste(rownames(Row_annot_2),Row_annot_2$Birth,Row_annot_2$Layer, sep="_")
write.table(df_test, "Input_heatmap.txt")

###look at canonical markers and housekeeping genes
exp_NORM_keep
#prepare to subset
exp_NORM_keep_melt<-melt(exp_NORM_keep,id.vars = 'GeneName')
exp_NORM_keep_melt_sum<- ddply(exp_NORM_keep_melt, c("GeneName"), summarise,
                                  sum = sum(value),sd=sd(value))

exp_Norm_summarized <- merge(exp_NORM_keep,exp_NORM_keep_melt_sum, by="GeneName")
#Remove Gm genes
exp_Norm_summarized_noGm <- exp_Norm_summarized[!grepl("^Gm", exp_Norm_summarized$GeneName),]#now 25603

sub_HK <- subset(exp_Norm_summarized_noGm, c(sum > 500 & sd < 50)) #33
##still disparity between E16 and P0 samples
HK_G<-sub_HK$GeneName
CanG<-c('Sox9','Gfap','Alhd1l1','Rbfox3','Apoe','S100b','Lum')
HK_C_G<-c(HK_G, CanG)

exp_NORM_keep_submat<-exp_NORM_keep[exp_NORM_keep$GeneName %in% HK_C_G, ] 
cols<-c("GeneName","E16_1.11", "E16_1.2", "E16_1.3", "E16_1.4", "E16_1.5", 
        "E16_1.6","E16_1.8","E16_2.11", "E16_2.2","E16_2.3","E16_2.4","E16_2.5","E16_2.6", 
       "E16_2.7","E16_2.8","P0_1.11","P0_1.1","P0_1.2","P0_1.3","P0_1.4","P0_1.5","P0_1.6",  
       "P0_1.7","P0_1.8","P0_1.9","P0_2.11","P0_2.2","P0_2.3","P0_2.4","P0_2.5","P0_2.6",  
      "P0_2.8","P0_3.11","P0_3.2","P0_3.3","P0_3.4","P0_3.5","P0_4.11","P0_4.2", 
        "P0_4.4","P0_4.5","P0_4.6","P0_4.7","P0_4.8","P0_4.9")
exp_NORM_keep_submat_ordered<- exp_NORM_keep_submat[match(HK_C_G,exp_NORM_keep_submat$GeneName), cols]
exp_NORM_keep_submat_ordered_noNA<-na.omit(exp_NORM_keep_submat_ordered)
rownames(exp_NORM_keep_submat_ordered_noNA)<-exp_NORM_keep_submat_ordered_noNA$GeneName
exp_norm_sub_mat<-as.matrix(exp_NORM_keep_submat_ordered_noNA[ ,-1])

lmat = rbind(c(0,4),c(0,3),c(2,1))
lwid = c(0.5,4)
lhei = c(1,0.5,4)

library(gplots)
library(viridis)

col_breaks <- c(seq(0,10,0.5),seq(11,50,1),seq(51,100,5),seq(101,500,25))
setwd('/Users/erica_1/Desktop/May_2023visit/Ines_PFA')
pdf(file='Po_E16_HK.pdf', width=6, height=12,bg="white")
par(mar=c(9,6,6,3)+0.2, pin=c(0,0)) 
#bottom, left, top, and right
heatmap.2(exp_norm_sub_mat,    # data matrix
          trace="none", 
          labRow = rownames(exp_norm_sub_mat),
          margins =c(5,5),     # widens margins around plot
          col=viridis,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          #labCol = c("Penumbra 1","Penumbra 2", "Penumbra 3", "Penumbra 4",
          #           "Away 1","Away 2","Away 3","Away 4"),
          dendrogram="none",     # only draw a row dendrogram
          Colv=F,
          Rowv=F,
          cexCol = 1.1,
          adjCol = c(0.9,0.6),    #(left-right, top-bottom) 
          adjRow = c(0.35,0.2),
          keysize = 5,
          key.xlab = "TPM",
          key.title = NA,
          key.par = list(mar=c(1,3,1,6)),  #c(bottom, left, top, right)
          #densadj = 0.2,
          density.info="density",
          denscol = "white",
          lmat = lmat,
          lwid = lwid,
          lhei = lhei)                # turn off column clustering
dev.off()  # close the PDF device


##stacked violin plots 
GOI<-c('Gata3','Rplp2','Tecr','Trpm6','Alox8','Nsfl1c', 'Gm43660','Calm3','Lsm2','Rapgef3','Comt','Chpt1')
All_noCC<-subset(meta_ordered_2, Birth %in% c('E16','P0') & Layer %in% c('U',"L")) #43 cells
All_noCC_order<-All_noCC[order(All_noCC$Birth,All_noCC$Layer),]
keep_all<-c(All_noCC_order$cells,'external_gene_name')
test<-exp_NORM[ ,keep_all, drop = FALSE]
exp_NORM_keep_all_melt<-melt(test,id.vars = "external_gene_name")
exp_NORM_keep_all_melt_sub<- subset(exp_NORM_keep_all_melt, external_gene_name %in% GOI)
Vplot_input<-merge(exp_NORM_keep_all_melt_sub,All_noCC_order,by.x='variable',by.y='cells')
Vplot_input$BirthLayer<-paste(Vplot_input$Birth, "_", Vplot_input$Layer)
Vplot_input2<-Vplot_input[ ,-c(4:6)]
colnames(Vplot_input2)<-c("Cell","Feat","Expr","Idents")
Vplot_input2$Feat_f<-factor(Vplot_input2$Feat, levels=GOI)
library(cowplot)
ggplot(Vplot_input2, aes(factor(Idents), Expr, fill = Feat)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(Feat_f), scales = "free", switch = "y") +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0)) +
  xlab("Birth_layer") + ylab("Expression Level")+
  scale_fill_manual(values=sample(wes_palette("BottleRocket2", 12, type = "continuous")))
write.table(Vplot_input2, "Vplot_input.txt")
