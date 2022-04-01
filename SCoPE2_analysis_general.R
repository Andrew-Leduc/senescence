source("/Users/andrewleduc/Desktop/GitHub/senescence/functions_parameters.R")

# User specific

# Reference channel number (1-11, or 1-16)
ref_channel<-2

# Add your cell type labels, must match those used in experimental design --> these are specified in annotation file
your_labels<-c("sc","neg")
your_control_label<-"neg"

# Import ------------------------------------------------------------------

# Marker genes

genes <- read.csv('genes.csv',header = F)


# For converting between uniprot and gene name

# convert_mouse <- read.fasta('uniprotMusmusculusSP.fasta',set.attributes = T,whole.header = T)
# convert_mouse <- names(convert_mouse)
# parse_row<-grep("GN=",convert_mouse, fixed=T)
# split_prot<-str_split(convert_mouse[parse_row], pattern = fixed("GN="))
# gene<-unlist(split_prot)[seq(2,2*length(split_prot),2)]
# prot <- unlist(split_prot)[seq(1,2*length(split_prot),2)]
# prot_parse <- grep("|",prot, fixed=T)
# gene_parse <- grep(" ",gene, fixed=T)
# split_gene<-str_split(gene[parse_row], pattern = fixed(" "))
# split_gene<-unlist(split_gene)[seq(1,3*length(split_gene),3)]
# split_prot<-str_split(prot[parse_row], pattern = fixed("|"))
# split_prot<-unlist(split_prot)[seq(2,3*length(split_prot),3)]
# convert_mouse  <- as.data.frame(cbind(split_prot,split_gene))
# convert_mouse$split_gene <- toupper(convert_mouse$split_gene)


# Load raw data 

ev<-read.delim("/Volumes/GoogleDrive/My Drive/MS/Collaborators/Senescence/DART/ev_updated.txt", header = TRUE)

#Here just cleaning up protein names from the fasta to uniprot accession (may not be neccisary)
# Parse sp|Q00000|HUMAN_XXX into just Uniprot accession: Q00000  
  
  parse_row<-grep("|",ev$Leading.razor.protein, fixed=T)
  
  if(length(parse_row)>0){
    split_prot<-str_split(ev$Leading.razor.protein[parse_row], pattern = fixed("|"))
    split_prot2<-unlist(split_prot)[seq(2,3*length(split_prot),3)]
    ev$Leading.razor.protein[parse_row]<-split_prot2
  }
  parse_row<-grep("-",ev$Leading.razor.protein, fixed=T)
  if(length(parse_row)>0){
    split_prot<-str_split(ev$Leading.razor.protein[parse_row], pattern = fixed("-"))
    split_prot2<-unlist(split_prot)[seq(1,2*length(split_prot),2)]
    ev$Leading.razor.protein[parse_row]<-split_prot2
  }
  
  
# Load experimental design and batches
  
  design<-read.csv("/Users/andrewleduc/Desktop/GitHub/senescence/annotation_18plex.csv")

  batch<-read.csv("/Users/andrewleduc/Desktop/GitHub/senescence/batch.csv")

  # Attach batch data to protein data
  ev[,colnames(batch)[-1]]<-NA
  for(X in batch$set){

    ev$lcbatch[ev$Raw.file==X] <- as.character(batch$lcbatch[batch$set%in%X])
    ev$sortday[ev$Raw.file==X] <- as.character(batch$sortday[batch$set%in%X])
    ev$digest[ev$Raw.file==X] <- as.character(batch$digest[batch$set%in%X])

  }


# Create unique peptide+charge column:
ev$modseq<-paste0(ev$Modified.sequence,ev$Charge)


# Add X in front of experiment names because R doesn't like column names starting with numbers
ev$Raw.file<-paste0("X",ev$Raw.file)
design$Set<-paste0("X",design$Set)


# Which columns hold the TMT Reporter ion (RI) data
ri.index<-which(colnames(ev)%in%paste0("Reporter.intensity.",1:18))

# Make sure all runs are described in design, if not, print and remove them:
not.described<-unique(ev$Raw.file)[ !unique(ev$Raw.file) %in% paste0(design$Set) ]
ev<-ev[!ev$Raw.file%in%not.described,]
unique(ev$Raw.file)

# Filter out reverse hits, contaminants, and contaminated spectra...
ev<-ev[-which(ev$Reverse=="+"),]
if(length(grep("REV", ev$Leading.razor.protein))>0){ ev<-ev[-grep("REV", ev$Leading.razor.protein),] }
if(length(grep("CON", ev$Leading.razor.protein))>0){ ev<-ev[-grep("CON", ev$Leading.razor.protein),] }
if(length(which(ev$Potential.contaminant=="+"))>0){ ev<-ev[-which(ev$Potential.contaminant=="+"),] }
ev<-ev[!is.na(ev$PIF),]
ev<-ev[ev$PIF>0.8,]

# Remove peptides that are more the 10% the intensity of the carrier in the single cell runs (only)
ev<-as.data.frame(ev)
ev$mrri<-0
ev$mrri <- rowMeans(ev[, ri.index[6:length(ri.index)]] / ev[, ri.index[1]], na.rm = T)
ev<-ev[ev$mrri < 0.05, ]

#ev <- ev %>% filter(Intensity > median(ev$Intensity,na.rm = T))


# Filter by PEP or FDR: CHOOSE ONE

ev<-ev[ev$dart_qval<0.01, ]
#ev<-ev[ev$PEP<0.03, ]
# ev<-ev[calc_fdr(ev$dart_PEP)<0.01, ]
#ev<-ev[calc_fdr(ev$PEP)<0.01, ]





# Normalize single cell runs to normalization channel
ev<-as.data.frame(ev)
ev[, ri.index] <- ev[, ri.index] / ev[, ri.index[ref_channel]]




#Here is just some data wrangling so all appropriate annotations can be matched to single cells 


# Organize data into a more convenient data structure:
# Create empty data frame
ev.melt<-melt(ev[0, c("Raw.file","modseq","Leading.razor.protein","lcbatch","sortday","digest", colnames(ev)[ri.index]) ],
              id.vars = c("Raw.file","modseq","Leading.razor.protein","lcbatch","sortday","digest"))

colnames(ev.melt)<-c("Raw.file","sequence","protein","lcbatch","sortday","digest","celltype","quantitation")



# Record mapping of cell type to Channel:
ct.v<-c()
qt.v<-c()

# Create a unique ID string
unique.id.numeric<-1:length(ri.index)
unique.id<-paste0("i",unique.id.numeric)

RI_keep<-ri.index
ev <- ev %>% filter(is.na(Raw.file)==F)
# Give each sample a unique identifier
for(X in unique(ev$Raw.file)){

  # Subset data by X'th experiment
  ev.t<-ev[ev$Raw.file%in%X, ]

  # Name the RI columns by what sample type they are: carrier, single cell, unused, etc...
  colnames(ev.t)[ri.index]<-paste0(as.character(unlist(design[design$Set==X,-1])),"-", unique.id)

    # Melt it! and combine with other experimental sets
    ev.t.melt<-melt(ev.t[, c("Raw.file","modseq","Leading.razor.protein","lcbatch","sortday","digest", colnames(ev.t)[RI_keep]) ],
                    id.vars = c("Raw.file","modseq","Leading.razor.protein","lcbatch","sortday","digest"));

    # Record mapping of cell type to Channel:
    ct.v<-c(ct.v, unique.id[which(ri.index%in%RI_keep)] )
    qt.v<-c(qt.v, colnames(ev)[RI_keep] )

    colnames(ev.t.melt)<-c("Raw.file","sequence","protein","lcbatch","sortday","digest","celltype","quantitation")

    ev.melt<-rbind(ev.melt, ev.t.melt)

  # Update unique ID string
  unique.id.numeric<-unique.id.numeric + length(ri.index)
  unique.id<-paste0("i", unique.id.numeric)

}

c2q<-data.frame(ct.v, qt.v); colnames(c2q)<-c("id","channel")
c2q$id <- as.character(c2q$id )

# Grab the unique number associate to each and every cell, carrier channel, and empty channel
ev.melt$id<-unlist(strsplit(as.character(ev.melt$celltype),"-"))[seq(2,length(unlist(strsplit(as.character(ev.melt$celltype),"-"))),2)]
ev.melt$celltype<-unlist(strsplit(as.character(ev.melt$celltype),"-"))[seq(1,length(unlist(strsplit(as.character(ev.melt$celltype),"-"))),2)]
ev.melt$id<-as.factor(ev.melt$id)

# Remove duplicate observations of peptides from a single experiment
ev.melt<-remove.duplicates(ev.melt,c("sequence","id") )
ev.melt<-ev.melt[!is.na(ev.melt$protein), ]

# Create additional meta data matrices
ev.melt.uniqueID<-remove.duplicates(ev.melt,"id")
ev.melt.uniqueID <- left_join(ev.melt.uniqueID,c2q, by = 'id')
ev.melt.pep<-remove.duplicates(ev.melt, c("sequence","protein") )

# Create data frame of peptides x cells, populated by quantitation
ev.unmelt<-dcast(ev.melt, sequence ~ id, value.var = "quantitation", fill=NA)

# Also create matrix of same shape
ev.matrix<-as.matrix(ev.unmelt[,-1]); row.names(ev.matrix)<-ev.unmelt$sequence

# Replace all 0s with NA
ev.matrix[ev.matrix==0]<-NA
ev.matrix[ev.matrix==Inf]<-NA
ev.matrix[ev.matrix==-Inf]<-NA



# Divide matrix into single cells (including intentional blanks) and carriers
sc_cols<-unique(ev.melt$id[(ev.melt$celltype%in%c(your_labels))])
ev.matrix.sc<-ev.matrix[, colnames(ev.matrix)%in%sc_cols]




# This is Filter single cells,  Here we are commuting the agreement between peptides mapping to a protein
# For each protien with atleast 3 peptides we get a score, then the averge of that for all proteins is the 
# value for each cell


sc.melt<-ev.melt

xd<-as_tibble( sc.melt )

xd <- xd %>% dplyr::group_by(id) %>% dplyr::mutate(med_per_c = median(quantitation, na.rm=T)); length(unique(xd$id))

length(unique(xd$id))

xd$quantitation[(xd$quantitation)==Inf]<-NA
xd$quantitation[(xd$quantitation)==0]<-NA

xd <- xd %>% dplyr::mutate_if(is.factor, as.character)

xd1 <- xd %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(norm_q1 = quantitation / median(quantitation, na.rm=T))

xd2 <- xd1 %>%
  dplyr::group_by(sequence, Raw.file) %>%
  dplyr::mutate(norm_q = quantitation / median(norm_q1, na.rm=T))

xd3<- xd2 %>%
  dplyr::filter(celltype%in%c(your_labels))

xd4<- xd3 %>%
  dplyr::group_by(protein, id) %>%
  mutate(cvq = cv(norm_q))

xd5<- xd4 %>%
  dplyr::group_by(protein, id) %>%
  dplyr::mutate(cvn = cvna(norm_q))

xd6<- xd5 %>%
  filter(cvn > 2)


xd7<-xd6 %>% group_by(id) %>% dplyr::mutate(cvm=median(cvq, na.rm=T))

xdf<-xd7
look <- xdf %>% filter(control == 'ctl')
look <- look %>% distinct(Raw.file,.keep_all = T)

#print("Number of unique proteins used in calculation:", length(unique(xdf$protein)))

median(xdf$cvm)
# Filter out variable wells and controls
cvPar <- 0.47
sc_kept<-unique( xdf$id[xdf$celltype!=your_control_label & xdf$cvm < cvPar])
sc0_kept<-unique( xdf$id[xdf$celltype==your_control_label & xdf$cvm > cvPar])


# Which wells to keep
keep_these<-unique( xdf$id)

sc_total<-unique( xdf$id[xdf$celltype!=your_control_label])
sc0_total<-unique( xdf$id[xdf$celltype==your_control_label])
scrate<-round(length(sc_kept) / length(sc_total),2)*100

ev.matrix.sc.f<-ev.matrix.sc[,colnames(ev.matrix.sc)%in%sc_kept]; dim(ev.matrix.sc.f)
ev.matrix.sc.f[ev.matrix.sc.f==Inf]<-NA
ev.matrix.sc.f[ev.matrix.sc.f==-Inf]<-NA
ev.matrix.sc.f[ev.matrix.sc.f==0]<-NA

xdf$control<-"sc"
xdf$control[xdf$celltype==your_control_label]<-"ctl"

my_col3<-c( "black", "purple2")


ggplot(data=xdf, aes(x=cvm,fill=control)) + geom_density(aes( alpha=0.5), adjust=4) + theme_pubr() +
  scale_fill_manual(values=my_col3[c(1,2)]) + 
  xlab("Quantification variability") + ylab("Fraction of cells") + rremove("y.ticks") + rremove("y.text") +
  font("xylab", size=30) +
  font("x.text", size=20) +
  coord_cartesian(xlim=c(.0,1))+
  #xlim(c(-0.15, 0.35)) +
  # annotate("text", x=0.27, y= 14, label=paste0(scrate,"% single cells passed"), size=8, color=my_col3[c(2)])+
  # annotate("text", x=0.27, y= 12.5, label=paste0(sc0rate,"% control wells passed"), size=8, color=my_col3[c(1)])+
  annotate("text", x=0.26, y= 14, label=paste0(length(sc_kept),"cells"), size=10, color=my_col3[c(2)])+
  annotate("text", x=0.64, y= 12, label=paste0(length(sc0_kept)," Ctr -"), size=10, color=my_col3[c(1)])+
  annotate("text", x=0.63, y= 14, label=paste0(length(sc_total) -length(sc_kept)," cells"), size=10, color=my_col3[c(2)])+
  annotate("text", x=0.245, y= 12, label=paste0(length(sc0_total) - length(sc0_kept)," Ctr -"), size=10, color=my_col3[c(1)])+
  #annotate("text", x=0.25, y= 3, label="Macrophage-like", size=6) +
  rremove("legend") + geom_vline(xintercept=0.46, lty=2, size=2, color="gray50") + theme(plot.margin = margin(0, 1, 0, 0, "cm"))




# Data transformations ----------------------------------------------------

# matrix is now peptide(rows) x single cell(columns)
## So first the data is normalized by the reference, then we divide each column (a single cell) by median
#this normalizes differences from single cells being different sizes

#then we normalize each row by mean, this makes it so more abundant peptides dont have larger values
#than less abundant peptides because we are only interested in relative protein changes, because
#absolute protein changes are less significant because absolute quant is more difficult

#then we colapse peptide measurements to protein level by taking median of all peptides mapping to a protien
#then we impute data with knn (this step is optional if you dont want to impute)

# Then we batch correct data (we can talk about what type you may want to apply to data)

#Then we compute PCA





# Perform normalizations / transformations in multiple steps with visual sanity checks:
b.t<-"FD"
xlim.t<-c(-2,2)
par(mfrow=c(3,3))

# Original data, normalized to reference channel, filtered for failed wells:
t0<-ev.matrix.sc.f

hist(c(t0), breaks=b.t, xlim=xlim.t)

# Column then row normalize by median or mean (see source functions):
t1<-cr_norm(t0)
hist(c(t1), breaks=b.t, xlim=xlim.t)


# Filter for missing data:
t2<-filt.mat.rc(t1, na.row, na.col)
hist(c(t2), breaks=b.t, xlim=xlim.t)


# Log2 transform:
t3<-log2(t2)
t3[t3==Inf]<-NA
t3[t3==-Inf]<-NA
t3[t3==0]<-NA
hist(c(t3), breaks=b.t, xlim=xlim.t)


# # Collapse to protein level by median:
t3m<-data.frame(t3)
t3m$pep<-rownames(t3)
t3m$prot <- ev.melt.pep$protein[match(t3m$pep, ev.melt.pep$sequence)]
t3m<-melt(t3m, variable.names = c("pep", "prot"))
colnames(t3m) <-c("pep","prot","id","quantitation")
t3m2<- t3m %>% group_by(prot, id) %>% dplyr::summarize(qp = median(quantitation, na.rm=T))
t4m<-dcast(t3m2, prot ~ id, value.var = "qp", fill=NA)
t4<-as.matrix(t4m[,-1]); row.names(t4)<-t4m[,1]
hist(c(t4), breaks=b.t, xlim=xlim.t)

# Re-column and row normalize:
t4b<-cr_norm_log(t4)
hist(c(t4b), breaks=b.t, xlim=xlim.t)

# Assign to a final variable name:
ev.matrix.sc.f.n<-t4b

dim(t4b)

## Impute single celldata
imp.input<-ev.matrix.sc.f.n
sc.imp <- hknn(imp.input, k.t)
t5<-sc.imp
sum(is.na(sc.imp))
dim(sc.imp)

sc.imp[(is.na(sc.imp))]<-0



# Define the batches and model:
batch.covs <- ev.melt.uniqueID$Raw.file[match(colnames(sc.imp), ev.melt.uniqueID$id)]
mod<-data.frame(ev.melt.uniqueID$celltype[match(colnames(sc.imp), ev.melt.uniqueID$id)]); colnames(mod)<-"celltype"
#mod<-model.matrix(~as.factor(celltype), data=mod)

batch.covsRI <- ev.melt.uniqueID$channel[match(colnames(sc.imp), ev.melt.uniqueID$id)]

#Cell Type Info Labels
#mod<-data.frame(ev.melt.uniqueID$celltype[match(colnames(sc.imp), ev.melt.uniqueID$id)]); colnames(mod)<-"celltype"
mod<-model.matrix(~as.factor(celltype), data=mod)
matrix.sc.batch <- ComBat(sc.imp, batch=batch.covs)

#matrix.sc.batch <- removeBatchEffect(sc.imp, batch=factor(batch.covs), batch2=factor(batch.covsRI))#,design=mod)

t6<-matrix.sc.batch

# visual sanity checks post-imputation:
hist(c(t5), breaks=b.t, xlim=xlim.t)
hist(c(t6), breaks=b.t, xlim=xlim.t)

par(mfrow=c(1,1))

t7 <- cr_norm_log(t6)
t7[is.na(t4b)==T] <- NA




# PCA ------------------------------------------------------------------------

# PC's to display:
PCx<-"PC1"
PCy<-"PC2"

# Data to use:


mat.sc.imp<-cr_norm_log(t6)
# Dot product of each protein correlation vector with itself
r1<-cor(t(matrix.sc.batch))
rsum<-rowSums(r1^2)

# Calculate the weighted data matrix:

X.m <- mat.sc.imp
X.m <- diag(rsum) %*%  X.m
pca.imp.cor <- cor(X.m,use = 'pairwise.complete.obs')

# PCA
sc.pca<-eigen(pca.imp.cor)
scx<-as.data.frame(sc.pca$vectors)
colnames(scx)<-paste0("PC",1:ncol(scx))
scx$cells<-colnames(pca.imp.cor)

# Percent of variance explained by each principle component
pca_var <- sc.pca$values
percent_var<- pca_var/sum(pca_var)*100
plot(1:length(percent_var), percent_var, xlab="PC", ylab="% of variance explained")

# Map meta data
pca.melt <- melt(scx); colnames(pca.melt)<-c("id","pc","value")
# Re map ...
pca.display <- dcast(pca.melt, id ~ pc, value.var = "value", fill=NA)
pca.display <- pca.display %>% left_join(ev.melt.uniqueID,by = 'id')


ggscatter(pca.display,x =PCx, y = PCy ,  size = 4, alpha=0.5) +
  xlab(paste0(PCx,"  (", round(percent_var[1],0),"%)")) +
  ylab(paste0(PCy,"  (", round(percent_var[2],0),"%)")) +
  font("ylab",size=30) +
  font("xlab",size=30) +
  font("xy.text", size=20) +
  #rremove("legend") +
  #scale_color_manual(values = my_colors[2:3]) +
  #annotate("text", x=-0.025, y=0.21,label=your_control_label, color=my_colors[2], size=10)  +
  #annotate("text", x=0.03, y=0.21, label="Monocyte", color=my_colors[3], size=10) +
  #annotate("text", x=-.03, y=0.12, label=paste0(dim(mat.sc.imp)[1], " proteins"), size=8) +
  #annotate("text", x=-0.03, y=0.07, label=paste0(dim(mat.sc.imp)[2], " cells"), size=8) +
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
                       high = "red") #+ theme(plot.title = element_text(size = 20))


c1 <- pca.display %>% filter(PC1 > .02)
c1 <- c1 %>% filter(PC2 < 0)
c2 <- pca.display %>% filter(PC1 > .02)
c2 <- c2 %>% filter(PC2 > 0)
c3 <- pca.display %>% filter(PC1 < .0) # this is the AT1 AT2 basel population


### Here we are getting things in an accessible format for doing differential espression or Gene set enrichment

data <- melt(mat.sc.imp)
data$Condition <- NA
data$Condition[data$Var2 %in% c1$id] <- 'c1'
data$Condition[data$Var2 %in% c2$id] <- 'c2'
data$Condition[data$Var2 %in% c3$id] <- 'c3'
data <-data %>% filter(is.na(Condition)==F)


colnames(data)<- c('Protein','SampID','Intensity','Condition')

GSEA_out <- Population_GSEA(data,GO_db,3,1,converter_df)
GSEA_out <- GSEA_out %>% filter(is.na(pVal)==F)
GSEA_out$qval <- p.adjust(GSEA_out$pVal,method = 'BH')
GSEA_out <- GSEA_out %>% filter(qval < .01)
GSEA_out_ <-GSEA_out %>% filter(numberOfMatches > 3)



## Diff protein on on AT1, AT2, Basel

##### <- make PCA on subcluster, this involves re-normalizing data

# filtering for genes with 20 percent present values 
t7_filt <- filt.mat.rc(t7,.8,.99)
dim(t7_filt)

#select for C3 population (basel,AT1,AT2)
t7_filt <- as.data.frame(t7_filt) %>% dplyr::select(c3$id)

## Here is the renormalization
t7_filt <- cr_norm_log(as.matrix(t7_filt))

X.m <- as.matrix(t7_filt)
pca.imp.cor <- cor(X.m,use = 'pairwise.complete.obs')
sc.pca<-eigen(pca.imp.cor)
scx<-as.data.frame(sc.pca$vectors)
colnames(scx)<-paste0("PC",1:ncol(scx))
scx$cells<-colnames(pca.imp.cor)

# Percent of variance explained by each principle component
pca_var <- sc.pca$values
percent_var<- pca_var/sum(pca_var)*100
plot(1:length(percent_var), percent_var, xlab="PC", ylab="% of variance explained")

# Map meta data
pca.melt <- melt(scx); colnames(pca.melt)<-c("id","pc","value")
# Re map ...
pca.display <- dcast(pca.melt, id ~ pc, value.var = "value", fill=NA)
pca.display <- pca.display %>% left_join(ev.melt.uniqueID,by = 'id')

ggscatter(pca.display,x =PCx,y = PCy ,  size = 4, alpha=0.5) +
  xlab(paste0(PCx,"  (", round(percent_var[1],0),"%)")) +
  ylab(paste0(PCy,"  (", round(percent_var[2],0),"%)")) +
  font("ylab",size=30) +
  font("xlab",size=30) +
  font("xy.text", size=20) 
  #rremove("legend") +
  #scale_color_manual(values = my_colors[2:3]) +
  #annotate("text", x=-0.025, y=0.21,label=your_control_label, color=my_colors[2], size=10)  +
  #annotate("text", x=0.03, y=0.21, label="Monocyte", color=my_colors[3], size=10) +
  #annotate("text", x=-.03, y=0.12, label=paste0(dim(mat.sc.imp)[1], " proteins"), size=8) +
  #annotate("text", x=-0.03, y=0.07, label=paste0(dim(mat.sc.imp)[2], " cells"), size=8) +
  #scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
  #                      high = "red") #+ theme(plot.title = element_text(size = 20))


### Define new clusters manually

c1 <- pca.display %>% filter(PC1 < 0)
c1 <- c1 %>% filter(PC2 > -.025)
c2 <- pca.display %>% filter(PC1 > .0)
c3 <- pca.display %>% filter(PC1 < .0) # this is the AT1 AT2 basel population
c3 <- c3 %>% filter(PC2 < -.025)

#colnames(converter_df) <- c('mouse_uni','shared_gene','human_uni')
#write.csv(converter_df,'/Users/andrewleduc/Desktop/convert_uniprot_mouse_human.csv')


#again get data in format that allows us to easily do diff exp

data <- melt(t7_filt)
data$Condition <- NA
data$Condition[data$Var2 %in% c1$id] <- 'c1'
data$Condition[data$Var2 %in% c2$id] <- 'c2'
data$Condition[data$Var2 %in% c3$id] <- 'c3'
data <-data %>% filter(is.na(Condition)==F)


colnames(data)<- c('Protein','SampID','Intensity','Condition')

converter_df <- read.csv('convert_uniprot_mouse_human.csv')
GO_db <- read.delim('GOA_db.txt',sep = " ")


## GSEA code
GSEA_out <- Population_GSEA(data,GO_db,3,1,converter_df)
GSEA_out <- GSEA_out %>% filter(is.na(pVal)==F)
GSEA_out$qval <- p.adjust(GSEA_out$pVal,method = 'BH')
GSEA_out <- GSEA_out %>% filter(qval < .01)
GSEA_out_ <-GSEA_out %>% filter(numberOfMatches > 3)



diff_prot_output <- data.frame(GO_term = character(), pVal = numeric(),numberOfMatches= numeric(),fractionOfDB_Observed = numeric(),stringsAsFactors = FALSE,
                          Cond1med_int = numeric(),Cond2med_int = numeric(),Cond3med_int = numeric())





