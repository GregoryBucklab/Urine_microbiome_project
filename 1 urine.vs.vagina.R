########### library ###########
library(stringr)
library(GUniFrac)
library(vegan)
library(ggplot2)
library(ALDEx2)
library(tidyr)
library(DiscriMiner)
library(compositions)

dif_abundance <- function(reads_table,metadata, pvalue_th = 0.05, fold_change_th = 1, paired_test = F, order_reverse = F, style = 1, order = NA) {
   # style =1
   # reads_table= cts
   # metadata=as.character(metadata_virus_PTB$Delivery)
   # paired_test = F
   # order_reverse = F
   # style =2 
   # order = c('TB','PTB')
   # fold_change_th = 1
   
   conds <- metadata
   
   x <- aldex.clr(reads_table, conds, mc.samples=128, denom="all", verbose=F)
   
   # paired Wilcoxon Rank Sum test and Welch's t-test
   x.tt <- aldex.ttest(x, paired.test= paired_test)
   
   if (sum(x.tt$wi.eBH == 'NaN') > 5) {
      x.tt$wi.eBH = x.tt$we.eBH
      print('we.eBH instead of wi.eBH')
   }
   
   if (min(x.tt$wi.eBH) > 0.05) {
      print('No taxon has significant abundance change')
      return(NA)
   }
   
   x.effect <- aldex.effect(x)
   
   x.all <- data.frame(cbind(x.tt,x.effect))
   
   abundance_change_das <- x.all$diff.btw # diff.btw is a vector containing the per-feature median difference between condition A and B
   
   if (max(abs(abundance_change_das)) < fold_change_th) {
      print('No taxon has significant abundance change')
      return(NA)
   }
   
   if (order_reverse == T) {
      abundance_change_das = -abundance_change_das
   }
   
   x.all <- cbind(abundance_change_das,x.all)
   
   `-Log10(adj-pvalue)` <- -log10(x.all$wi.eBH)     # use wi.eBH as the adj-pvalue
   #x.all$abundance_das <- abundance_das
   x.all$`-Log10(adj-pvalue)` <- `-Log10(adj-pvalue)`
   
   # draw figure
   das <- x.all[(x.all$`-Log10(adj-pvalue)` >= -log10(pvalue_th) & (x.all$abundance_change_das >=fold_change_th | x.all$abundance_change_das <=-fold_change_th)),]
   
   if (nrow(das)==0) {
      print('No taxon has significant abundance change')
      return(NA)
   }
   
   das$Species <- row.names(das)
   das <- das[order(das$abundance_change_das),] 
   
   metadata = as.factor(metadata)
   lev = levels(metadata)
   
   if (order_reverse == T) {
      das$Color <- ifelse(das$abundance_change_das < 0, paste0("Enriched in ",lev[2]), paste0("Enriched in ",lev[1]))  # above / below avg flag
   } else {
      das$Color <- ifelse(das$abundance_change_das < 0, paste0("Enriched in ",lev[1]), paste0("Enriched in ",lev[2]))  # above / below avg flag
      
   }
   
   das$Species <- factor(das$Species, levels = das$Species)  # convert to factor to retain sorted order in plot.
   
   das$abundance_change_das[das$abundance_change_das == Inf] = 10
   das$abundance_change_das[das$abundance_change_das == -Inf] = -10
   das$abundance_change_das[das$abundance_change_das <= -10] = -10
   das$abundance_change_das[das$abundance_change_das >= 10] = 10
   
   
   if (style == 1) {
      theme_set(theme_bw())  
      
      p <- ggplot(das, aes(Species, abundance_change_das)) + 
         geom_point(aes(col=Color, size=`-Log10(adj-pvalue)`)) + 
         coord_flip() +          # convert x y axis
         labs(x = 'Taxa', y = "Median difference in clr values")+ 
         theme(axis.title = element_text(size = 7), 
               axis.text = element_text(size = 7), 
               legend.text = element_text(size = 7), 
               legend.title = element_text(size = 7))  
      
      
   } else {
      taxa_list = row.names(das)
      
      keep = which(row.names(reads_table) %in% taxa_list)
      
      reads_table_abundance = get_abundance_table(reads_table, mc.cores = 8)
      reads_table_2 <- reads_table_abundance[keep,]
      reads_table_2 <- as.data.frame(t(reads_table_2))
      
      #   colnames(reads_table_2)=taxa_list
      reads_table_3 = gather(reads_table_2)
      reads_table_3$Type = rep(conds, length(taxa_list))
      colnames(reads_table_3) = c('Taxa', 'Abundance','Type')
      
      reads_table_3$Taxa = str_replace_all(reads_table_3$Taxa,'.*__','')
      
      if (!is.na(order)[1]) {
         reads_table_3$Type = factor(reads_table_3$Type, levels = order)
      }
      
      p <- ggplot(reads_table_3, aes(x=Taxa, y=Abundance,fill=Type)) +
         geom_boxplot(outlier.shape = NA)+
         ylab("Abundance (%)")+
         theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust=1))
   }
   return(c(list(p = p), list(data = x.all)))
   
}


###################### input raw data ##################################################
# get reads of vaginal and urine samples
setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/Chris')
metadata_v = read.table('data_vagina.txt',header = F, sep = "|")
metadata_u = read.table('data_urine.txt',header = F, sep = "|")

setwd('/Users/binzhu/secure/momspi-data/home/metadata/urine_project/')
metadata_all = read.table('metadata.txt',header = T)
sample_pair = read.table('MU1W2MV1D_urineChris_022520.txt',header = F)
   
setwd('/Users/binzhu/Desktop/Urine')

###################### data reform to a reads table & get total reads, control and metadata of participants ###########################################
### distribution of samples
data = data.frame(pID = unique(metadata_all$PID), number = NA)
data$number = sapply(1:nrow(data), function(j) (data$number[j] = sum(metadata_all$PID == data$pID[j])))

data2 = data.frame(pID = unique(data$number), number = NA)
data2$number = sapply(1:nrow(data2), function(j) (data2$number[j] = sum(data$number == data2$pID[j])))
colnames(data2) = c('Sample number','Participants number')
data2$`Case number`
ggplot(data2, aes(`Sample number`, `Participants number`)) + 
   geom_bar(stat="identity")
ggsave('Case_number.pdf',width=2, height=3)


###
metadata_v=data_v
metadata_u=data_u

sum(metadata_v$V4)
sum(metadata_u$V4)

metadata_v$V3 = str_replace_all(metadata_v$V3,'AT','')
metadata_u$V3 = str_replace_all(metadata_u$V3,'AT','')

metadata_v_2 = as.data.frame(str_c(metadata_v$V2,metadata_v$V3,sep = " "))
metadata_v = cbind(metadata_v$V1,metadata_v_2,metadata_v$V4)

rm(metadata_v_2)

sample_name <- unique(metadata_v$`metadata_v$V1`)
species_name <- unique(metadata_v$`str_c(metadata_v$V2, metadata_v$V3, sep = " ")`)

reads_table_v <- matrix(0, ncol = length(sample_name), nrow = length(species_name))     # create reads table   sample name as columns and species name as rows
row.names(reads_table_v) = species_name 
colnames(reads_table_v) = sample_name

for (a in 1: dim(metadata_v)[1]) {
   
   column_num <- which(sample_name == metadata_v[a,1])
   row_num <- which(species_name == metadata_v[a,2])
   
   reads_table_v[row_num,column_num] = as.numeric(as.character(metadata_v[a,3]))
   
}
reads_table_v = as.data.frame(reads_table_v)

# get reads of urine samples

metadata_u_2 = as.data.frame(str_c(metadata_u$V2,metadata_u$V3,sep = " "))
metadata_u = cbind(metadata_u$V1,metadata_u_2,metadata_u$V4)
rm(metadata_u_2)

sample_name <- unique(metadata_u$`metadata_u$V1`)
species_name <- unique(metadata_u$`str_c(metadata_u$V2, metadata_u$V3, sep = " ")`)

reads_table_u <- matrix(0, ncol = length(sample_name), nrow = length(species_name))     # create reads table   sample name as columns and species name as rows
row.names(reads_table_u) = species_name 
colnames(reads_table_u) = sample_name

for (a in 1: dim(metadata_u)[1]) {
   
   column_num <- which(sample_name == metadata_u[a,1])
   row_num <- which(species_name == metadata_u[a,2])
   
   reads_table_u[row_num,column_num] = as.numeric(as.character(metadata_u[a,3]))
}

reads_table_u = as.data.frame(reads_table_u)


sample_pair$V1 = str_replace_all(sample_pair$V1, 'MU1W', 'MU1D')
reads_table_u = reads_table_u[sample_pair$V1]
reads_table_v = reads_table_v[sample_pair$V2]

metadata_all = metadata_all[match(colnames(reads_table_v), metadata_all$sampleID),]
metadata_all$urine_reads = colSums(reads_table_u)
metadata_all$vagina_reads = colSums(reads_table_v)

### remove duplicated participants
{
   #   keep = colSums(reads_table_v) >=10000 & colSums(reads_table_u) >=10000
   #   keep = colnames(reads_table_v)[keep]
   #   metadata_all = metadata_all[metadata_all$sampleID %in% keep,]
   #   
   #   metadata_all = metadata_all[order(metadata_all$urine_reads, decreasing = T),]
   #   metadata_all = metadata_all[order(metadata_all$PID, decreasing = F),]
   #   metadata_all = metadata_all[!duplicated(metadata_all$PID),]
   #   
   #   keep = colnames(reads_table_v) %in% metadata_all$sampleID
   #   sum(keep)
   #   reads_table_v = reads_table_v[,keep]
   #   reads_table_u = reads_table_u[,keep]
}


colnames(reads_table_v) = paste0('Vagina_',c(1:ncol(reads_table_v)))
colnames(reads_table_u) = paste0('Urine_',c(1:ncol(reads_table_v)))
metadata_all$sampleID = paste0('Sample_',c(1:ncol(reads_table_v)))

metadata_all$PID = as.numeric(factor(metadata_all$PID))
metadata_all$PID = paste0('Participant_',metadata_all$PID)
write.csv(metadata_all,'metadata_all.csv')



### total reads rarefaction curve
rare_v <- matrix(data = 0, ncol =2, nrow = 100)
m=0
for (b in seq(1000, 100000, by=1000)) {
   keep <- matrix(ncol = dim(reads_table_v)[2],)
   for (a in 1: dim(reads_table_v)[2]) {
      keep[,a] = sum(reads_table_v[,a]) >= b    
   }
   m = m+1
   rare_v[m,2] = sum(keep)
   rare_v[m,1] = b
}
rare_v <- as.data.frame(rare_v)
colnames(rare_v) <- c('Total reads threshold','Sample number')
rare_v$Site = rep('Vagina', times = nrow(rare_v))


rare_u <- matrix(data = 0, ncol =2, nrow = 100)
m=0
for (b in seq(1000, 100000, by=1000)) {
   keep <- matrix(ncol = dim(reads_table_u)[2],)
   for (a in 1: dim(reads_table_u)[2]) {
      keep[,a] = sum(reads_table_u[,a]) >= b    
   }
   m = m+1
   rare_u[m,2] = sum(keep)
   rare_u[m,1] = b
}
rare_u <- as.data.frame(rare_u)
colnames(rare_u) <- c('Total reads threshold','Sample number')
rare_u$Site = rep('Urine', times = nrow(rare_u))

rare <- rbind(rare_v,rare_u)

ggplot(data=rare, aes(x=`Total reads threshold`, y=`Sample number`, group = Site,color = Site)) +
   geom_line()+
   geom_point()

ggsave('Total_reads_threshold.pdf', width=4, height=3)


### sample total reads threshold
keep <- matrix(ncol = dim(reads_table_v)[2],)
for (a in 1: dim(reads_table_v)[2]) {
   keep[,a] = sum(reads_table_v[,a]) >=10000    # input
}
reads_table_u <- reads_table_u[,keep]
reads_table_v <- reads_table_v[,keep]
metadata_all = metadata_all[keep,]

keep <- matrix(ncol = dim(reads_table_u)[2],)
for (a in 1: dim(reads_table_u)[2]) {
   keep[,a] = sum(reads_table_u[,a]) >=10000    # input
}

reads_table_u <- reads_table_u[,keep]
reads_table_v <- reads_table_v[,keep]
metadata_all = metadata_all[keep,]

### prepare two reads tables
species_list = unique(c(row.names(reads_table_v),row.names(reads_table_u)))
species_list <- species_list[order(species_list)]
reads_table_u <- reads_table_u[order(row.names(reads_table_u)), ]
reads_table_v <- reads_table_v[order(row.names(reads_table_v)), ]

reads_table_u_new <- as.data.frame(matrix(data=0,ncol=ncol(reads_table_u),nrow=length(species_list)))
reads_table_v_new <- as.data.frame(matrix(data=0,ncol=ncol(reads_table_u),nrow=length(species_list)))
row.names(reads_table_u_new) = species_list
colnames(reads_table_u_new) = colnames(reads_table_u)
row.names(reads_table_v_new) = species_list
colnames(reads_table_v_new) = colnames(reads_table_v)


keep <- species_list %in% row.names(reads_table_v)
reads_table_v_new[keep,] = reads_table_v

keep <- species_list %in% row.names(reads_table_u)
reads_table_u_new[keep,] = reads_table_u

reads_table_all <- cbind(reads_table_v_new,reads_table_u_new)



##### difference associated with ga age #####
reads_table = t(reads_table_v)
reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
reads_table <- reads_table$otu.tab.rff
reads_table <- as.data.frame(reads_table)

pvalue <- adonis2(reads_table ~ ga_at_col, data = metadata_all, method = "bray")  
pvalue

reads_table = t(reads_table_u)
reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
reads_table <- reads_table$otu.tab.rff
reads_table <- as.data.frame(reads_table)

pvalue <- adonis2(reads_table ~ ga_at_col, data = metadata_all, method = "bray")  
pvalue


##### STIRRUPS no th #####
# get samples metadata
sample_list_v = as.matrix(colnames(reads_table_v_new))
sample_list_u = as.matrix(colnames(reads_table_u_new))

sample_list = as.data.frame(matrix(data = NA, ncol=3, nrow= 2*nrow(sample_list_u)))
colnames(sample_list) = c('Sample','Pair','Site')
sample_list$Sample = rbind(sample_list_v,sample_list_u)
sample_list$Pair = rbind(sample_list_v,sample_list_v)
sample_list$Site = c(matrix(data = 'Vagina' , nrow = nrow(sample_list_u)),matrix(data = 'Urine' , nrow = nrow(sample_list_u)))

reads_table_all = reads_table_all[!(str_detect(row.names(reads_table_all),'No Hit')),]
urine_vagna_analyses(sample_list, reads_table_all, reads_table_u_new , reads_table_v_new)
data = dif_abundance(reads_table_all,sample_list$Site, pvalue_th = 0.05, fold_change_th = 1, paired_test = T, order_reverse = F, style = 1, order = NA)
data$p
ggsave("das.pdf", width=5, height=3)
ggsave("das2.pdf", width=8, height=4)

reads_table = t(reads_table_all)
reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
reads_table <- reads_table$otu.tab.rff
reads_table <- as.data.frame(reads_table)

metadata = rbind(metadata_all,metadata_all)
metadata$site = c(rep('vagina', 75), rep('urine', 75))
pvalue <- adonis2(reads_table ~ ga_at_col*site, data = metadata, method = "bray")  
pvalue

#####  present ###################################
reads_table = reads_table_all
keep = str_detect(row.names(reads_table),'Hit')  # remove a taxa named 'No hit'
reads_table = reads_table[!keep,]

reads_table_v_new <- as.data.frame(reads_table[,(1:(ncol(reads_table)/2))])   # split the U and V combined reads table to U and V tables
reads_table_u_new <- as.data.frame(reads_table[,(ncol(reads_table)/2+1) : ncol(reads_table)])

reads_table_abundance <- matrix(data =0, ncol = ncol(reads_table),nrow = nrow(reads_table))  # get abundance table

for (a in 1:ncol(reads_table)) {
   reads_table_abundance[,a] <- reads_table[,a] / colSums(reads_table)[a]
}
reads_table_abundance <- as.data.frame(reads_table_abundance)
row.names(reads_table_abundance) = row.names(reads_table)
colnames(reads_table_abundance) = colnames(reads_table)

# get samples metadata
sample_list_v = as.matrix(colnames(reads_table_v_new))
sample_list_u = as.matrix(colnames(reads_table_u_new))

sample_list = as.data.frame(matrix(data = NA, ncol=3, nrow= 2*nrow(sample_list_u)))
colnames(sample_list) = c('Sample','Pair','Site')
sample_list$Sample = rbind(sample_list_v,sample_list_u)
sample_list$Pair = rbind(sample_list_v,sample_list_v)
sample_list$Site = c(matrix(data = 'Vagina' , nrow = nrow(sample_list_u)),matrix(data = 'Urine' , nrow = nrow(sample_list_u)))


# find species in two sites # any taxon with a relative abundance no less than 0.01% in a sample was considered present in a sample.
present <- as.data.frame(matrix(data = NA, nrow = nrow(reads_table_u_new), ncol = ncol(reads_table_u_new)))
colnames(present) <- colnames(reads_table_u_new)
row.names(present) <- row.names(reads_table_u_new)

for (a in 1: nrow(present)) {
   for (b in 1: ncol(present)) {
      c = reads_table_abundance[a,b] > 0.0001    # input
      d = reads_table_abundance[a,b+(ncol(reads_table_u_new))] > 0.0001    # input
      
      if (c == T & d == T) {
         present[a,b] = 'Both'
      } else if (c == T & d == F) {
         present[a,b] = 'Vagina'
      } else if (c == F & d == T) {
         present[a,b] = 'Urine'
      } else {
         present[a,b] = 'None'
      }
   }
}

# change the format of the present table for plotting
present <- as.data.frame(t(present))

present_new <- as.data.frame(matrix(data =0, nrow = 4, ncol = ncol(present)))
colnames(present_new) <- colnames(present)
row.names(present_new) <- c('Both','Vagina','Urine','None')

present_new_2 <- as.data.frame(matrix(data =0, nrow = 3, ncol = ncol(present)))
colnames(present_new_2) <- colnames(present)
row.names(present_new_2) <- c('Vagina','Urine','Ratio')


for (a in 1: nrow(present)) {
   for (b in 1: ncol(present)) {
      if (present[a,b] == 'Both') {
         present_new[1,b] = present_new[1,b]+1
         present_new_2[1,b] = present_new_2[1,b]+1
         present_new_2[2,b] = present_new_2[2,b]+1
      } else if (present[a,b] == 'Vagina') {
         present_new[2,b] = present_new[2,b]+1
         present_new_2[1,b] = present_new_2[1,b]+1
      } else if (present[a,b] == 'Urine') {
         present_new[3,b] = present_new[3,b]+1
         present_new_2[2,b] = present_new_2[2,b]+1
      } else {
         present_new[4,b] = present_new[4,b]+1
      }
   }
}
present_new <- as.data.frame(t(present_new))
present_new <- present_new[order(present_new$Both),]
present_new_name <- row.names(present_new)

present_new <- gather(present_new)

present_name <- rep( present_new_name , 4)
present_new$Species <- present_name

present_new$Species<- as.factor(present_new$Species)

present_new$Species <- factor(present_new$Species, levels=present_new_name)

colnames(present_new)[1] = 'Present'

ggplot(data=present_new, aes(x=Species, y=value, fill = `Present`)) +
   geom_bar(stat="identity") +
   labs(x = 'Species', y = "Number of sample pairs")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 1), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16)
   ) + coord_flip()+ scale_fill_manual(values=c("#9F9F9F", "#FFFFFF", "#000000",'#D6D6D6'))

ggsave('present.pdf', width=3.5, height=8)

present_new_2 = as.data.frame(t(present_new_2))
present_new_2$Ratio = present_new_2$Urine / present_new_2$Vagina
write.csv(present_new_2,'present.csv', quote = F)

# detect outliers
out = boxplot.stats(present_new_2$Ratio[present_new_2$Ratio != Inf & !is.na(present_new_2$Ratio)])$out
out <- which(present_new_2$Ratio %in% out)
out <- row.names(present_new_2)[out]
write.table(out,'outlier_boxplot.stats.txt', quote = F, sep = '\t')

sum(rowSums(reads_table_v)[row.names(reads_table_v) %in% out]) / sum(rowSums(reads_table_v))

###################### STIRRUPS with taxa th ###################################
reads_table = reads_table_all
keep = str_detect(row.names(reads_table),'Hit')  # remove a taxa named 'No hit'
reads_table = reads_table[!keep,]

# get abundance table
reads_table_abundance <- matrix(data =0, ncol = ncol(reads_table),nrow = nrow(reads_table))

for (a in 1:ncol(reads_table)) {
   reads_table_abundance[,a] <- reads_table[,a] / colSums(reads_table)[a]
}
reads_table_abundance <- as.data.frame(reads_table_abundance)
row.names(reads_table_abundance) = row.names(reads_table)
colnames(reads_table_abundance) = colnames(reads_table)

# species threshold
keep <- matrix(, ncol = nrow(reads_table_abundance))

for (a in 1: nrow(reads_table_abundance)) {
   c = sum(reads_table_abundance[a,] >= 0.001) / ncol(reads_table_abundance) >= 0.05  # only those with at least 0.1% (or 0.01%) relative abundance in at least 5% (or 15%), respectively, of the samples were kept in the analysis
   
   d = sum(reads_table_abundance[a,] >= 0.0001) / ncol(reads_table_abundance) >= 0.15
   keep[a] = c|d
}
sum(keep)
reads_table <- reads_table[keep,]
write.csv(reads_table,'STIRRUPs_taxa_th.csv')

reads_table_v_new <- as.data.frame(reads_table[,(1:(ncol(reads_table)/2))])   # split the U and V combined reads table to U and V tables
reads_table_u_new <- as.data.frame(reads_table[,(ncol(reads_table)/2+1) : ncol(reads_table)])


# get samples metadata
sample_list_v = as.matrix(colnames(reads_table_v_new))
sample_list_u = as.matrix(colnames(reads_table_u_new))

sample_list = as.data.frame(matrix(data = NA, ncol=3, nrow= 2*nrow(sample_list_u)))
colnames(sample_list) = c('Sample','Pair','Site')
sample_list$Sample = rbind(sample_list_v,sample_list_u)
sample_list$Pair = rbind(sample_list_v,sample_list_v)
sample_list$Site = c(matrix(data = 'Vagina' , nrow = nrow(sample_list_u)),matrix(data = 'Urine' , nrow = nrow(sample_list_u)))

### run ###
urine_vagna_analyses(sample_list, reads_table, reads_table_u_new , reads_table_v_new)
data = dif_abundance(reads_table,sample_list$Site, pvalue_th = 0.05, fold_change_th = 1, paired_test = T, order_reverse = F, style = 1, order = NA)
data$p
data2 = data$data
ggsave("das.pdf", width=5, height=3)
ggsave("das2.pdf", width=8, height=4)

##### remove outlier taxa and no taxa th #####
reads_table = reads_table_all
keep = row.names(reads_table_all) %in% out
sum(keep)
reads_table = reads_table[!keep,]

write.csv(reads_table , 'STIRRUPs_rm_outliers.csv')
reads_table_v_new <- as.data.frame(reads_table[,(1:(ncol(reads_table)/2))])   # split the U and V combined reads table to U and V tables
reads_table_u_new <- as.data.frame(reads_table[,(ncol(reads_table)/2+1) : ncol(reads_table)])

# get samples metadata
sample_list_v = as.matrix(colnames(reads_table_v_new))
sample_list_u = as.matrix(colnames(reads_table_u_new))

sample_list = as.data.frame(matrix(data = NA, ncol=3, nrow= 2*nrow(sample_list_u)))
colnames(sample_list) = c('Sample','Pair','Site')
sample_list$Sample = rbind(sample_list_v,sample_list_u)
sample_list$Pair = rbind(sample_list_v,sample_list_v)
sample_list$Site = c(matrix(data = 'Vagina' , nrow = nrow(sample_list_u)),matrix(data = 'Urine' , nrow = nrow(sample_list_u)))

urine_vagna_analyses(sample_list, reads_table, reads_table_u_new , reads_table_v_new)
data = dif_abundance(reads_table,sample_list$Site, pvalue_th = 0.05, fold_change_th = 1, paired_test = T, order_reverse = F, style = 1, order = NA)
data$p
data2 = data$data
ggsave("das.pdf", width=5, height=1.2)
ggsave("das2.pdf", width=8, height=4)













##### remove outlier taxa and taxa th #####
reads_table = reads_table_all
keep = str_detect(row.names(reads_table),'Hit')  # remove a taxa named 'No hit'
reads_table = reads_table[!keep,]

# get abundance table
reads_table_abundance <- matrix(data =0, ncol = ncol(reads_table),nrow = nrow(reads_table))

for (a in 1:ncol(reads_table)) {
   reads_table_abundance[,a] <- reads_table[,a] / colSums(reads_table)[a]
}
reads_table_abundance <- as.data.frame(reads_table_abundance)
row.names(reads_table_abundance) = row.names(reads_table)
colnames(reads_table_abundance) = colnames(reads_table)

# species threshold
keep <- matrix(, ncol = nrow(reads_table_abundance))

for (a in 1: nrow(reads_table_abundance)) {
   c = sum(reads_table_abundance[a,] >= 0.001) / ncol(reads_table_abundance) >= 0.05  # only those with at least 0.1% (or 0.01%) relative abundance in at least 5% (or 15%), respectively, of the samples were kept in the analysis
   
   d = sum(reads_table_abundance[a,] >= 0.0001) / ncol(reads_table_abundance) >= 0.15
   keep[a] = c|d
}
sum(keep)
reads_table <- reads_table[keep,]

keep = row.names(reads_table) %in% out
sum(keep)
write.table(row.names(reads_table)[keep],'outlier_boxplot.stats_rm_after_taxa_th.txt', quote = F, sep = '\t')

reads_table = reads_table[!keep,]
write.csv(reads_table , 'STIRRUPs_taxa_th_rm_outliers.csv')

reads_table_v_new <- as.data.frame(reads_table[,(1:(ncol(reads_table)/2))])   # split the U and V combined reads table to U and V tables
reads_table_u_new <- as.data.frame(reads_table[,(ncol(reads_table)/2+1) : ncol(reads_table)])

# get samples metadata
sample_list_v = as.matrix(colnames(reads_table_v_new))
sample_list_u = as.matrix(colnames(reads_table_u_new))

sample_list = as.data.frame(matrix(data = NA, ncol=3, nrow= 2*nrow(sample_list_u)))
colnames(sample_list) = c('Sample','Pair','Site')
sample_list$Sample = rbind(sample_list_v,sample_list_u)
sample_list$Pair = rbind(sample_list_v,sample_list_v)
sample_list$Site = c(matrix(data = 'Vagina' , nrow = nrow(sample_list_u)),matrix(data = 'Urine' , nrow = nrow(sample_list_u)))

urine_vagna_analyses(sample_list, reads_table, reads_table_u_new , reads_table_v_new)
data = dif_abundance(reads_table,sample_list$Site, pvalue_th = 0.05, fold_change_th = 1, paired_test = T, order_reverse = F, style = 1, order = NA)
data$p
data2 = data$data
ggsave("das.pdf", width=5, height=1)
ggsave("das2.pdf", width=8, height=4)












##### lm compare #####
data = read.csv('linear_regression_all.csv', header = T)
data = data[,c(4,10,16,22,28)]
data = data[-1,]

keep = 30 # top n abundant taxa 
data = data[c(1:keep),]

median(as.numeric(data[,1]))
median(as.numeric(data[,2]))
median(as.numeric(data[,3]))
median(as.numeric(data[,4]))
median(as.numeric(data[,5]))


data_2 = reads_table_all
data_2 = rowSums(data_2)
data_2 = order(data_2, decreasing = T)
data_2 = data_2[1:keep]
sum(rowSums(reads_table_all)[data_2]) / sum(rowSums(reads_table_all))

plot_data <- gather(data)
plot_data = plot_data[plot_data$value!="",]
plot_data$value = as.numeric(plot_data$value)
plot_data$key = factor(plot_data$key, levels = colnames(data))

ggplot(plot_data, aes(x=key, y=value)) +
   geom_boxplot() +
   geom_jitter(shape=16, position=position_jitter(0.2))+
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16),
         axis.text.x = element_text(angle = 60, hjust = 1))
ggsave('Rvalue_of_linear_regression_comparison.pdf', width=3, height=6)

# significance of group sample distance, adonis test
group_level = unique(plot_data$key)
n = length(group_level)

group_dis_sig <- matrix(data = NA, ncol=n, nrow = n)
colnames(group_dis_sig) = group_level
row.names(group_dis_sig) = group_level

for (a in 1:(n-1)) {
   for (b in (a+1):n) {
      keep = plot_data$key == group_level[a] | plot_data$key == group_level[b]
      plot_data_2 = plot_data[keep,]
      pvalue <- wilcox.test(value ~ key,plot_data_2)$p.value
      group_dis_sig[a,b] = pvalue
      group_dis_sig[b,a] = pvalue
   }
}
write.csv(group_dis_sig,'Rvalue_of_linear_regression_comparison.csv', row.names = T, quote = F)


data = read.csv('linear_regression_all.csv', header = T)
data = data[,c(2,8,14,20,26)]
data = data[-1,]

data = sapply(1:ncol(data), function(j) (data[,j] = -log10(as.numeric(data[,j]))))
data = as.data.frame(data)

keep = 30 # top n abundant taxa 
data = data[c(1:keep),]

median(as.numeric(data[,1]))
median(as.numeric(data[,2]))
median(as.numeric(data[,3]))
median(as.numeric(data[,4]))
median(as.numeric(data[,5]))

data_2 = reads_table_all
data_2 = rowSums(data_2)
data_2 = order(data_2, decreasing = T)
data_2 = data_2[1:keep]
sum(rowSums(reads_table_all)[data_2]) / sum(rowSums(reads_table_all))

plot_data <- gather(data)
plot_data = plot_data[plot_data$value!="",]
plot_data$value = as.numeric(plot_data$value)
plot_data$key = factor(plot_data$key, levels = colnames(data))

ggplot(plot_data, aes(x=key, y=value)) +
   geom_boxplot() +
   geom_jitter(shape=16, position=position_jitter(0.2))+
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16),
         axis.text.x = element_text(angle = 60, hjust = 1))
ggsave('Pvalue_of_linear_regression_comparison.pdf', width=3, height=6)

# significance of group sample distance, adonis test
group_level = unique(plot_data$key)
n = length(group_level)

group_dis_sig <- matrix(data = NA, ncol=n, nrow = n)
colnames(group_dis_sig) = group_level
row.names(group_dis_sig) = group_level

for (a in 1:(n-1)) {
   for (b in (a+1):n) {
      keep = plot_data$key == group_level[a] | plot_data$key == group_level[b]
      plot_data_2 = plot_data[keep,]
      pvalue <- wilcox.test(value ~ key,plot_data_2)$p.value
      group_dis_sig[a,b] = pvalue
      group_dis_sig[b,a] = pvalue
   }
}
write.csv(group_dis_sig,'Pvalue_of_linear_regression_comparison.csv', row.names = T, quote = F)

##### lm genus level #####
reads_table = reads_table_all
keep = row.names(reads_table_all) %in% out
sum(keep)
reads_table = reads_table[!keep,]
row.names(reads_table) = str_replace_all(row.names(reads_table), ' ','_')

# convert to genus level  
genus_list = row.names(reads_table)
genus_list = str_replace(genus_list,'Lactobacillus_iners','Lactobacillusxxxiners')
genus_list = str_replace(genus_list,'Lactobacillus_jensenii','Lactobacillusxxxjensenii')
genus_list = str_replace(genus_list,'Lactobacillus_crispatus_cluster','Lactobacillusxxxcrispatus')
genus_list = str_replace(genus_list,'Lactobacillus_gasseri_cluster','Lactobacillusxxxgasseri')
genus_list = str_replace_all(genus_list, '_.*','')
genus_list[genus_list == "Lactobacillus"] = "Other_Lactobacillus"
genus_list = str_replace(genus_list,'xxx','_')

genus_list_2 = unique(genus_list)
genus_list_2 = genus_list_2[-which(genus_list_2 == 'No')]

reads_table_all_genus = as.data.frame(matrix(data = 0, ncol = ncol(reads_table), nrow = length(genus_list_2)))
colnames(reads_table_all_genus) = colnames(reads_table)
row.names(reads_table_all_genus) = genus_list_2

for (a in 1: nrow(reads_table_all_genus)) {
   n = which(genus_list == genus_list_2[a])
   reads_table_all_genus[a,] = colSums(reads_table[which(genus_list == genus_list_2[a]),])
}

out_liers = read.csv('outlier_boxplot.stats_genus.csv')
keep = row.names(reads_table_all_genus) %in% out_liers$out_genus
sum(keep)
reads_table_all_genus = reads_table_all_genus[!keep,]
reads_table = reads_table_all_genus

reads_table_v_new <- as.data.frame(reads_table[,(1:(ncol(reads_table)/2))])   # split the U and V combined reads table to U and V tables
reads_table_u_new <- as.data.frame(reads_table[,(ncol(reads_table)/2+1) : ncol(reads_table)])

# get samples metadata
sample_list_v = as.matrix(colnames(reads_table_v_new))
sample_list_u = as.matrix(colnames(reads_table_u_new))

sample_list = as.data.frame(matrix(data = NA, ncol=3, nrow= 2*nrow(sample_list_u)))
colnames(sample_list) = c('Sample','Pair','Site')
sample_list$Sample = rbind(sample_list_v,sample_list_u)
sample_list$Pair = rbind(sample_list_v,sample_list_v)
sample_list$Site = c(matrix(data = 'Vagina' , nrow = nrow(sample_list_u)),matrix(data = 'Urine' , nrow = nrow(sample_list_u)))

urine_vagna_analyses(sample_list, reads_table, reads_table_u_new , reads_table_v_new)
