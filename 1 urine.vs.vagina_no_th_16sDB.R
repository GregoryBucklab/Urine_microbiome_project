########### library ###########
library(stringr)
library(GUniFrac)
library(vegan)
library(ggplot2)
library(ALDEx2)
library(tidyr)
library(DiscriMiner)
library(compositions)


###################### prepare samples ##################################################
# get reads of vaginal samples
setwd('/Users/binzhu/secure/momspi-data/home/metadata/urine_project')
metadata_v = read.table('allVaginal.16SDB.V1-V3.uclust99.summary.txt',header = F)
metadata_u = read.table('allUrine.16SDB.V1-V3.uclust99.summary.txt',header = F)

setwd('/Users/binzhu/Desktop/Urine')

data_v = metadata_v
data_u = metadata_u

###################### data arrangement ###########################################
metadata_v = metadata_v[-1,]
metadata_u = metadata_u[-1,]

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


colnames(reads_table_v) = paste0('Vagina_',c(1:ncol(reads_table_v)))
colnames(reads_table_u) = paste0('Urine_',c(1:ncol(reads_table_v)))



##### total reads rarefaction curve ###########################
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


##### sample total reads threshold ####################
keep <- matrix(ncol = dim(reads_table_v)[2],)
for (a in 1: dim(reads_table_v)[2]) {
   keep[,a] = sum(reads_table_v[,a]) >=10000    # input
}
reads_table_u <- reads_table_u[,keep]
reads_table_v <- reads_table_v[,keep]

keep <- matrix(ncol = dim(reads_table_u)[2],)
for (a in 1: dim(reads_table_u)[2]) {
   keep[,a] = sum(reads_table_u[,a]) >=10000    # input
}

reads_table_u <- reads_table_u[,keep]
reads_table_v <- reads_table_v[,keep]

##### prepare two reads tables #########
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

reads_table <- cbind(reads_table_v_new,reads_table_u_new)
write.csv(reads_table,'16S_no_th_u.csv') 

###################### no threshold  ##################################################
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