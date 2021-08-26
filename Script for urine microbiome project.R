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
# get reads of vaginal and urine samples

setwd('/Users/binzhu/secure/godel/gpfs_fs/bccl/bzhu/Chris')
metadata_v = read.table('data_vagina.txt',header = F, sep = "|")
metadata_u = read.table('data_urine.txt',header = F, sep = "|")

setwd('/Users/binzhu/Desktop/Urine')



###################### data reform to a reads table ###########################################
sum(metadata_v$V4)
sum(metadata_u$V4)

metadata_v$V3 = str_replace_all(metadata_v$V3,'AT','')
metadata_u$V3 = str_replace_all(metadata_u$V3,'AT','')

metadata_v_2 = as.data.frame(str_c(metadata_v$V2,metadata_v$V3,sep = " "))
metadata_v = cbind(metadata_v$V1,metadata_v_2,metadata_v$V4)

rm(metadata_v_2)

sample_name <- unique(metadata_v$`metadata_v$V1`)  # get unique sample names
species_name <- unique(metadata_v$`str_c(metadata_v$V2, metadata_v$V3, sep = " ")`)    # get unique taxa names

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

colnames(reads_table_v) = paste0('Vagina_',c(1:84))
colnames(reads_table_u) = paste0('Urine_',c(1:84))

write.csv(reads_table_v,'reads_table_v.csv')
write.csv(reads_table_u,'reads_table_u.csv')

##### test sample number after different total reads threshold ###########################
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


##### sample total reads threshold = 10,000 (this high threshold is for establishing better linear model shown below) ####################
keep <- matrix(ncol = dim(reads_table_v)[2],)
for (a in 1: dim(reads_table_v)[2]) {
   keep[,a] = sum(reads_table_v[,a]) >=10000    # input
}
reads_table_u <- reads_table_u[,keep]   # paired urine and vaginal samples are removed if total reads of vagina < 10,000
reads_table_v <- reads_table_v[,keep]

keep <- matrix(ncol = dim(reads_table_u)[2],)
for (a in 1: dim(reads_table_u)[2]) {
   keep[,a] = sum(reads_table_u[,a]) >=10000    # input
}

reads_table_u <- reads_table_u[,keep]   # paired urine and vaginal samples are removed if total reads of urine < 10,000
reads_table_v <- reads_table_v[,keep]

##### prepare two reads tables #########
species_list = unique(c(row.names(reads_table_v),row.names(reads_table_u)))   # get taxa overlap in two reads table, unique taxa are removed.
species_list <- species_list[order(species_list)]
reads_table_u <- reads_table_u[order(row.names(reads_table_u)), ]
reads_table_v <- reads_table_v[order(row.names(reads_table_v)), ]

reads_table_u_new <- as.data.frame(matrix(data=0,ncol=ncol(reads_table_u),nrow=length(species_list)))
reads_table_v_new <- as.data.frame(matrix(data=0,ncol=ncol(reads_table_u),nrow=length(species_list)))
row.names(reads_table_u_new) = species_list
colnames(reads_table_u_new) = colnames(reads_table_u)
row.names(reads_table_v_new) = species_list
colnames(reads_table_v_new) = colnames(reads_table_v)


keep <- species_list %in% row.names(reads_table_v)   # remove unique taxa in V
reads_table_v_new[keep,] = reads_table_v

keep <- species_list %in% row.names(reads_table_u)   # remove unique taxa in U
reads_table_u_new[keep,] = reads_table_u 

reads_table <- cbind(reads_table_v_new,reads_table_u_new)

###################### no taxa threshold rarefaction ##################################################

reads_table_rare <- t(reads_table)
reads_table_rare = Rarefy(reads_table_rare, depth = min(rowSums(reads_table_rare)))
reads_table_rare <- reads_table_rare$otu.tab.rff
reads_table_rare <- as.data.frame(reads_table_rare)

# get samples metadata
sample_list_v = as.matrix(colnames(reads_table_v_new))
sample_list_u = as.matrix(colnames(reads_table_u_new))

sample_list = as.data.frame(matrix(data = NA, ncol=3, nrow= 2*nrow(sample_list_u)))
colnames(sample_list) = c('Sample','Pair','Site')
sample_list$Sample = rbind(sample_list_v,sample_list_u)
sample_list$Pair = rbind(sample_list_v,sample_list_v)
sample_list$Site = c(matrix(data = 'Vagina' , nrow = nrow(sample_list_u)),matrix(data = 'Urine' , nrow = nrow(sample_list_u)))


##### alpha diversity ##################################################

alpha.shannon_diversity <- data.frame(diversity(reads_table_rare))
alpha.evenness <- alpha.shannon_diversity/log(specnumber(reads_table_rare))
alpha.ovserved_OTU <- data.frame(colSums(t(reads_table_rare) != 0))

sample_list$alpha.shannon <- alpha.shannon_diversity$diversity.reads_table_rare.
sample_list$alpha.evenness <- alpha.evenness$diversity.reads_table_rare.
sample_list$alpha.ovserved_OTU <- alpha.ovserved_OTU$colSums.t.reads_table_rare.....0.

# creat images for alpha diversity
geom_boxplot(outlier.colour="black", outlier.shape=16,
             outlier.size=2, notch=FALSE)

ggplot(sample_list, aes(x=Site, y=alpha.shannon ,color = Site)) +
   geom_boxplot() +
   geom_jitter(shape=16, position=position_jitter(0.2))+
   labs(x = NULL, y = "Shannon index")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16))+ theme_bw()
ggsave('no_th_shannon.pdf', width=3.5, height=3.5)

ggplot(sample_list, aes(x=Site, y=alpha.shannon)) +
   geom_line(aes(group = Pair), size = 0.5, alpha = 0.7) +
   geom_point(aes(color = Site), alpha = 0.5) +
   labs(x = NULL, y = "Shannon index")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16))+ theme_bw()
ggsave("no_th_shannon_2.pdf", width=5.5, height=5.5)

ggplot(sample_list, aes(x=Site, y=alpha.evenness ,color = Site)) +
   geom_boxplot() +
   geom_jitter(shape=16, position=position_jitter(0.2))+
   labs(x = NULL, y = "Evenness")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16))+ theme_bw()
ggsave('no_th_evenness.pdf', width=3.5, height=3.5)

ggplot(sample_list, aes(x=Site, y=alpha.evenness)) +
   geom_line(aes(group = Pair), size = 0.5, alpha = 0.7) +
   geom_point(aes(color = Site), alpha = 0.5) +
   labs(x = NULL, y = "Evenness")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16))+ theme_bw()
ggsave("no_th_evenness_2.pdf", width=5.5, height=5.5)

ggplot(sample_list, aes(x=Site, y=alpha.ovserved_OTU ,color = Site)) +
   geom_boxplot() +
   geom_jitter(shape=16, position=position_jitter(0.2))+
   labs(x = NULL, y = "Observed Species")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16))+ theme_bw()
ggsave('no_th_otu.pdf', width=3.5, height=3.5)

ggplot(sample_list, aes(x=Site, y=alpha.ovserved_OTU)) +
   geom_line(aes(group = Pair), size = 0.5, alpha = 0.7) +
   geom_point(aes(color = Site), alpha = 0.5) +
   labs(x = NULL, y = "Observed species")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16))+ theme_bw()
ggsave("no_th_otu_2.pdf", width=5.5, height=5.5)

# Compute paired Wilcoxon test for the significance of alpha diversity
Urine <- subset(sample_list,  Site == "Urine", alpha.shannon,
                drop = TRUE)
Vagina <- subset(sample_list,  Site == "Vagina", alpha.shannon,
                 drop = TRUE)
no_th_pvalue.alpha.shannon <- wilcox.test(Urine, Vagina, paired = TRUE)$p.value
no_th_pvalue.alpha.shannon

Urine <- subset(sample_list,  Site == "Urine", alpha.evenness,
                drop = TRUE)
Vagina <- subset(sample_list,  Site == "Vagina", alpha.evenness,
                 drop = TRUE)
no_th_pvalue.alpha.evenness <- wilcox.test(Urine, Vagina, paired = TRUE)$p.value
no_th_pvalue.alpha.evenness 

Urine <- subset(sample_list,  Site == "Urine", alpha.ovserved_OTU,
                drop = TRUE)
Vagina <- subset(sample_list,  Site == "Vagina", alpha.ovserved_OTU,
                 drop = TRUE)
no_th_pvalue.alpha.ovserved_OTU <- wilcox.test(Urine, Vagina, paired = TRUE)$p.value
no_th_pvalue.alpha.ovserved_OTU


pvalue = rbind(no_th_pvalue.alpha.shannon,no_th_pvalue.alpha.evenness,no_th_pvalue.alpha.ovserved_OTU)

##### beta diversity ##################################################
# Bray_Curtis distance
Bray_Curtis <- as.matrix(vegdist(reads_table_rare, METHODS="bray", binary=FALSE))
Bray_Curtis <- as.data.frame(Bray_Curtis)

# Bray_Curtis distance between paired samples
x = nrow(sample_list)/2
pair_dis = as.data.frame(matrix(data = NA, nrow = x, ncol =3))
pair_dis$V1 = sample_list[(x+1):(x*2),1]
pair_dis$V2 = sample_list[(x+1):(x*2),2]

for (a in 1: x) {
   bray_x = which(row.names(Bray_Curtis) == as.character(pair_dis[a,1]))
   bray_y = which(colnames(Bray_Curtis) == as.character(pair_dis[a,2]))
   pair_dis[a,3] = Bray_Curtis[bray_x,bray_y]
}
pair_dis <- pair_dis[order(pair_dis$V3),]
colnames(pair_dis) <- c('Sample_name_u', 'Sample_name_v','Bray_Curtis')

pair_dis$Sample_name_u <- factor(pair_dis$Sample_name_u, levels = pair_dis$Sample_name_u)  # convert to factor to retain sorted order in plot.


ggplot(data=pair_dis, aes(x=Sample_name_u, y=Bray_Curtis)) +
   geom_bar(stat="identity") +
   labs(x = 'Paired samples', y = "Bray-Curtis distance")+ 
   theme(axis.title = element_text(size = 12), 
         axis.text.x=element_blank(),
         axis.text.y = element_text(size = 9),
         legend.text = element_text(size = 12), 
         legend.title = element_text(size = 12))   
ggsave("no_th_paired_distance.pdf", width=3.5, height=5.5)


# distance comparison 
vagina_bray_distantc = Bray_Curtis[(1:(nrow(Bray_Curtis)/2)),(1:(ncol(Bray_Curtis)/2))]
vagina_bray_distantc <- gather(vagina_bray_distantc)

distance_comparison <- matrix(data = NA, nrow =(nrow(pair_dis)+nrow(vagina_bray_distantc)), ncol =2)
distance_comparison[1:nrow(pair_dis),1] = pair_dis$Bray_Curtis
distance_comparison[1:nrow(pair_dis),2] = rep('Vagina_Paired urine',nrow(pair_dis))
distance_comparison[(nrow(pair_dis)+1):nrow(distance_comparison),1] = vagina_bray_distantc$value
distance_comparison[(nrow(pair_dis)+1):nrow(distance_comparison),2] = rep('Vagina_Other vagina',nrow(vagina_bray_distantc))
distance_comparison<- as.data.frame(distance_comparison)
colnames(distance_comparison) = c('Bray-Curtis Distance', 'Type')
distance_comparison$`Bray-Curtis Distance` = as.numeric(as.character(distance_comparison$`Bray-Curtis Distance`))
distance_comparison$Type <- as.character(distance_comparison$Type)

ggplot(distance_comparison, aes(x=Type, y=`Bray-Curtis Distance`,fill = Type)) +
   geom_boxplot() +
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16),
         axis.text.x = element_text(angle = 60, hjust = 1))
ggsave('no_th_Bray_curtis_distance_compare.pdf', width=5, height=5)

pvalue.wil1 = wilcox.test(`Bray-Curtis Distance` ~ Type, data=distance_comparison)$p.value 
no_th_pvalue.Bray_Curtis_comparision = pvalue.wil1
pvalue = rbind(pvalue,no_th_pvalue.Bray_Curtis_comparision )

# Running Nonmetric Multidimensional Scaling (NMDS) Ordination
NMDS <-
   metaMDS(Bray_Curtis,
           distance = "bray",
           k = 2,
           maxit = 999, 
           trymax = 50,
           wascores = TRUE)

mds_data <- as.data.frame(NMDS$points)
mds_data$Site <- sample_list$Site
mds_data$Pair <- sample_list$Pair
ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Site)) +
   geom_point()
ggsave("no_th_NMDS1.pdf", width=5.5, height=5.5)

ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Site)) +
   geom_point()+
   geom_line(aes(group = Pair), size = 0.3, alpha = 0.7, color = 'black')
ggsave("no_th_NMDS2.pdf", width=5.5, height=5.5)

# bar plot for beta diversity
vagina2vagina = as.matrix(Bray_Curtis[1:(length(Bray_Curtis)/2),1:(length(Bray_Curtis)/2)])
vagina2vagina = as.vector(vagina2vagina)
vagina2vagina = vagina2vagina[vagina2vagina!=0]
vagina = as.data.frame(cbind(vagina2vagina,'Vagina'))
colnames(vagina) = c('Bray_Curtis_distance','Body_site')

Urine2vagina = as.matrix(Bray_Curtis[1:(length(Bray_Curtis)/2),(length(Bray_Curtis)/2+1):length(Bray_Curtis)])
Urine2vagina = as.vector(Urine2vagina)
urine = as.data.frame(cbind(Urine2vagina,'Urine'))
colnames(urine) = c('Bray_Curtis_distance','Body_site')

beta_plot = rbind(vagina,urine)
beta_plot$Body_site <- as.factor(beta_plot$Body_site)
beta_plot$Bray_Curtis_distance <- as.numeric(as.character(beta_plot$Bray_Curtis_distance))

ggplot(beta_plot, aes(x=Body_site, y=Bray_Curtis_distance,fill=Body_site)) +
   geom_boxplot()  +
   labs(x = NULL, y = "Distance to Vagina samples")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16))

ggsave("no_th_Bray_Curtis_bar.pdf", width=4, height=5.5)

# Compute perMANOVA test for the significance of beta diversity
# adonis2 : permutational multivariate analysis of variance using distance matrices the centroid and/or the spread of the objects is different between the groups.
pvalue.adonis.site <- adonis2(reads_table_rare ~ Site, data = sample_list, method = "bray")  
no_th_pvalue.adonis.site <- pvalue.adonis.site[1,5]
pvalue <- rbind(pvalue,no_th_pvalue.adonis.site)
pvalue.adonis.pair <- adonis2(reads_table_rare ~ Pair, data = sample_list, method = "bray")  
no_th_pvalue.adonis.pair <- pvalue.adonis.pair[1,5]
pvalue <- rbind(pvalue,no_th_pvalue.adonis.pair)

write.csv(pvalue,'pvalue.csv')

# vagitype
vagitype_v = matrix(data = NA, ncol =1, nrow = ncol(reads_table_v_new))
colnames(vagitype_v) = 'Vagitype'

for (a in 1:ncol(reads_table_v_new)) {
   n= which(reads_table_v_new[,a] == max(reads_table_v_new[,a]), arr.ind=TRUE)
   vagitype_v[a,1]= row.names(reads_table_v_new)[n]
}

vagitype_u= matrix(data = NA, ncol =1, nrow = ncol(reads_table_u_new))
colnames(vagitype_u) = 'Vagitype'

for (a in 1:ncol(reads_table_u_new)) {
   n= which(reads_table_u_new[,a] == max(reads_table_u_new[,a]), arr.ind=TRUE)
   vagitype_u[a,1]= row.names(reads_table_u_new)[n]
}

vagitype = rbind(vagitype_v,vagitype_u)
vagitype_unique = unique(vagitype)
for (a in 1:length(vagitype_unique)) {
   if (sum(vagitype == vagitype_unique[a]) < 5) {
      vagitype[vagitype == vagitype_unique[a]] = 'Others'
   }
}

mds_data = cbind(mds_data,vagitype)

mds_data$Vagitype <- as.factor(mds_data$Vagitype)
mds_data$Vagitype <- factor(mds_data$Vagitype, levels = c('Lactobacillus iners ','Enterobacteriaceae cluster31 ',
                                                          'Lachnospiraceae BVAB1 ','Gardnerella vaginalis ',
                                                          'Atopobium vaginae ','Lactobacillus crispatus_cluster ',
                                                          'Lactobacillus gasseri_cluster ','Others'))


ggplot(mds_data, aes(MDS1, MDS2))  +
   geom_line(aes(group = Pair), size = 0.3, alpha = 0.7) +
   geom_point(size=2,aes(shape = Site, color = Vagitype)) +
   scale_color_manual(values=c('lightblue','deepskyblue','orange','red','brown','yellow','pink','grey')) +
   coord_fixed()+ 
   theme(
      axis.title.x = element_text( size=16),
      axis.title.y = element_text( size=16),
      legend.text = element_text(size=16),
      legend.title = element_text(size=16),
      plot.title = element_text(hjust = 0.5, size = 14)
   ) 
ggsave("no_th_vagitype.pdf", width=8, height=5.5)

same_vagitype = matrix(data = NA, nrow= nrow(mds_data)/2, ncol =2)
same_vagitype[,1] = mds_data$Vagitype[1:(nrow(mds_data)/2)]
same_vagitype[,2] = mds_data$Vagitype[(nrow(mds_data)/2+1):nrow(mds_data)]
no_th_number_vagitype_2 = sum(same_vagitype[,1] == same_vagitype[,2])
pvalue <- rbind(pvalue,no_th_number_vagitype_2)

write.csv(pvalue,'pvalue.csv')





###################### species threshold & present ###################################
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

reads_table <- reads_table[keep,]
write.csv(reads_table,'reads_table.csv')

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

for (a in 1: nrow(present)) {
   for (b in 1: ncol(present)) {
      if (present[a,b] == 'Both') {
         present_new[1,b] = present_new[1,b]+1
      } else if (present[a,b] == 'Vagina') {
         present_new[2,b] = present_new[2,b]+1
      } else if (present[a,b] == 'Urine') {
         present_new[3,b] = present_new[3,b]+1
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
         axis.text = element_text(size = 5), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16)
   ) + coord_flip() 
ggsave('present_0.01.pdf', width=7, height=7)

##### rarefaction # that is to normalize total reads of every sample to the same level for diversity analysis ##################################################

reads_table_rare <- t(reads_table)
reads_table_rare = Rarefy(reads_table_rare, depth = min(rowSums(reads_table_rare)))
reads_table_rare <- reads_table_rare$otu.tab.rff
reads_table_rare <- as.data.frame(reads_table_rare)



##### alpha diversity ##################################################

alpha.shannon_diversity <- data.frame(diversity(reads_table_rare))  # see 'vegan' package for diversity analysi
alpha.evenness <- alpha.shannon_diversity/log(specnumber(reads_table_rare))
alpha.ovserved_OTU <- data.frame(colSums(t(reads_table_rare) != 0))

sample_list$alpha.shannon <- alpha.shannon_diversity$diversity.reads_table_rare.
sample_list$alpha.evenness <- alpha.evenness$diversity.reads_table_rare.
sample_list$alpha.ovserved_OTU <- alpha.ovserved_OTU$colSums.t.reads_table_rare.....0.

# creat images for alpha diversity
geom_boxplot(outlier.colour="black", outlier.shape=16,
             outlier.size=2, notch=FALSE)

ggplot(sample_list, aes(x=Site, y=alpha.shannon ,color = Site)) +
   geom_boxplot() +
   geom_jitter(shape=16, position=position_jitter(0.2))+
   labs(x = NULL, y = "Shannon index")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16))+ theme_bw()
ggsave('shannon.pdf', width=3.5, height=3.5)

ggplot(sample_list, aes(x=Site, y=alpha.shannon)) +
   geom_line(aes(group = Pair), size = 0.5, alpha = 0.7) +
   geom_point(aes(color = Site), alpha = 0.5) +
   labs(x = NULL, y = "Shannon index")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16))+ theme_bw()
ggsave("shannon_2.pdf", width=5.5, height=5.5)

ggplot(sample_list, aes(x=Site, y=alpha.evenness ,color = Site)) +
   geom_boxplot() +
   geom_jitter(shape=16, position=position_jitter(0.2))+
   labs(x = NULL, y = "Evenness")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16))+ theme_bw()
ggsave('evenness.pdf', width=3.5, height=3.5)

ggplot(sample_list, aes(x=Site, y=alpha.evenness)) +
   geom_line(aes(group = Pair), size = 0.5, alpha = 0.7) +
   geom_point(aes(color = Site), alpha = 0.5) +
   labs(x = NULL, y = "Evenness")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16))+ theme_bw()
ggsave("evenness_2.pdf", width=5.5, height=5.5)

ggplot(sample_list, aes(x=Site, y=alpha.ovserved_OTU ,color = Site)) +
   geom_boxplot() +
   geom_jitter(shape=16, position=position_jitter(0.2))+
   labs(x = NULL, y = "Observed Species")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16))+ theme_bw()
ggsave('otu.pdf', width=3.5, height=3.5)

ggplot(sample_list, aes(x=Site, y=alpha.ovserved_OTU)) +
   geom_line(aes(group = Pair), size = 0.5, alpha = 0.7) +
   geom_point(aes(color = Site), alpha = 0.5) +
   labs(x = NULL, y = "Observed species")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16))+ theme_bw()
ggsave("otu_2.pdf", width=5.5, height=5.5)

# Compute paired Wilcoxon test for the significance of alpha diversity
Urine <- subset(sample_list,  Site == "Urine", alpha.shannon,
                drop = TRUE)
Vagina <- subset(sample_list,  Site == "Vagina", alpha.shannon,
                 drop = TRUE)
pvalue.alpha.shannon <- wilcox.test(Urine, Vagina, paired = TRUE)$p.value
pvalue.alpha.shannon

Urine <- subset(sample_list,  Site == "Urine", alpha.evenness,
                drop = TRUE)
Vagina <- subset(sample_list,  Site == "Vagina", alpha.evenness,
                 drop = TRUE)
pvalue.alpha.evenness <- wilcox.test(Urine, Vagina, paired = TRUE)$p.value
pvalue.alpha.evenness 

Urine <- subset(sample_list,  Site == "Urine", alpha.ovserved_OTU,
                drop = TRUE)
Vagina <- subset(sample_list,  Site == "Vagina", alpha.ovserved_OTU,
                 drop = TRUE)
pvalue.alpha.ovserved_OTU <- wilcox.test(Urine, Vagina, paired = TRUE)$p.value
pvalue.alpha.ovserved_OTU


pvalue = rbind(pvalue,pvalue.alpha.shannon,pvalue.alpha.evenness,pvalue.alpha.ovserved_OTU)

##### beta diversity ##################################################
# Bray_Curtis distance
Bray_Curtis <- as.matrix(vegdist(reads_table_rare, METHODS="bray", binary=FALSE))   # see 'vegan' package for diversity analysis
Bray_Curtis <- as.data.frame(Bray_Curtis)

# Bray_Curtis distance between paired samples
x = nrow(sample_list)/2
pair_dis = as.data.frame(matrix(data = NA, nrow = x, ncol =3))
pair_dis$V1 = sample_list[(x+1):(x*2),1]
pair_dis$V2 = sample_list[(x+1):(x*2),2]

for (a in 1: x) {
   bray_x = which(row.names(Bray_Curtis) == as.character(pair_dis[a,1]))
   bray_y = which(colnames(Bray_Curtis) == as.character(pair_dis[a,2]))
   pair_dis[a,3] = Bray_Curtis[bray_x,bray_y]
}
pair_dis <- pair_dis[order(pair_dis$V3),]
colnames(pair_dis) <- c('Sample_name_u', 'Sample_name_v','Bray_Curtis')

pair_dis$Sample_name_u <- factor(pair_dis$Sample_name_u, levels = pair_dis$Sample_name_u)  # convert to factor to retain sorted order in plot.


ggplot(data=pair_dis, aes(x=Sample_name_u, y=Bray_Curtis)) +
   geom_bar(stat="identity") +
   labs(x = 'Paired samples', y = "Bray-Curtis distance")+ 
   theme(axis.title = element_text(size = 12), 
         axis.text.x=element_blank(),
         axis.text.y = element_text(size = 9),
         legend.text = element_text(size = 12), 
         legend.title = element_text(size = 12))   
ggsave("paired_distance.pdf", width=3.5, height=5.5)


# The Bray-Curtis distance of a vaginal/urine microbiome to its paired urine/vaginal microbiome was compared to the distance of this vaginal/urine microbiome to vaginal/urine microbiomes of other participants and the significance was calculated using the Wilcoxon test.
vagina_bray_distantc = Bray_Curtis[(1:(nrow(Bray_Curtis)/2)),(1:(ncol(Bray_Curtis)/2))]
vagina_bray_distantc <- gather(vagina_bray_distantc)

distance_comparison <- matrix(data = NA, nrow =(nrow(pair_dis)+nrow(vagina_bray_distantc)), ncol =2)
distance_comparison[1:nrow(pair_dis),1] = pair_dis$Bray_Curtis
distance_comparison[1:nrow(pair_dis),2] = rep('Vagina_Paired urine',nrow(pair_dis))
distance_comparison[(nrow(pair_dis)+1):nrow(distance_comparison),1] = vagina_bray_distantc$value
distance_comparison[(nrow(pair_dis)+1):nrow(distance_comparison),2] = rep('Vagina_Other vagina',nrow(vagina_bray_distantc))
distance_comparison<- as.data.frame(distance_comparison)
colnames(distance_comparison) = c('Bray-Curtis Distance', 'Type')
distance_comparison$`Bray-Curtis Distance` = as.numeric(as.character(distance_comparison$`Bray-Curtis Distance`))
distance_comparison$Type <- as.character(distance_comparison$Type)

ggplot(distance_comparison, aes(x=Type, y=`Bray-Curtis Distance`,fill = Type)) +
   geom_boxplot(outlier.shape = NA) +
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16),
         axis.text.x = element_text(angle = 60, hjust = 1))+ theme_bw()
ggsave('Bray_curtis_distance_compare.pdf', width=5, height=5)

pvalue.wil1 = wilcox.test(`Bray-Curtis Distance` ~ Type, data=distance_comparison)$p.value 
pvalue.Bray_Curtis_comparision = pvalue.wil1
pvalue = rbind(pvalue,pvalue.wil1 )


urine_bray_distantc = Bray_Curtis[(((nrow(Bray_Curtis)/2)+1):nrow(Bray_Curtis)),(((ncol(Bray_Curtis)/2)+1):ncol(Bray_Curtis))]
urine_bray_distantc <- gather(urine_bray_distantc)

distance_comparison <- matrix(data = NA, nrow =(nrow(pair_dis)+nrow(urine_bray_distantc)), ncol =2)
distance_comparison[1:nrow(pair_dis),1] = pair_dis$Bray_Curtis
distance_comparison[1:nrow(pair_dis),2] = rep('Urine_Paired vagina',nrow(pair_dis))
distance_comparison[(nrow(pair_dis)+1):nrow(distance_comparison),1] = urine_bray_distantc$value
distance_comparison[(nrow(pair_dis)+1):nrow(distance_comparison),2] = rep('Urine_Other urine',nrow(urine_bray_distantc))
distance_comparison<- as.data.frame(distance_comparison)
colnames(distance_comparison) = c('Bray-Curtis Distance', 'Type')
distance_comparison$`Bray-Curtis Distance` = as.numeric(as.character(distance_comparison$`Bray-Curtis Distance`))
distance_comparison$Type <- as.character(distance_comparison$Type)

ggplot(distance_comparison, aes(x=Type, y=`Bray-Curtis Distance`,fill = Type)) +
   geom_boxplot(outlier.shape = NA) +
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16),
         axis.text.x = element_text(angle = 60, hjust = 1))+ theme_bw()
ggsave('Bray_curtis_distance_compare_urine.pdf', width=5, height=5)

pvalue.wil2 = wilcox.test(`Bray-Curtis Distance` ~ Type, data=distance_comparison)$p.value 
pvalue.Bray_Curtis_comparision = pvalue.wil2
pvalue = rbind(pvalue,pvalue.wil2 )

# Running Nonmetric Multidimensional Scaling (NMDS) Ordination   # see 'vegan' package for diversity analysis
NMDS <-
   metaMDS(Bray_Curtis,
           distance = "bray",
           k = 2,
           maxit = 999, 
           trymax = 50,
           wascores = TRUE)

mds_data <- as.data.frame(NMDS$points)
mds_data$Site <- sample_list$Site
mds_data$Pair <- sample_list$Pair
ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Site)) +
   geom_point()
ggsave("NMDS1.pdf", width=5.5, height=5.5)

ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Site)) +
   geom_point()+
   geom_line(aes(group = Pair), size = 0.3, alpha = 0.7, color = 'black')
ggsave("NMDS2.pdf", width=5.5, height=5.5)

# bar plot for beta diversity
vagina2vagina = as.matrix(Bray_Curtis[1:(length(Bray_Curtis)/2),1:(length(Bray_Curtis)/2)])
vagina2vagina = as.vector(vagina2vagina)
vagina2vagina = vagina2vagina[vagina2vagina!=0]
vagina = as.data.frame(cbind(vagina2vagina,'Vagina'))
colnames(vagina) = c('Bray_Curtis_distance','Body_site')

Urine2vagina = as.matrix(Bray_Curtis[1:(length(Bray_Curtis)/2),(length(Bray_Curtis)/2+1):length(Bray_Curtis)])
Urine2vagina = as.vector(Urine2vagina)
urine = as.data.frame(cbind(Urine2vagina,'Urine'))
colnames(urine) = c('Bray_Curtis_distance','Body_site')

beta_plot = rbind(vagina,urine)
beta_plot$Body_site <- as.factor(beta_plot$Body_site)
beta_plot$Bray_Curtis_distance <- as.numeric(as.character(beta_plot$Bray_Curtis_distance))

ggplot(beta_plot, aes(x=Body_site, y=Bray_Curtis_distance,fill=Body_site)) +
   geom_boxplot()  +
   labs(x = NULL, y = "Distance to Vagina samples")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16))

ggsave("Bray_Curtis_bar.pdf", width=4, height=5.5)

# Compute perMANOVA test for the significance of beta diversity
# adonis2 : permutational multivariate analysis of variance using distance matrices the centroid and/or the spread of the objects is different between the groups.
pvalue.adonis.site <- adonis2(reads_table_rare ~ Site, data = sample_list, method = "bray")  
pvalue.adonis.site <- pvalue.adonis.site[1,5]
pvalue <- rbind(pvalue,pvalue.adonis.site)
pvalue.adonis.pair <- adonis2(reads_table_rare ~ Pair, data = sample_list, method = "bray")  
pvalue.adonis.pair <- pvalue.adonis.pair[1,5]
pvalue <- rbind(pvalue,pvalue.adonis.pair)

write.csv(pvalue,'pvalue.csv')

# the same NMDS figure color coded by vagitype 
# Samples were assigned to vagitypes according to the species with an relative abundance of 30% or greater
vagitype_v = matrix(data = NA, ncol =1, nrow = ncol(reads_table_v_new))
colnames(vagitype_v) = 'Vagitype'

for (a in 1:ncol(reads_table_v_new)) {
   n= which(reads_table_v_new[,a] == max(reads_table_v_new[,a]), arr.ind=TRUE)
   vagitype_v[a,1]= row.names(reads_table_v_new)[n]
}

vagitype_u= matrix(data = NA, ncol =1, nrow = ncol(reads_table_u_new))
colnames(vagitype_u) = 'Vagitype'

for (a in 1:ncol(reads_table_u_new)) {
   n= which(reads_table_u_new[,a] == max(reads_table_u_new[,a]), arr.ind=TRUE)
   vagitype_u[a,1]= row.names(reads_table_u_new)[n]
}

vagitype = rbind(vagitype_v,vagitype_u)
vagitype_unique = unique(vagitype)
for (a in 1:length(vagitype_unique)) {
   if (sum(vagitype == vagitype_unique[a]) < 5) {
      vagitype[vagitype == vagitype_unique[a]] = 'Others'
   }
}

mds_data = cbind(mds_data,vagitype)

mds_data$Vagitype <- as.factor(mds_data$Vagitype)
mds_data$Vagitype <- factor(mds_data$Vagitype, levels = c('Lactobacillus iners ','Enterobacteriaceae cluster31 ',
                                                          'Lachnospiraceae BVAB1 ','Gardnerella vaginalis ',
                                                          'Atopobium vaginae ','Lactobacillus crispatus_cluster ',
                                                          'Lactobacillus gasseri_cluster ','Others'))


ggplot(mds_data, aes(MDS1, MDS2))  +
   geom_line(aes(group = Pair), size = 0.3, alpha = 0.7) +
   geom_point(size=3,aes(shape = Site, color = Vagitype)) +
   scale_color_manual(values=c('lightblue','deepskyblue','orange','red','brown','yellow','pink','grey')) +
   coord_fixed()+ 
   theme(
      axis.title.x = element_text( size=16),
      axis.title.y = element_text( size=16),
      legend.text = element_text(size=16),
      legend.title = element_text(size=16),
      plot.title = element_text(hjust = 0.5, size = 14)
   ) + theme_bw()
ggsave("vagitype.pdf", width=7, height=4)

same_vagitype = matrix(data = NA, nrow= nrow(mds_data)/2, ncol =2)
same_vagitype[,1] = mds_data$Vagitype[1:(nrow(mds_data)/2)]
same_vagitype[,2] = mds_data$Vagitype[(nrow(mds_data)/2+1):nrow(mds_data)]
number_vagitype_2 = sum(same_vagitype[,1] == same_vagitype[,2])
pvalue <- rbind(pvalue,number_vagitype_2)

write.csv(pvalue,'pvalue.csv')



##### species with differential abundance #############################
reads_table_das <- as.data.frame((reads_table)) # see 'ALDEx2' package for details from line 594 to 610

conds <- sample_list$Site

x <- aldex.clr(reads_table_das, conds, mc.samples=128, denom="all", verbose=F)

# paired Wilcoxon Rank Sum test and Welch's t-test
x.tt <- aldex.ttest(x, paired.test=TRUE)

x.effect <- aldex.effect(x)

x.all <- data.frame(cbind(x.tt,x.effect))

abundance_change_das <- x.all$diff.btw # diff.btw is a vector containing the per-feature median difference between condition A and B


x.all <- cbind(abundance_change_das,x.all)


# Reform the 'ALDEx2' results for plotting
`-Log10(adj-pvalue)` <- -log10(x.all$wi.eBH)     # use wi.eBH as the adj-pvalue
x.all$abundance_das <- abundance_das
x.all$`-Log10(adj-pvalue)` <- `-Log10(adj-pvalue)`

# draw figure
das <- x.all[(x.all$`-Log10(adj-pvalue)` >=1.301 & (x.all$abundance_change_das >=1 | x.all$abundance_change_das <=-1)),]

das$Species <- row.names(das)
das <- das[order(das$abundance_change_das),]    
das$Color <- ifelse(das$abundance_change_das < 0, "Enriched in Urine", "Enriched in Vagina")  # above / below avg flag
das$Species <- factor(das$Species, levels = das$Species)  # convert to factor to retain sorted order in plot.

das$abundance_change_das[das$abundance_change_das == Inf] = 5
das$abundance_change_das[das$abundance_change_das == -Inf] = -5
das$abundance_change_das[das$abundance_change_das <= -5] = -5
das$abundance_change_das[das$abundance_change_das >= 5] = 5

theme_set(theme_bw())  

ggplot(das, aes(Species, abundance_change_das)) + 
   geom_point(aes(col=Color, size=`-Log10(adj-pvalue)`)) + 
   coord_flip() +          # convert x y axis
   labs(x = 'Species', y = "Median difference in clr values")+ 
   theme(axis.title = element_text(size = 12), 
         axis.text = element_text(size = 12), 
         legend.text = element_text(size = 12), 
         legend.title = element_text(size = 12))    
ggsave("das.pdf", width=8, height=3.5)




##### draw bar figure ########################

# find top n species
reads_table_v_new_bar <- reads_table_v_new

reads_table_v_new_bar$total_reads <- rowSums(reads_table_v_new_bar)
reads_table_v_new_bar = reads_table_v_new_bar[order(reads_table_v_new_bar$total_reads,decreasing = T),]
reads_table_v_new_bar$total_reads <- NULL
other_reads = as.data.frame(t(colSums(reads_table_v_new_bar[10:nrow(reads_table_v_new_bar),])))
reads_table_v_new_bar <- reads_table_v_new_bar[-c(10:nrow(reads_table_v_new_bar)), ]
reads_table_v_new_bar <- rbind(reads_table_v_new_bar,other_reads)
row.names(reads_table_v_new_bar)[10] = 'Others'

reads_table_v_new_bar_abundance <- matrix(data =0, ncol = ncol(reads_table_v_new_bar),nrow = nrow(reads_table_v_new_bar))

for (a in 1:ncol(reads_table_v_new)) {
   reads_table_v_new_bar_abundance[,a] <- reads_table_v_new_bar[,a] / colSums(reads_table_v_new_bar)[a]
   
}
row.names(reads_table_v_new_bar_abundance) = row.names(reads_table_v_new_bar)
colnames(reads_table_v_new_bar_abundance) = colnames(reads_table_v_new_bar) 

reads_table_v_new_bar <- as.data.frame(reads_table_v_new_bar_abundance)
reads_table_v_new_bar <- as.data.frame(t(reads_table_v_new_bar))
reads_table_v_new_bar = reads_table_v_new_bar[order(reads_table_v_new_bar$`Lactobacillus iners`,decreasing = T),]
reads_table_v_new_bar <- as.data.frame(t(reads_table_v_new_bar))
reads_table_v_new_bar <- reads_table_v_new_bar/(colSums(reads_table_v_new_bar))

plot_v <- gather(reads_table_v_new_bar)
plot_v_name <- row.names(reads_table_v_new_bar)
plot_v_name <- rep( plot_v_name , ncol(reads_table_v_new))
plot_v$Species <- plot_v_name

plot_v$Species <- as.factor(plot_v$Species)
plot_v$Species <- factor(plot_v$Species, levels = plot_v$Species[1:10])

plot_v$key <- as.factor(plot_v$key)
plot_v$key <- factor(plot_v$key, levels = colnames(reads_table_v_new_bar))

ggplot(data=plot_v, aes(x=key, y=value, fill=Species)) +
   geom_bar(stat="identity") +
   scale_fill_manual(values=c('#A8C5D3','#FC9800','#C80B0B','#FDFFBA','#55ff42','#F7A0A0','#D7A7F1',
                              '#3264B8','#bfbfbf','#3F3F3F')) +
   labs(x = 'Vaginal samples', y = "Abundance")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank())

ggsave("bar_v.pdf", width=9, height=3.5)


# urine
# find top n species
reads_table_u_new_bar <- reads_table_u_new

reads_table_u_new_bar$total_reads <- rowSums(reads_table_u_new_bar)
reads_table_u_new_bar = reads_table_u_new_bar[order(reads_table_u_new_bar$total_reads,decreasing = T),]
reads_table_u_new_bar$total_reads <- NULL
other_reads = as.data.frame(t(colSums(reads_table_u_new_bar[10:nrow(reads_table_u_new_bar),])))
reads_table_u_new_bar <- reads_table_u_new_bar[-c(10:nrow(reads_table_u_new_bar)), ]
reads_table_u_new_bar <- rbind(reads_table_u_new_bar,other_reads)
row.names(reads_table_u_new_bar)[10] = 'Others'

reads_table_u_new_bar_abundance <- matrix(data =0, ncol = ncol(reads_table_u_new_bar),nrow = nrow(reads_table_u_new_bar))

for (a in 1:ncol(reads_table_u_new)) {
   reads_table_u_new_bar_abundance[,a] <- reads_table_u_new_bar[,a] / colSums(reads_table_u_new_bar)[a]
   
}
row.names(reads_table_u_new_bar_abundance) = row.names(reads_table_u_new_bar)
colnames(reads_table_u_new_bar_abundance) = colnames(reads_table_u_new_bar) 

reads_table_u_new_bar <- as.data.frame(reads_table_u_new_bar_abundance)
reads_table_u_new_bar <- as.data.frame(t(reads_table_u_new_bar))
reads_table_u_new_bar = reads_table_u_new_bar[order(reads_table_u_new_bar$`Lactobacillus iners`,decreasing = T),]
reads_table_u_new_bar <- as.data.frame(t(reads_table_u_new_bar))
reads_table_u_new_bar <- reads_table_u_new_bar/(colSums(reads_table_u_new_bar))


# urine with marched order
plot_u <- gather(reads_table_u_new_bar)
plot_u_name <- row.names(reads_table_u_new_bar)
plot_u_name <- rep( plot_u_name , ncol(reads_table_u_new))
plot_u$Species <- plot_u_name

plot_u$Species <- as.factor(plot_u$Species)
plot_u$Species <- factor(plot_u$Species, levels = plot_u$Species[1:10])

plot_u$key <- as.factor(plot_u$key)
order <- colnames(reads_table_v_new_bar)
order_urine <- vector(mode = "character", length = length(order))
sample_list_2 <- sample_list[(length(order)+1):nrow(sample_list),]
for (a in 1: length(order)) {
   n = which(sample_list_2[,2] == order[a])
   order_urine[a] = sample_list_2[n,1]
}
plot_u$key <- factor(plot_u$key, levels = order_urine)

ggplot(data=plot_u, aes(x=key, y=value, fill=Species)) +
   geom_bar(stat="identity") +
   scale_fill_manual(values=c('#A8C5D3','#3be5f7','#FC9800','#C80B0B','#D7A7F1','#55ff42','#FDFFBA',
                              '#3264B8','#F7A0A0','#3F3F3F')) +
   labs(x = 'Urine samples', y = "Abundance")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank())
ggsave("bar_u_vagina_order.pdf", width=9, height=3.5)


##### create new abundance table ##################
reads_table_u_abundance <- matrix(data =0, ncol = ncol(reads_table_u_new),nrow = nrow(reads_table_u_new))

for (a in 1:ncol(reads_table_u_new)) {
   reads_table_u_abundance[,a] <- reads_table_u_new[,a] / colSums(reads_table_u_new)[a]
}
reads_table_u_abundance <- as.data.frame(reads_table_u_abundance)
row.names(reads_table_u_abundance) = row.names(reads_table_u_new)
colnames(reads_table_u_abundance) = colnames(reads_table_v_new)

reads_table_v_abundance <- matrix(data =0, ncol = ncol(reads_table_v_new),nrow = nrow(reads_table_v_new))

for (a in 1:ncol(reads_table_v_new)) {
   reads_table_v_abundance[,a] <- reads_table_v_new[,a] / colSums(reads_table_v_new)[a]
}
reads_table_v_abundance <- as.data.frame(reads_table_v_abundance)
row.names(reads_table_v_abundance) = row.names(reads_table_v_new)
colnames(reads_table_v_abundance) = colnames(reads_table_v_new)

total_abundance_v = as.data.frame(rowSums(reads_table_v_abundance))
total_abundance_v <- total_abundance_v/sum(total_abundance_v$`rowSums(reads_table_v_abundance)`)

total_abundance_u = as.data.frame(rowSums(reads_table_u_abundance))
total_abundance_u <- total_abundance_u/sum(total_abundance_u$`rowSums(reads_table_u_abundance)`)


##### linear regression  ################

line_list <- matrix(data = NA, nrow = nrow(reads_table_v_abundance), ncol =7)
line_list <- as.data.frame(line_list)
colnames(line_list) = c('Species','p-value','R-value',
                        'linear_model_function','p-value_remove_low_abundance',
                        'R-value_remove_low_abundance','linear_model_function_remove_low_abundance')
row.names(line_list)= row.names(reads_table_v_abundance)
line_list$Species = row.names(reads_table_v_abundance)

for (a in 1:nrow(line_list)) {
   # no treatment
   line <- as.data.frame(t(rbind(reads_table_u_abundance[a,],reads_table_v_abundance[a,])))
   
   line <- log10(line)
   
   colnames(line) = c('Urine','Vagina')
   
   # linear regression
   line2=line
   line2[line2 == '-Inf'] = -6  # all zeros in the reads table are replaced by 10E-6 for linear regression analysis
   linearMod <- lm(Urine ~ Vagina, data=line2)  #  function 'lm' is for linear regression
   
   b = format(round(linearMod$coefficients[2], 3), nsmall = 3)  
   c = format(round(linearMod$coefficients[1], 3), nsmall = 3)  
   
   if (linearMod$coefficients[1] > 0){  # output the formula of linear regression
      line_list[a,4]= paste0('y=',b,'x+',c)
   } else if (c < 0) {
      line_list[a,4]= paste0('y=',b,'x',c)
   } else {
      line_list[a,4]= paste0('y=',b,'x')
   }
   
   r_squared <- summary(linearMod)$r.squared  # get R-value
   line_list[a,3]= r_squared
   pvalue <- as.numeric(summary(linearMod)$coefficients[,4][2])  # get p-value
   line_list[a,2]= pvalue
   
   # linear regression after removing low abundance (just a try, results not present in the publication)
   keep = line2$Urine > -3 | line2$Vagina > -3 
   line3 = line2[keep,]
   
   if (sum(keep) >=5) {
      linearMod <- lm(Urine ~ Vagina, data=line3) 
      
      b = format(round(linearMod$coefficients[2], 3), nsmall = 3)
      c = format(round(linearMod$coefficients[1], 3), nsmall = 3)
      
      if (linearMod$coefficients[1] > 0){
         line_list[a,7]= paste0('y=',b,'x+',c)
      } else if (c < 0) {
         line_list[a,7]= paste0('y=',b,'x',c)
      } else {
         line_list[a,7]= paste0('y=',b,'x')
      }
      
      r_squared2 <- summary(linearMod)$r.squared
      line_list[a,6]= r_squared2
      pvalue2 <- as.numeric(summary(linearMod)$coefficients[,4][2])
      line_list[a,5]= pvalue2
   }
   
   
   # regression plotting of four PTB-associated taxa
   if (line_list[a,1] %in% c('Sneathia amnii ','Lachnospiraceae BVAB1 ','TM7 OTU-H1 ','Prevotella cluster2 ')) {
      ggplot(line, aes(Vagina, Urine) )  +
         geom_point(size=1) +
         xlab(paste0("Urine abundance")) +
         ylab(paste0("Vagina abundance")) + ggtitle(paste0(line_list[a,1], '\n',"R = ", r_squared, '\n',"pvalue = ", pvalue))  + 
         geom_smooth(method='lm', se=T) + 
         theme(axis.title.x = element_text( size=16),
               axis.title.y = element_text( size=16)) 
      ggsave(paste0('linear_1_',line_list[a,1],'.pdf'), width=5.5, height=5.5)
   }
   
}


line_list$adj_p = p.adjust(line_list$`p-value`)  #  Benjamini & Hochberg FDR (false discovery rate)
line_list$adj_p_remove = p.adjust(line_list$`p-value_remove_low_abundance`)  #  Benjamini & Hochberg 


line_list = cbind(line_list,total_abundance_v, total_abundance_u)
colnames(line_list)[10] = 'Abundance in the vaginal microbiome' 
colnames(line_list)[11] = 'Abundance in the urinary microbiome' 

write.csv(line_list,'linear_regression.csv')



###################### threshold to remove urine-specific taxa ###########
reads_table_u_new2 <- reads_table_u_new
reads_table_v_new2 <- reads_table_v_new

# get six Enterobacteriaceae
reads_table_u_abundance <- matrix(data =0, ncol = ncol(reads_table_u_new),nrow = nrow(reads_table_u_new))

for (a in 1:ncol(reads_table_u_new)) {
   reads_table_u_abundance[,a] <- reads_table_u_new[,a] / colSums(reads_table_u_new)[a]
}
reads_table_u_abundance <- as.data.frame(reads_table_u_abundance)
row.names(reads_table_u_abundance) = row.names(reads_table_u_new)
colnames(reads_table_u_abundance) = colnames(reads_table_v_new)

keep = row.names(reads_table_u_abundance) == 'Enterobacteriaceae cluster31 '
a = reads_table_u_abundance[keep,]
keep = a[1,] >= 0.8
a = a[,keep]

#if ratio of total abundance is larger than 10 folder, the species is removed 
abundance_ratio <- matrix(data = NA, nrow = nrow(reads_table_u_new), ncol =3)

abundance_ratio[,1] <- row.names(reads_table_u_new)
abundance_ratio[,2] <- rowSums(reads_table_u_abundance)/ rowSums(reads_table_v_abundance)
abundance_ratio[,3] <- rowSums(reads_table_u_abundance) +rowSums(reads_table_v_abundance)

abundance_ratio <- as.data.frame(abundance_ratio)
abundance_ratio$V2 <- as.numeric(as.character(abundance_ratio$V2))
abundance_ratio$V3 <- as.numeric(as.character(abundance_ratio$V3))

keep = abundance_ratio$V2 < 10

# Remove species present (0.01% threshold) urine / vagina >= 5 
present_urine <- present == 'Urine' | present == 'Both'
present_urine <- colSums(present_urine)
present_vagina <- present == 'Vagina' | present == 'Both'
present_vagina <- colSums(present_vagina)
present_ratio <- present_urine/present_vagina
present_ratio <- as.data.frame(present_ratio)

keep2 = present_ratio$present_ratio < 2.5 

#keep3 = keep & keep2
keep3 = keep2

reads_table_v_new2 <- reads_table_v_new2[keep3,]
reads_table_u_new2 <- reads_table_u_new2[keep3,]

########## new reads and abundance table ############
# sample total reads threshold
keep <- matrix(ncol = dim(reads_table_v_new2)[2],)
for (a in 1: dim(reads_table_v_new2)[2]) {
   keep[,a] = sum(reads_table_v_new2[,a]) >=10000    # input
}
reads_table_u_new2 <- reads_table_u_new2[,keep]
reads_table_v_new2 <- reads_table_v_new2[,keep]

keep <- matrix(ncol = dim(reads_table_u_new2)[2],)
for (a in 1: dim(reads_table_u_new2)[2]) {
   keep[,a] = sum(reads_table_u_new2[,a]) >=10000    # input
}

reads_table_u_new2 <- reads_table_u_new2[,keep]
reads_table_v_new2 <- reads_table_v_new2[,keep]

reads_table = cbind(reads_table_v_new2, reads_table_u_new2)

# get abundance table
reads_table_abundance <- matrix(data =0, ncol = ncol(reads_table),nrow = nrow(reads_table))

for (a in 1:ncol(reads_table)) {
   reads_table_abundance[,a] <- reads_table[,a] / colSums(reads_table)[a]
}
reads_table_abundance <- as.data.frame(reads_table_abundance)
row.names(reads_table_abundance) = row.names(reads_table)
colnames(reads_table_abundance) = colnames(reads_table)

reads_table_v_abundance <- matrix(data =0, ncol = ncol(reads_table_v_new2),nrow = nrow(reads_table_v_new2))

for (a in 1:ncol(reads_table_v_new2)) {
   reads_table_v_abundance[,a] <- reads_table_v_new2[,a] / colSums(reads_table_v_new2)[a]
}
reads_table_v_abundance <- as.data.frame(reads_table_v_abundance)
row.names(reads_table_v_abundance) = row.names(reads_table_v_new2)
colnames(reads_table_v_abundance) = colnames(reads_table_v_new2)

reads_table_u_abundance <- matrix(data =0, ncol = ncol(reads_table_u_new2),nrow = nrow(reads_table_u_new2))

for (a in 1:ncol(reads_table_u_new2)) {
   reads_table_u_abundance[,a] <- reads_table_u_new2[,a] / colSums(reads_table_u_new2)[a]
}
reads_table_u_abundance <- as.data.frame(reads_table_u_abundance)
row.names(reads_table_u_abundance) = row.names(reads_table_u_new2)
colnames(reads_table_u_abundance) = colnames(reads_table_v_new2)

total_abundance_v = as.data.frame(rowSums(reads_table_v_abundance))
total_abundance_v <- total_abundance_v/sum(total_abundance_v$`rowSums(reads_table_v_abundance)`)


total_abundance_u = as.data.frame(rowSums(reads_table_u_abundance))
total_abundance_u <- total_abundance_u/sum(total_abundance_u$`rowSums(reads_table_u_abundance)`)

# get samples metadata
sample_list_v = as.matrix(colnames(reads_table_v_new2))
sample_list_u = as.matrix(colnames(reads_table_u_new2))

sample_list = as.data.frame(matrix(data = NA, ncol=3, nrow= 2*nrow(sample_list_u)))
colnames(sample_list) = c('Sample','Pair','Site')
sample_list$Sample = rbind(sample_list_v,sample_list_u)
sample_list$Pair = rbind(sample_list_v,sample_list_v)
sample_list$Site = c(matrix(data = 'Vagina' , nrow = nrow(sample_list_u)),matrix(data = 'Urine' , nrow = nrow(sample_list_u)))

reads_table_u_new <- reads_table_u_new2
reads_table_v_new <- reads_table_v_new2
##### rarefaction ##################################################

reads_table_rare <- t(reads_table)
reads_table_rare = Rarefy(reads_table_rare, depth = min(rowSums(reads_table_rare)))
reads_table_rare <- reads_table_rare$otu.tab.rff
reads_table_rare <- as.data.frame(reads_table_rare)



##### alpha diversity ##################################################

alpha.shannon_diversity <- data.frame(diversity(reads_table_rare))
alpha.evenness <- alpha.shannon_diversity/log(specnumber(reads_table_rare))
alpha.ovserved_OTU <- data.frame(colSums(t(reads_table_rare) != 0))

sample_list$alpha.shannon <- alpha.shannon_diversity$diversity.reads_table_rare.
sample_list$alpha.evenness <- alpha.evenness$diversity.reads_table_rare.
sample_list$alpha.ovserved_OTU <- alpha.ovserved_OTU$colSums.t.reads_table_rare.....0.
#colnames(sample_list)=c('Sample','Pair','Site','Shannon','Evenness','Observed OTU')

# creat images for alpha diversity

geom_boxplot(outlier.colour="black", outlier.shape=16,
             outlier.size=2, notch=FALSE)

ggplot(sample_list, aes(x=Site, y=alpha.shannon ,color = Site)) +
   geom_boxplot() +
   geom_jitter(shape=16, position=position_jitter(0.2))+
   labs(x = NULL, y = "Shannon index")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16))+ theme_bw()
ggsave('shannon__2.pdf', width=3.5, height=3.5)

ggplot(sample_list, aes(x=Site, y=alpha.shannon)) +
   geom_line(aes(group = Pair), size = 0.1, alpha = 0.7) +
   geom_point(aes(color = Site), alpha = 0.5) +
   labs(x = NULL, y = "Shannon index")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16))+ theme_bw()
ggsave("shannon_2__2.pdf", width=5.5, height=5.5)

ggplot(sample_list, aes(x=Site, y=alpha.evenness ,color = Site)) +
   geom_boxplot() +
   geom_jitter(shape=16, position=position_jitter(0.2))+
   labs(x = NULL, y = "Evenness")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16))+ theme_bw()
ggsave('evenness__2.pdf', width=3.5, height=3.5)

ggplot(sample_list, aes(x=Site, y=alpha.evenness)) +
   geom_line(aes(group = Pair), size = 0.1, alpha = 0.7) +
   geom_point(aes(color = Site), alpha = 0.5) +
   labs(x = NULL, y = "Evenness")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16))+ theme_bw()
ggsave("evenness_2__2.pdf", width=5.5, height=5.5)

ggplot(sample_list, aes(x=Site, y=alpha.ovserved_OTU ,color = Site)) +
   geom_boxplot() +
   geom_jitter(shape=16, position=position_jitter(0.2))+
   labs(x = NULL, y = "Observed Species")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16))+ theme_bw()
ggsave('otu__2.pdf', width=3.5, height=3.5)

ggplot(sample_list, aes(x=Site, y=alpha.ovserved_OTU)) +
   geom_line(aes(group = Pair), size = 0.1, alpha = 0.7) +
   geom_point(aes(color = Site), alpha = 0.5) +
   labs(x = NULL, y = "Observed species")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16))+ theme_bw()
ggsave("otu_2__2.pdf", width=5.5, height=5.5)

# Compute paired Wilcoxon test for the significance of alpha diversity
Urine <- subset(sample_list,  Site == "Urine", alpha.shannon,
                drop = TRUE)
Vagina <- subset(sample_list,  Site == "Vagina", alpha.shannon,
                 drop = TRUE)
pvalue.alpha.shannon <- wilcox.test(Urine, Vagina, paired = TRUE)$p.value
pvalue.alpha.shannon

Urine <- subset(sample_list,  Site == "Urine", alpha.evenness,
                drop = TRUE)
Vagina <- subset(sample_list,  Site == "Vagina", alpha.evenness,
                 drop = TRUE)
pvalue.alpha.evenness <- wilcox.test(Urine, Vagina, paired = TRUE)$p.value
pvalue.alpha.evenness 

Urine <- subset(sample_list,  Site == "Urine", alpha.ovserved_OTU,
                drop = TRUE)
Vagina <- subset(sample_list,  Site == "Vagina", alpha.ovserved_OTU,
                 drop = TRUE)
pvalue.alpha.ovserved_OTU <- wilcox.test(Urine, Vagina, paired = TRUE)$p.value
pvalue.alpha.ovserved_OTU


pvalue = rbind(pvalue.alpha.shannon,pvalue.alpha.evenness,pvalue.alpha.ovserved_OTU)

##### beta diversity ##################################################
# Bray_Curtis distance
Bray_Curtis <- as.matrix(vegdist(reads_table_rare, METHODS="bray", binary=FALSE))
Bray_Curtis <- as.data.frame(Bray_Curtis)

# Bray_Curtis distance between paired samples
x = nrow(sample_list)/2
pair_dis = as.data.frame(matrix(data = NA, nrow = x, ncol =3))
pair_dis$V1 = sample_list[(x+1):(x*2),1]
pair_dis$V2 = sample_list[(x+1):(x*2),2]

for (a in 1: x) {
   bray_x = which(row.names(Bray_Curtis) == as.character(pair_dis[a,1]))
   bray_y = which(colnames(Bray_Curtis) == as.character(pair_dis[a,2]))
   pair_dis[a,3] = Bray_Curtis[bray_x,bray_y]
}
pair_dis <- pair_dis[order(pair_dis$V3),]
colnames(pair_dis) <- c('Sample_name_u', 'Sample_name_v','Bray_Curtis')

pair_dis$Sample_name_u <- factor(pair_dis$Sample_name_u, levels = pair_dis$Sample_name_u)  # convert to factor to retain sorted order in plot.


ggplot(data=pair_dis, aes(x=Sample_name_u, y=Bray_Curtis)) +
   geom_bar(stat="identity") +
   labs(x = 'Paired samples', y = "Bray-Curtis distance between paired samples")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text.x=element_blank(),
         axis.text.y = element_text(size = 9),
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16))  + theme_bw() 

ggsave("paired_distance_2.pdf", width=3.5, height=5.5)


# distance comparison 
vagina_bray_distantc = Bray_Curtis[(1:(nrow(Bray_Curtis)/2)),(1:(ncol(Bray_Curtis)/2))]
vagina_bray_distantc <- gather(vagina_bray_distantc)

distance_comparison_2 <- matrix(data = NA, nrow =(nrow(pair_dis)+nrow(vagina_bray_distantc)), ncol =2)
distance_comparison_2[1:nrow(pair_dis),1] = pair_dis$Bray_Curtis
distance_comparison_2[1:nrow(pair_dis),2] = rep('Vagina_Paired urine',nrow(pair_dis))
distance_comparison_2[(nrow(pair_dis)+1):nrow(distance_comparison_2),1] = vagina_bray_distantc$value
distance_comparison_2[(nrow(pair_dis)+1):nrow(distance_comparison_2),2] = rep('Vagina_Other vagina',nrow(vagina_bray_distantc))
distance_comparison_2<- as.data.frame(distance_comparison_2)
colnames(distance_comparison_2) = c('Bray-Curtis Distance', 'Type')
distance_comparison_2$`Bray-Curtis Distance` = as.numeric(as.character(distance_comparison_2$`Bray-Curtis Distance`))
distance_comparison_2$Type <- as.character(distance_comparison_2$Type)

ggplot(distance_comparison_2, aes(x=Type, y=`Bray-Curtis Distance`,fill = Type)) +
   geom_boxplot() +
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16),
         axis.text.x = element_text(angle = 60, hjust = 1))+ theme_bw()
ggsave('Bray_curtis_distance_compare_2.pdf', width=5, height=5)

pvalue.wil1 = wilcox.test(`Bray-Curtis Distance` ~ Type, data=distance_comparison_2)$p.value 
pvalue.Bray_Curtis_comparision = pvalue.wil1
pvalue = rbind(pvalue,pvalue.wil1 )

distance_comparison$Threshold = 'N'
distance_comparison_2$Threshold = 'Y'
distance_comparison_3 = rbind(distance_comparison, distance_comparison_2)
distance_comparison_3$Type = paste(distance_comparison_3$Type, distance_comparison_3$Threshold, sep = '_')
distance_comparison_3 = distance_comparison_3[distance_comparison_3$Type == 'Vagina_Paired urine_N'
                                              | distance_comparison_3$Type == 'Vagina_Paired urine_Y',]

ggplot(distance_comparison_3, aes(x=Type, y=`Bray-Curtis Distance`,fill = Type)) +
   geom_boxplot() +
   geom_jitter(size = 1) +
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16),
         axis.text.x = element_text(angle = 60, hjust = 1))+ theme_bw()
ggsave('Bray_curtis_distance_compare_3.pdf', width=5, height=5)

pvalue.wil1 = wilcox.test(`Bray-Curtis Distance` ~ Type, data=distance_comparison_3)$p.value 
pvalue.Bray_Curtis_comparision = pvalue.wil1
pvalue = rbind(pvalue,pvalue.wil1 )

# Running Nonmetric Multidimensional Scaling (NMDS) Ordination
NMDS <-
   metaMDS(Bray_Curtis,
           distance = "bray",
           k = 2,
           maxit = 999, 
           trymax = 50,
           wascores = TRUE)

mds_data <- as.data.frame(NMDS$points)
mds_data$Site <- sample_list$Site
mds_data$Pair <- sample_list$Pair
ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Site)) +
   geom_point()+
   stat_ellipse(type = "t")
ggsave("NMDS1_2.pdf", width=5.5, height=5.5)

ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Site)) +
   geom_point()+
   geom_line(aes(group = Pair), size = 0.3, alpha = 0.7, color = 'black')
ggsave("NMDS2_2.pdf", width=5.5, height=5.5)

# bar plot for beta diversity
vagina2vagina = as.matrix(Bray_Curtis[1:(length(Bray_Curtis)/2),1:(length(Bray_Curtis)/2)])
vagina2vagina = as.vector(vagina2vagina)
vagina2vagina = vagina2vagina[vagina2vagina!=0]
vagina = as.data.frame(cbind(vagina2vagina,'Vagina'))
colnames(vagina) = c('Bray_Curtis_distance','Body_site')

Urine2vagina = as.matrix(Bray_Curtis[1:(length(Bray_Curtis)/2),(length(Bray_Curtis)/2+1):length(Bray_Curtis)])
Urine2vagina = as.vector(Urine2vagina)
urine = as.data.frame(cbind(Urine2vagina,'Urine'))
colnames(urine) = c('Bray_Curtis_distance','Body_site')

beta_plot = rbind(vagina,urine)
beta_plot$Body_site <- as.factor(beta_plot$Body_site)
beta_plot$Bray_Curtis_distance <- as.numeric(as.character(beta_plot$Bray_Curtis_distance))

ggplot(beta_plot, aes(x=Body_site, y=Bray_Curtis_distance,fill=Body_site)) +
   geom_boxplot()  +
   labs(x = NULL, y = "Distance to Vagina samples")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16))+ theme_bw()

ggsave("Bray_Curtis_bar_2.pdf", width=4, height=5.5)

# vagitype
vagitype_v = matrix(data = NA, ncol =1, nrow = ncol(reads_table_v_new))
colnames(vagitype_v) = 'Vagitype'

for (a in 1:ncol(reads_table_v_new)) {
   n= which(reads_table_v_new[,a] == max(reads_table_v_new[,a]), arr.ind=TRUE)
   vagitype_v[a,1]= row.names(reads_table_v_new)[n]
}

vagitype_u= matrix(data = NA, ncol =1, nrow = ncol(reads_table_u_new))
colnames(vagitype_u) = 'Vagitype'

for (a in 1:ncol(reads_table_u_new)) {
   n= which(reads_table_u_new[,a] == max(reads_table_u_new[,a]), arr.ind=TRUE)
   vagitype_u[a,1]= row.names(reads_table_u_new)[n]
}

vagitype = rbind(vagitype_v,vagitype_u)
vagitype_unique = unique(vagitype)
for (a in 1:length(vagitype_unique)) {
   if (sum(vagitype == vagitype_unique[a]) < 5) {
      vagitype[vagitype == vagitype_unique[a]] = 'Others'
   }
}

mds_data = cbind(mds_data,vagitype)
mds_data$Vagitype <- as.factor(mds_data$Vagitype)
mds_data$Vagitype <- factor(mds_data$Vagitype, levels = c('Lactobacillus iners ',
                                                          'Lachnospiraceae BVAB1 ','Gardnerella vaginalis ',
                                                          'Atopobium vaginae ','Lactobacillus crispatus_cluster ',
                                                          'Lactobacillus gasseri_cluster ','Others'))


ggplot(mds_data, aes(MDS1, MDS2))  +
   geom_line(aes(group = Pair), size = 0.3, alpha = 0.7) +
   geom_point(size=2,aes(shape = Site, color = Vagitype)) +
   scale_color_manual(values=c('lightblue','orange','red','brown','yellow','pink','grey')) +
   coord_fixed()+ 
   theme(
      axis.title.x = element_text( size=16),
      axis.title.y = element_text( size=16),
      legend.text = element_text(size=16),
      legend.title = element_text(size=16),
      plot.title = element_text(hjust = 0.5, size = 14)+ theme_bw()
   ) + theme_bw()
ggsave("vagitype2.pdf", width=5, height=3.5)

same_vagitype = matrix(data = NA, nrow= nrow(mds_data)/2, ncol =2)
same_vagitype[,1] = mds_data$Vagitype[1:(nrow(mds_data)/2)]
same_vagitype[,2] = mds_data$Vagitype[(nrow(mds_data)/2+1):nrow(mds_data)]
number_vagitype_2 = sum(same_vagitype[,1] == same_vagitype[,2])
pvalue <- rbind(pvalue,number_vagitype_2)



# Compute perMANOVA test for the significance of beta diversity
# adonis2 : permutational multivariate analysis of variance using distance matrices the centroid and/or the spread of the objects is different between the groups.
pvalue.adonis.site <- adonis2(reads_table_rare ~ Site, data = sample_list, method = "bray")  
pvalue.adonis.site <- pvalue.adonis.site[1,5]
pvalue <- rbind(pvalue,pvalue.adonis.site)
pvalue.adonis.pair <- adonis2(reads_table_rare ~ Pair, data = sample_list, method = "bray")  
pvalue.adonis.pair <- pvalue.adonis.pair[1,5]
pvalue <- rbind(pvalue,pvalue.adonis.pair)

write.csv(pvalue,'pvalue_2.csv')
















##### species with differential abundance #############################
reads_table_das <- as.data.frame((reads_table))

conds <- sample_list$Site

x <- aldex.clr(reads_table_das, conds, mc.samples=128, denom="all", verbose=F)

# paired Wilcoxon Rank Sum test and Welch's t-test
x.tt <- aldex.ttest(x, paired.test=TRUE)

x.effect <- aldex.effect(x)

x.all <- data.frame(cbind(x.tt,x.effect))

abundance_change_das <- x.all$diff.btw 

x.all <- cbind(abundance_change_das,x.all)


`-Log10(adj-pvalue)` <- -log10(x.all$wi.eBH)     # use wi.eBH as the pvalue
x.all$abundance_das <- abundance_das
x.all$`-Log10(adj-pvalue)` <- `-Log10(adj-pvalue)`


# draw figure
das <- x.all[(x.all$`-Log10(adj-pvalue)` >=1.301 & (x.all$abundance_change_das >=1 | x.all$abundance_change_das <=-1)),]

das$Species <- row.names(das)
das <- das[order(das$abundance_change_das),]    #13. ♦ dif.btw - median difference in clr values between S and NS groups
das$Color <- ifelse(das$abundance_change_das < 0, "Enriched in Urine", "Enriched in Vagina")  # above / below avg flag
das$Species <- factor(das$Species, levels = das$Species)  # convert to factor to retain sorted order in plot.

das$abundance_change_das[das$abundance_change_das == Inf] = 5
das$abundance_change_das[das$abundance_change_das == -Inf] = -5
das$abundance_change_das[das$abundance_change_das <= -5] = -5
das$abundance_change_das[das$abundance_change_das >= 5] = 5

theme_set(theme_bw())  # pre-set the bw theme.

ggplot(das, aes(Species, abundance_change_das)) + 
   geom_point(aes(col=Color, size=`-Log10(adj-pvalue)`)) + 
   coord_flip() +          # convert x y axis
   labs(x = 'Species', y = "Median difference in clr values")+ 
   theme(axis.title = element_text(size = 12), 
         axis.text = element_text(size = 12), 
         legend.text = element_text(size = 12), 
         legend.title = element_text(size = 12))  + theme_bw()
ggsave("das_2_2.pdf", width=5, height=1.1)
ggsave("das_2_2_2.pdf", width=6, height=5)


##### draw bar figure ########################

# find top n species
reads_table_v_new_bar <- reads_table_v_new

reads_table_v_new_bar$total_reads <- rowSums(reads_table_v_new_bar)
reads_table_v_new_bar = reads_table_v_new_bar[order(reads_table_v_new_bar$total_reads,decreasing = T),]
reads_table_v_new_bar$total_reads <- NULL
other_reads = as.data.frame(t(colSums(reads_table_v_new_bar[10:nrow(reads_table_v_new_bar),])))
reads_table_v_new_bar <- reads_table_v_new_bar[-c(10:nrow(reads_table_v_new_bar)), ]
reads_table_v_new_bar <- rbind(reads_table_v_new_bar,other_reads)
row.names(reads_table_v_new_bar)[10] = 'Others'

reads_table_v_new_bar_abundance <- matrix(data =0, ncol = ncol(reads_table_v_new_bar),nrow = nrow(reads_table_v_new_bar))

for (a in 1:ncol(reads_table_v_new)) {
   reads_table_v_new_bar_abundance[,a] <- reads_table_v_new_bar[,a] / colSums(reads_table_v_new_bar)[a]
   
}
row.names(reads_table_v_new_bar_abundance) = row.names(reads_table_v_new_bar)
colnames(reads_table_v_new_bar_abundance) = colnames(reads_table_v_new_bar) 

reads_table_v_new_bar <- as.data.frame(reads_table_v_new_bar_abundance)
reads_table_v_new_bar <- as.data.frame(t(reads_table_v_new_bar))
reads_table_v_new_bar = reads_table_v_new_bar[order(reads_table_v_new_bar$`Lactobacillus iners`,decreasing = T),]
reads_table_v_new_bar <- as.data.frame(t(reads_table_v_new_bar))
reads_table_v_new_bar <- reads_table_v_new_bar/(colSums(reads_table_v_new_bar))

plot_v <- gather(reads_table_v_new_bar)
plot_v_name <- row.names(reads_table_v_new_bar)
plot_v_name <- rep( plot_v_name , ncol(reads_table_v_new))
plot_v$Species <- plot_v_name

plot_v$Species <- as.factor(plot_v$Species)
plot_v$Species <- factor(plot_v$Species, levels = plot_v$Species[1:10])

plot_v$key <- as.factor(plot_v$key)
plot_v$key <- factor(plot_v$key, levels = colnames(reads_table_v_new_bar))

ggplot(data=plot_v, aes(x=key, y=value, fill=Species)) +
   geom_bar(stat="identity") +
   scale_fill_manual(values=c('#A8C5D3','#FC9800','#C80B0B','#FDFFBA','#55ff42','#F7A0A0',
                              '#D7A7F1','#3264B8','#bfbfbf','#3F3F3F')) +
   labs(x = 'Vaginal samples', y = "Abundance")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank())

ggsave("bar_v_2.pdf", width=9, height=3.5)


# urine
# find top n species
reads_table_u_new_bar <- reads_table_u_new

reads_table_u_new_bar$total_reads <- rowSums(reads_table_u_new_bar)
reads_table_u_new_bar = reads_table_u_new_bar[order(reads_table_u_new_bar$total_reads,decreasing = T),]
reads_table_u_new_bar$total_reads <- NULL
other_reads = as.data.frame(t(colSums(reads_table_u_new_bar[10:nrow(reads_table_u_new_bar),])))
reads_table_u_new_bar <- reads_table_u_new_bar[-c(10:nrow(reads_table_u_new_bar)), ]
reads_table_u_new_bar <- rbind(reads_table_u_new_bar,other_reads)
row.names(reads_table_u_new_bar)[10] = 'Others'

reads_table_u_new_bar_abundance <- matrix(data =0, ncol = ncol(reads_table_u_new_bar),nrow = nrow(reads_table_u_new_bar))

for (a in 1:ncol(reads_table_u_new)) {
   reads_table_u_new_bar_abundance[,a] <- reads_table_u_new_bar[,a] / colSums(reads_table_u_new_bar)[a]
   
}
row.names(reads_table_u_new_bar_abundance) = row.names(reads_table_u_new_bar)
colnames(reads_table_u_new_bar_abundance) = colnames(reads_table_u_new_bar) 

reads_table_u_new_bar <- as.data.frame(reads_table_u_new_bar_abundance)
reads_table_u_new_bar <- as.data.frame(t(reads_table_u_new_bar))
reads_table_u_new_bar = reads_table_u_new_bar[order(reads_table_u_new_bar$`Lactobacillus iners`,decreasing = T),]
reads_table_u_new_bar <- as.data.frame(t(reads_table_u_new_bar))
reads_table_u_new_bar <- reads_table_u_new_bar/(colSums(reads_table_u_new_bar))


# urine with marched order
plot_u <- gather(reads_table_u_new_bar)
plot_u_name <- row.names(reads_table_u_new_bar)
plot_u_name <- rep( plot_u_name , ncol(reads_table_u_new))
plot_u$Species <- plot_u_name

plot_u$Species <- as.factor(plot_u$Species)
plot_u$Species <- factor(plot_u$Species, levels = plot_u$Species[1:10])

plot_u$key <- as.factor(plot_u$key)
order <- colnames(reads_table_v_new_bar)
order_urine <- vector(mode = "character", length = length(order))
sample_list_2 <- sample_list[(length(order)+1):nrow(sample_list),]
for (a in 1: length(order)) {
   n = which(sample_list_2[,2] == order[a])
   order_urine[a] = sample_list_2[n,1]
}
plot_u$key <- factor(plot_u$key, levels = order_urine)

ggplot(data=plot_u, aes(x=key, y=value, fill=Species)) +
   geom_bar(stat="identity") +
   scale_fill_manual(values=c('#A8C5D3','#FC9800','#C80B0B','#D7A7F1','#55ff42','#FDFFBA',
                              '#3264B8','#F7A0A0','#bfbfbf','#3F3F3F')) +
   labs(x = 'Urine samples', y = "Abundance")+ 
   theme(axis.title = element_text(size = 16), 
         axis.text = element_text(size = 16), 
         legend.text = element_text(size = 16), 
         legend.title = element_text(size = 16),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank())
ggsave("bar_u_vagina_order_2.pdf", width=9, height=3.5)

##### linear regression  ################

line_list <- matrix(data = NA, nrow = nrow(reads_table_v_abundance), ncol =7)
line_list <- as.data.frame(line_list)
colnames(line_list) = c('Species','p-value','R-value',
                        'linear_model_function','p-value_remove_low_abundance',
                        'R-value_remove_low_abundance','linear_model_function_remove_low_abundance')
row.names(line_list)= row.names(reads_table_v_abundance)
line_list$Species = row.names(reads_table_v_abundance)

for (a in 1:nrow(line_list)) {
   # no treatment
   line <- as.data.frame(t(rbind(reads_table_u_abundance[a,],reads_table_v_abundance[a,])))
   
   line <- log10(line)
   
   colnames(line) = c('Urine','Vagina')
   
   # linear regression
   line2=line
   line2[line2 == '-Inf'] = -6
   linearMod <- lm(Urine ~ Vagina, data=line2) 
   
   b = format(round(linearMod$coefficients[2], 3), nsmall = 3)
   c = format(round(linearMod$coefficients[1], 3), nsmall = 3)
   
   if (linearMod$coefficients[1] > 0){
      line_list[a,4]= paste0('y=',b,'x+',c)
   } else if (c < 0) {
      line_list[a,4]= paste0('y=',b,'x',c)
   } else {
      line_list[a,4]= paste0('y=',b,'x')
   }
   
   r_squared <- summary(linearMod)$r.squared
   line_list[a,3]= r_squared
   pvalue <- as.numeric(summary(linearMod)$coefficients[,4][2])
   line_list[a,2]= pvalue
   
   # linear regression after removing low abundance
   keep = line2$Urine > -3 | line2$Vagina > -3 
   line3 = line2[keep,]
   
   if (sum(keep) >=5) {
      linearMod <- lm(Urine ~ Vagina, data=line3) 
      
      b = format(round(linearMod$coefficients[2], 3), nsmall = 3)
      c = format(round(linearMod$coefficients[1], 3), nsmall = 3)
      
      if (linearMod$coefficients[1] > 0){
         line_list[a,7]= paste0('y=',b,'x+',c)
      } else if (c < 0) {
         line_list[a,7]= paste0('y=',b,'x',c)
      } else {
         line_list[a,7]= paste0('y=',b,'x')
      }
      
      r_squared2 <- summary(linearMod)$r.squared
      line_list[a,6]= r_squared2
      pvalue2 <- as.numeric(summary(linearMod)$coefficients[,4][2])
      line_list[a,5]= pvalue2
   }
   
   
   if (line_list[a,1] %in% c('Sneathia amnii ','Lachnospiraceae BVAB1 ',
                             'TM7 OTU-H1 ','Prevotella cluster2 ',
                             'Gardnerella vaginalis ','Atopobium vaginae ')) {
      ggplot(line, aes(Vagina, Urine) )  +
         geom_point(size=1) +
         xlab(paste0("Urine abundance")) +
         ylab(paste0("Vagina abundance")) + ggtitle(paste0(line_list[a,1], '\n',"R = ", r_squared, '\n',"pvalue = ", pvalue))  + 
         geom_smooth(method='lm', se=T) + 
         theme(axis.title.x = element_text( size=16),
               axis.title.y = element_text( size=16)) 
      ggsave(paste0('linear_2_',line_list[a,1],'.pdf'), width=5.5, height=5.5)
   }
   
}


line_list$adj_p = p.adjust(line_list$`p-value`)  #  Benjamini & Hochberg 
line_list$adj_p_remove = p.adjust(line_list$`p-value_remove_low_abundance`)  #  Benjamini & Hochberg 


line_list = cbind(line_list,total_abundance_v, total_abundance_u)
colnames(line_list)[10] = 'Abundance in the vaginal microbiome' 
colnames(line_list)[11] = 'Abundance in the urinary microbiome' 

write.csv(line_list,'linear_regression_2.csv')


###################### metadata and Fig. 1f ###################################


##### prepare metadata  ##################################################
setwd('/Users/binzhu/secure/momspi-data/home/metadata/urine_project')
metadata_v = read.table('allVaginal.16SDB.V1-V3.uclust99.summary.txt',header = T, sep = "\t")
metadata_u = read.table('allUrine.16SDB.V1-V3.uclust99.summary.txt',header = T, sep = "\t")
metadata_r = read.delim('allRectal.16SDB.V1-V3.uclust99.summary.txt',header = T, sep = "\t")
metadata_b = read.delim('allBuccal.16SDB.V1-V3.uclust99.summary.txt',header = T, sep = "\t")

metadata = read.table('metadata.txt',header = T)
metadata2 = read.table('age.txt',header = F)
match_cv= read.table('MCKD2MV1D.txt',header = F)
match_rv= read.table('MRCD2MV1D.txt',header = F)
match_uv= read.table('MU1W2MV1D_urineChris_022520.txt',header = F)

setwd('/Users/binzhu/Desktop/Urine')

##### data arrangement ###########################################
colnames(metadata_v) = c('V1','V2','V3','V4','V5','V6')
colnames(metadata_u) = c('V1','V2','V3','V4','V5','V6')
colnames(metadata_r) = c('V1','V2','V3','V4','V5','V6')
colnames(metadata_b) = c('V1','V2','V3','V4','V5','V6')

metadata_v$V3 = str_replace_all(metadata_v$V3,'AT','')
metadata_u$V3 = str_replace_all(metadata_u$V3,'AT','')
metadata_r$V3 = str_replace_all(metadata_r$V3,'AT','')
metadata_b$V3 = str_replace_all(metadata_b$V3,'AT','')

metadata_v_2 = as.data.frame(str_c(metadata_v$V2,metadata_v$V3,sep = " "))
metadata_v = cbind(metadata_v$V1,metadata_v_2,metadata_v$V4)
rm(metadata_v_2)
#colnames(metadata_v) = c('V1','V2','V3')

metadata_u_2 = as.data.frame(str_c(metadata_u$V2,metadata_u$V3,sep = " "))
metadata_u = cbind(metadata_u$V1,metadata_u_2,metadata_u$V4)
rm(metadata_u_2)
#colnames(metadata_u) = c('V1','V2','V3')

metadata_r_2 = as.data.frame(str_c(metadata_r$V2,metadata_r$V3,sep = " "))
metadata_r = cbind(metadata_r$V1,metadata_r_2,metadata_r$V4)
rm(metadata_r_2)
#colnames(metadata_r) = c('V1','V2','V3')

metadata_b_2 = as.data.frame(str_c(metadata_b$V2,metadata_b$V3,sep = " "))
metadata_b = cbind(metadata_b$V1,metadata_b_2,metadata_b$V4)
rm(metadata_b_2)
#colnames(metadata_b) = c('V1','V2','V3')

# get reads table
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


sample_name <- unique(metadata_r$`metadata_r$V1`)
species_name <- unique(metadata_r$`str_c(metadata_r$V2, metadata_r$V3, sep = " ")`)

reads_table_r <- matrix(0, ncol = length(sample_name), nrow = length(species_name))     # create reads table   sample name as columns and species name as rows
row.names(reads_table_r) = species_name 
colnames(reads_table_r) = sample_name

for (a in 1: dim(metadata_r)[1]) {
   
   column_num <- which(sample_name == metadata_r[a,1])
   row_num <- which(species_name == metadata_r[a,2])
   
   reads_table_r[row_num,column_num] = as.numeric(as.character(metadata_r[a,3]))
}

reads_table_r = as.data.frame(reads_table_r)


sample_name <- unique(metadata_b$`metadata_b$V1`)
species_name <- unique(metadata_b$`str_c(metadata_b$V2, metadata_b$V3, sep = " ")`)

reads_table_b <- matrix(0, ncol = length(sample_name), nrow = length(species_name))     # create reads table   sample name as columns and species name as rows
row.names(reads_table_b) = species_name 
colnames(reads_table_b) = sample_name

for (a in 1: dim(metadata_b)[1]) {
   
   column_num <- which(sample_name == metadata_b[a,1])
   row_num <- which(species_name == metadata_b[a,2])
   
   reads_table_b[row_num,column_num] = as.numeric(as.character(metadata_b[a,3]))
}

reads_table_b = as.data.frame(reads_table_b)


##### get sample pair #####
sample_list = colnames(reads_table_v)
keep = (sample_list %in% match_cv$V2) & (sample_list %in% match_rv$V2) & (sample_list %in% match_uv$V2)
sum(keep)
sample_list = sample_list[keep]
match_cv = match_cv[match_cv$V2 %in% sample_list,]
match_rv = match_rv[match_rv$V2 %in% sample_list,]
match_uv = match_uv[match_uv$V2 %in% sample_list,]

sample_list = sample_list[order(sample_list)]
match_cv = match_cv[order(match_cv$V2),]
match_rv = match_rv[order(match_rv$V2),]
match_uv = match_uv[order(match_uv$V2),]

reads_table_v = reads_table_v[, colnames(reads_table_v) %in% sample_list]
reads_table_v = reads_table_v[sample_list]
reads_table_b = reads_table_b[, colnames(reads_table_b) %in% match_cv$V1]
reads_table_b = reads_table_b[match_cv$V1]
reads_table_r = reads_table_r[, colnames(reads_table_r) %in% match_rv$V1]
reads_table_r = reads_table_r[match_rv$V1]
match_uv$V1 = str_replace_all(match_uv$V1, 'MU1W','MU1D')
reads_table_u = reads_table_u[, colnames(reads_table_u) %in% match_uv$V1]
reads_table_u = reads_table_u[match_uv$V1]

##### sample total reads threshold ####################
keep = colSums(reads_table_b) >= 5000
sum(keep)
reads_table_b = reads_table_b[,keep]
reads_table_v = reads_table_v[,keep]
reads_table_r = reads_table_r[,keep]
reads_table_u = reads_table_u[,keep]

keep = colSums(reads_table_r) >= 5000
sum(keep)
reads_table_b = reads_table_b[,keep]
reads_table_v = reads_table_v[,keep]
reads_table_r = reads_table_r[,keep]
reads_table_u = reads_table_u[,keep]

keep = colSums(reads_table_u) >= 5000
sum(keep)
reads_table_b = reads_table_b[,keep]
reads_table_v = reads_table_v[,keep]
reads_table_r = reads_table_r[,keep]
reads_table_u = reads_table_u[,keep]

keep = colSums(reads_table_v) >= 5000
sum(keep)
reads_table_b = reads_table_b[,keep]
reads_table_v = reads_table_v[,keep]
reads_table_r = reads_table_r[,keep]
reads_table_u = reads_table_u[,keep]


metadata_all = c(rep('Vagina',ncol(reads_table_v)),rep('Urine',ncol(reads_table_u)),
                 rep('Rectum',ncol(reads_table_r)),rep('Cheek',ncol(reads_table_b)))


species_list = unique(c(row.names(reads_table_v),row.names(reads_table_u),
                        row.names(reads_table_r),row.names(reads_table_b)))
species_list <- species_list[order(species_list)]
reads_table_v <- reads_table_v[order(row.names(reads_table_v)), ]
reads_table_u <- reads_table_u[order(row.names(reads_table_u)), ]
reads_table_r <- reads_table_r[order(row.names(reads_table_r)), ]
reads_table_b <- reads_table_b[order(row.names(reads_table_b)), ]

reads_table_u_new <- as.data.frame(matrix(data=0,ncol=ncol(reads_table_u),nrow=length(species_list)))
reads_table_v_new <- as.data.frame(matrix(data=0,ncol=ncol(reads_table_v),nrow=length(species_list)))
reads_table_r_new <- as.data.frame(matrix(data=0,ncol=ncol(reads_table_r),nrow=length(species_list)))
reads_table_b_new <- as.data.frame(matrix(data=0,ncol=ncol(reads_table_b),nrow=length(species_list)))

row.names(reads_table_u_new) = species_list
colnames(reads_table_u_new) = colnames(reads_table_u)
row.names(reads_table_v_new) = species_list
colnames(reads_table_v_new) = colnames(reads_table_v)
row.names(reads_table_r_new) = species_list
colnames(reads_table_r_new) = colnames(reads_table_r)
row.names(reads_table_b_new) = species_list
colnames(reads_table_b_new) = colnames(reads_table_b)

keep <- species_list %in% row.names(reads_table_v)
reads_table_v_new[keep,] = reads_table_v
keep <- species_list %in% row.names(reads_table_u)
reads_table_u_new[keep,] = reads_table_u
keep <- species_list %in% row.names(reads_table_r)
reads_table_r_new[keep,] = reads_table_r
keep <- species_list %in% row.names(reads_table_b)
reads_table_b_new[keep,] = reads_table_b

reads_table <- cbind(reads_table_v_new,reads_table_u_new,reads_table_r_new,reads_table_b_new)


##### species threshold ###################################

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
   c = sum(reads_table_abundance[a,] >= 0.001) / ncol(reads_table_abundance) >= 0.05        # input
   d = sum(reads_table_abundance[a,] >= 0.0001) / ncol(reads_table_abundance) >= 0.15        # input
   keep[a] = c|d
}

reads_table <- reads_table[keep,]

##### rarefaction ##################################################

reads_table_rare <- t(reads_table)
reads_table_rare = Rarefy(reads_table_rare, depth = min(rowSums(reads_table_rare)))
reads_table_rare <- reads_table_rare$otu.tab.rff
reads_table_rare <- as.data.frame(reads_table_rare)

##### beta diversity ##################################################
# Bray_Curtis distance
Bray_Curtis <- as.matrix(vegdist(reads_table_rare, METHODS="bray", binary=FALSE))
Bray_Curtis <- as.data.frame(Bray_Curtis)

# Running Nonmetric Multidimensional Scaling (NMDS) Ordination
NMDS <-
   metaMDS(Bray_Curtis,
           distance = "bray",
           k = 2,
           maxit = 999, 
           trymax = 50,
           wascores = TRUE)

mds_data <- as.data.frame(NMDS$points)
mds_data$Site <- metadata_all

mds_data$Site = factor(mds_data$Site, levels = c('Urine','Rectum','Vagina','Cheek'))

ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Site)) +
   geom_point()+ theme_bw()
ggsave("NMDS_all.pdf", width=5.5, height=5.5)


############################### get metadata of 75 participants #########################
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

match_uv$V1 = str_replace_all(match_uv$V1,'MU1W','MU1D')
reads_table_u = reads_table_u[match_uv$V1]
reads_table_v = reads_table_v[match_uv$V2]

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

##### get metadata #####
keep = metadata$sampleID %in% colnames(reads_table_v)

metadata = metadata[keep,]

mean(metadata$ga_at_col)
sd(metadata$ga_at_col)


mean(metadata2$V2)
sd(metadata2$V2)



##### empty control #####
setwd('/Users/binzhu/secure/godel/gpfs_fs/momspi/16s/M16-001/stirrups')

file_list <- list.files(path = ".", pattern = '*Water*')
file_list <- file_list[str_detect(file_list,'hq_summary_97.txt')]


trials = c(1: length(file_list))

func_1 = function(trial) {
   if (file.info(file_list[trial])$size == 0) {
      return(NA)
   }
   
   reads_table = read.table(file_list[trial], sep = '|')
   
   reads_table$V3[reads_table$V3 == 'AT'] = ''
   reads_table$V3 = paste(reads_table$V2,reads_table$V3,sep='_')
   reads_table[,c(2,5,6)] = NULL
   
   return(reads_table)
   
}

data_3 = mclapply(trials, func_1, mc.cores = 8)

setwd('/Users/binzhu/Desktop/Urine')

data = as.data.frame(matrix(data = NA, ncol = 3, nrow = 0))
colnames(data) = c('SampleID','Taxa','No_of_Reads')

for (a in 1: length(data_3)) {
   data_2 = data_3[a]
   data_2 <- unlist(data_2)
   
   data_4 = as.data.frame(matrix(data = NA, ncol = 3, nrow = length(data_2)/3))
   colnames(data_4) = c('SampleID','Taxa','No_of_Reads')
   data_4[,1] = as.character(data_2[1:(length(data_2)/3)])
   data_4[,2] = as.character(data_2[((length(data_2)/3)+1):((length(data_2)/3)*2)])
   data_4[,3] = as.character(data_2[((length(data_2)/3)*2+1):((length(data_2)/3)*3)])
   
   data = rbind(data,data_4)
}

colnames(data) = c('V1','V2','V4')
sample_name <- unique(data$V1)
species_name <- unique(data$V2)

data$V2 <- as.character(data$V2)
data$V1 <- as.character(data$V1)
data$V4 <- as.numeric(as.character(data$V4))

reads_table <- matrix(0, ncol = length(sample_name), nrow = length(species_name))     # create reads table   sample name as columns and species name as rows
row.names(reads_table) = species_name 
colnames(reads_table) = sample_name

for (a in 1: nrow(data)) {
   
   column_num <- which(sample_name == data[a,1])
   row_num <- which(species_name == data[a,2])
   
   reads_table[row_num,column_num] =  data$V4[a]
   
}
reads_table <- as.data.frame(reads_table)
colnames(reads_table) = str_remove_all(colnames(reads_table) ,'M16-001-')
write.csv(reads_table,'Blank_control.csv')

















