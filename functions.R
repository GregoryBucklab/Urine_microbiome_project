urine_vagna_analyses = function(sample_list, reads_table, reads_table_u_new , reads_table_v_new) {
#  reads_table = reads_table_all

  ##### rarefaction # that is to normalize total reads of every sample to the same level for diversity analysis ##################################################
  
  reads_table_rare <- t(reads_table)
  reads_table_rare = Rarefy(reads_table_rare, depth = min(rowSums(reads_table_rare)))
  reads_table_rare <- reads_table_rare$otu.tab.rff
  reads_table_rare <- as.data.frame(reads_table_rare)
  
  
  
  ##### alpha diversity ##################################################
  pvalue = vector()
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
    geom_point()+
    stat_ellipse(type = "t")
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
  
  ggplot(mds_data, aes(MDS1, MDS2))  +
    geom_line(aes(group = Pair), size = 0.3, alpha = 0.7) +
    geom_point(size=3,aes(shape = Site, color = Vagitype)) +
    scale_color_manual(values=c('lightblue','orange','red','brown','yellow','pink','grey')) +
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
  x.all$abundance_change_das <- abundance_change_das
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
  ggsave("das.pdf", width=8, height=1.2)
  ggsave("das2.pdf", width=8, height=4)
  
  
  
  
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
  
  
  ##### create abundance table ##################
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
  
  
  
}
