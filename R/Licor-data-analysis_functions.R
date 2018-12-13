
cleandata <- function(data, csvname){

  #remove unwanted rows and columns
  the.data <- data[-c(1:13, 15), -c(2, 4:6, 22, 23, 25:36, 39:44, 50, 51,
                                    54:168, 170:172)]
  #set column names and remove from row 1
  names(the.data) <- as.character(unlist(the.data[1,]))
  the.data <- the.data[-1,]
  #set data as numeric 
  the.data[] <- lapply(the.data, function(x) as.numeric(as.character(x)))
  #remove bad data points from match errors
  the.data1 <- the.data[the.data$A > 5, ]
  #make clean and filtered csv 
  write.csv(the.data1, file = csvname, row.names = FALSE)
  return(the.data1)
}


graphPack <- function(cleandata){

  #graph stomatal conductance for data put into function 
  gsw <- ggplot(cleandata, aes_string(x="elapsed", y="gsw"))+geom_point()+
    scale_shape_manual(values=c(4))+ labs(x=expression(paste("Time (seconds)")),
                                          y=expression(paste("gsw")))+
    theme(axis.line.x=element_line(), axis.line.y=element_line(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.text=element_text(size=12),axis.title=element_text(size=12))+
    theme(legend.position=('none'))+geom_point()+ theme(legend.background = element_blank(),
                                                        legend.key = element_blank(),
                                                        legend.title=element_blank())
  
  #graph net photosynthesis for data put into function
  anet <- ggplot(cleandata, aes_string(x="elapsed", y="A"))+geom_point()+
    scale_shape_manual(values=c(4))+ labs(x=expression(paste("Time (seconds)")),
                                          y=expression(paste("A")))+
    theme(axis.line.x=element_line(), axis.line.y=element_line(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.text=element_text(size=12),axis.title=element_text(size=12))+
    theme(legend.position=('none'))+geom_point()+ theme( legend.background = element_blank(),
                                                         legend.key = element_blank(),
                                                         legend.title=element_blank())
  #return gsw and anet graphs
  return(list(gsw, anet))
}


analysis <- function(genotype1, genotype2){

  # This code takes the last 10 data points of each CO2 level for genotype 1
  fourppm_1 <- genotype1[which(genotype1$elapsed > 1140 &
                                 genotype1$elapsed < 1681), ]
  eightppm_1 <- genotype1[which(genotype1$elapsed > 2940 &
                                  genotype1$elapsed < 3481), ]
  oneppm_1 <- genotype1[which(genotype1$elapsed > 1140 &
                                genotype1$elapsed < 1681), ]
  # Average gsw and Anet for each set of 10 data points 
  meanfourgsw_1 <- mean(fourppm_1[,"gsw"])
  meaneightgsw_1 <- mean(eightppm_1[,"gsw"])
  meanonegsw_1 <- mean(oneppm_1[,"gsw"])
  
  meanfoura_1 <- mean(fourppm_1[,"A"])
  meaneighta_1 <- mean(eightppm_1[,"A"])
  meanonea_1 <- mean(oneppm_1[,"A"])
  
  
  # This code takes the last 10 data points of each CO2 level for genotype 2
  fourppm_2 <- genotype2[which(genotype2$elapsed > 1140 &
                                 genotype2$elapsed < 1681), ]
  eightppm_2 <- genotype2[which(genotype2$elapsed > 2940 &
                                  genotype2$elapsed < 3481), ]
  oneppm_2 <- genotype2[which(genotype2$elapsed > 1140 &
                                genotype2$elapsed < 1681), ]
  # Average gsw and Anet for each set of 10 data points 
  meanfourgsw_2 <- mean(fourppm_2[,"gsw"])
  meaneightgsw_2 <- mean(eightppm_2[,"gsw"])
  meanonegsw_2 <- mean(oneppm_2[,"gsw"])
  
  meanfoura_2 <- mean(fourppm_2[,"A"])
  meaneighta_2 <- mean(eightppm_2[,"A"])
  meanonea_2 <- mean(oneppm_2[,"A"])
  
  #Standard errors for gsw and Anet of genotype 1
  
  se_fourgsw_1 <- (sd(fourppm_1[,"gsw"]))/
    sqrt(length(fourppm_1$gsw))
  se_eightgsw_1 <- (sd(eightppm_1[,"gsw"]))/
    sqrt(length(eightppm_1$gsw))
  se_onegsw_1 <- (sd(oneppm_1[,"gsw"]))/
    sqrt(length(oneppm_1$gsw))
  
  se_foura_1 <- (sd(fourppm_1[,"A"]))/
    sqrt(length(fourppm_1$A))
  se_eighta_1 <- (sd(eightppm_1[,"A"]))/
    sqrt(length(eightppm_1$A))
  se_onea_1 <- (sd(oneppm_1[,"A"]))/
    sqrt(length(oneppm_1$A))
  
  #Standard errors for gsw and Anet of genotype 2
  se_fourgsw_2 <- (sd(fourppm_2[,"gsw"]))/
    sqrt(length(fourppm_2$gsw))
  se_eightgsw_2 <- (sd(eightppm_2[,"gsw"]))/
    sqrt(length(eightppm_2$gsw))
  se_onegsw_2 <- (sd(oneppm_2[,"gsw"]))/
    sqrt(length(oneppm_2$gsw))
  
  se_foura_2 <- (sd(fourppm_2[,"A"]))/
    sqrt(length(fourppm_2$A))
  se_eighta_2 <- (sd(eightppm_2[,"A"]))/
    sqrt(length(eightppm_2$A))
  se_onea_2 <- (sd(oneppm_2[,"A"]))/
    sqrt(length(oneppm_2$A))
  
 # ttest to compare gsw and Anet of genotype1 to genotype 2 at each CO2 steady state (400, 800, 100)
  ttest1 <- t.test(fourppm_1$gsw, fourppm_2$gsw)
  ttest2 <- t.test(eightppm_1$gsw, eightppm_2$gsw)
  ttest3 <- t.test(oneppm_1$gsw, oneppm_2$gsw)
  
  ttest4 <- t.test(fourppm_1$A, fourppm_2$A)
  ttest5 <- t.test(eightppm_1$A, eightppm_2$A)
  ttest6 <- t.test(oneppm_1$A, oneppm_2$A)
  
  #bar graphs of gsw and Anet for each genotype showing averages with error bars  
  #make a data frame to graph gsw for genotype 1 and 2
  gsw_ppm <- c("400", "800", "100", "400", "800", "100")
  gsw_id <- c("gen1_400ppm", "gen1_800ppm", "gen1_100ppm",
              "gen2_400ppm", "gen2_800ppm", "gen2_100ppm")
  gsw <- c(meanfourgsw_1, meaneightgsw_1, meanonegsw_1,
           meanfourgsw_2, meaneightgsw_2, meanonegsw_2)
  gsw_se <- c(se_fourgsw_1, se_eightgsw_1, se_onegsw_1,
              se_fourgsw_2, se_eightgsw_2, se_onegsw_2)
  gsw_graph.data <- data.frame(gsw_ppm, gsw_id, gsw, gsw_se)
  
  # produce graph of gsw values for genotype 1 and 2
  gsw_graph <- ggplot(data=gsw_graph.data, aes(x=gsw_id, y=gsw, fill= gsw_ppm)) +
    geom_bar(stat="identity")+
    scale_x_discrete(limits = c("gen1_400ppm", "gen1_800ppm", "gen1_100ppm",
                                "gen2_400ppm", "gen2_800ppm", "gen2_100ppm"))+
    geom_errorbar(aes(ymin=gsw-gsw_se, ymax=gsw+gsw_se), width=.2,
                  position=position_dodge(.9))+
    theme(axis.line.x=element_line(),axis.line.y=element_line(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.text=element_text(size=12),axis.title=element_text(size=12))+
    theme(legend.position=('none'))+
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  
  #make a data frame to graph Anet for genotype 1 and 2
  A_ppm <- c("400", "800", "100", "400", "800", "100")
  A_id <- c("gen1_400ppm", "gen1_800ppm", "gen1_100ppm", "gen2_400ppm",
            "gen2_800ppm", "gen2_100ppm")
  A <- c(meanfoura_1, meaneighta_1, meanonea_1, meanfoura_2, meaneighta_2,
         meanonea_2)
  A_se <- c(se_foura_1, se_eighta_1, se_onea_1, se_foura_2, se_eighta_2,
            se_onea_2)
  A_graph.data <- data.frame(A_ppm, A_id, A, A_se)
  
  # produce graph of Anet values for genotype 1 and 2
  a_graph <- ggplot(data=A_graph.data, aes(x=A_id, y=A, fill= A_ppm)) +
    geom_bar(stat="identity")+
    scale_x_discrete(limits = c("gen1_400ppm", "gen1_800ppm", "gen1_100ppm",
                                "gen2_400ppm", "gen2_800ppm", "gen2_100ppm"))+
    geom_errorbar(aes(ymin=A-A_se, ymax=A+A_se), width=.2,
                  position=position_dodge(.9))+ theme(axis.line.x=element_line(),
                                                      axis.line.y=element_line(),
                                                      panel.grid.major = element_blank(),
                                                      panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_blank(), axis.text=element_text(size=12),axis.title=element_text(size=12))+ theme(legend.position=('none'))+
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  
  return(list(ttest1, ttest2, ttest3, ttest4, ttest5, ttest6, gsw_graph, a_graph))
}
