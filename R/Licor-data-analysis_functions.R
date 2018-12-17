
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
  the.data1 <- the.data[which(the.data$A > 5), ]
  #make clean and filtered csv 
  write.csv(the.data1, file = csvname, row.names = FALSE)
}


graphPack <- function(cleandata){
  #graph stomatal conductance for data put into function 
  gsw <- ggplot(cleandata, aes_string(x = "elapsed", y = "gsw")) +
    geom_point() + scale_shape_manual(values = c(4)) +
    labs(x = expression(paste("Time (seconds)")),
         y = expression(paste("gsw"))) +
    theme(axis.line.x = element_line(), axis.line.y = element_line(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)) +
    theme(legend.position = ('none')) + geom_point() +
    theme(legend.background = element_blank(),
          legend.key = element_blank(),
          legend.title = element_blank())
  
  #graph net photosynthesis for data put into function
  anet <- ggplot(cleandata, aes_string(x = "elapsed", y = "A")) +
    geom_point() + scale_shape_manual(values = c(4)) +
    labs(x = expression(paste("Time (seconds)")),
         y = expression(paste("A"))) +
    theme(axis.line.x = element_line(), axis.line.y = element_line(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)) +
    theme(legend.position = ('none')) + geom_point() +
    theme( legend.background = element_blank(),
           legend.key = element_blank(),
           legend.title = element_blank())
  #return gsw and anet graphs
  return(list(gswgraph_gsw = gsw, agraph_anet = anet))
}


# function #3 comparison analysis 
#ttest and graph average steadystates

analysis <- function(genotype1, genotype2){
  #function for ten data points for each genotype
  points <- function(x){
  # This code takes the last 10 data points of each CO2 level for supplied data
  #points at 400ppm
  fourppm <- x[which(x$elapsed > 1140 &
                         x$elapsed < 1681), ]
  # seperates A and gsw data
  fourppma <- fourppm[,4]
  fourppmgsw <- fourppm[,9]
  #points at 800ppm
  eightppm <- x[which(x$elapsed > 2940 &
                          x$elapsed < 3481), ]
  #seperates A and gsw data
  eightppma <- eightppm[,4]
  eightppmgsw <- eightppm[,9]
    
  #points at 100ppm 
  oneppm <- x[which(x$elapsed > 6543 &
                        x$elapsed < 7142), ]
  #seperates A and gsw data
  oneppma <- oneppm[,4]
  oneppmgsw <- oneppm[,9]
    
  return(list(fourppma = fourppma, fourppmgsw = fourppmgsw,
              eightppma = eightppma, eightppmgsw = eightppmgsw,
              oneppma = oneppma, oneppmgsw = oneppmgsw))
  }
  
data <- points(genotype1)
data2 <- points(genotype2)
  
#name each value in list produced from function above
names(data) <- paste(c("fourppma_1", "fourppmgsw_1",
                       "eightppma_1","eightppmgsw_1",
                       "oneppma_1", "oneppmgsw_1"))
names(data2) <- paste(c("fourppma_2", "fourppmgsw_2",
                        "eightppma_2","eightppmgsw_2",
                        "oneppma_2", "oneppmgsw_2"))
  
#function to calculate mean and standard error for each CO2 level
mean_se <- function(y){
  # calculate mean
  mean <- lapply(y, mean)
  #calculate st deviation
  sd <- lapply(y, sd)
  # unlist data
  sd1_seperated <- unlist(sd, use.names = TRUE)
  mean_seperated <- unlist(mean, use.names = TRUE)
  #calculate st error
  se <- sd1_seperated/sqrt(10)
    
  return(list(mean_seperated, se))
  }
  
data_mean <- mean_se(data)
data2_mean <- mean_se(data2)
  
# function to test for significance at each CO2 level between the two genotypes 
ttests <- function(z,n){
  dataa_z <- z[c(1, 3, 5)]
  dataa_n <- n[c(1, 3, 5)]
    
  datagsw_z <- z[c(2, 4, 6)]
  datagsw_n <- n[c(2, 4, 6)]
    
  # ttest to compare gsw and Anet of genotype1 to genotype 2 at each CO2 steady state (400, 800, 100)
  ttest1 <- t.test(datagsw_z$fourppmgsw_1, datagsw_n$fourppmgsw_2)
  ttest2 <- t.test(datagsw_z$eightppmgsw_1, datagsw_n$eightppmgsw_2)
  ttest3 <- t.test(datagsw_z$oneppmgsw_1, datagsw_n$oneppmgsw_2)
    
  ttest4 <- t.test(dataa_z$fourppma_1, dataa_n$fourppma_2)
  ttest5 <- t.test(dataa_z$eightppma_1, dataa_n$eightppma_2)
  ttest6 <- t.test(dataa_z$oneppma_1, dataa_n$oneppma_2)
    
  return(list(ttest1 = ttest1, ttest2 = ttest2,
              ttest3 = ttest3, ttest4 = ttest4,
              ttest5 = ttest5, ttest6 = ttest6))
  
  }
  
ttestreturn <- ttests(data, data2)
  
# function to graph means and st errors calculated above 
graphaverages <- function(meangen1, meangen2){
  # unlist data
  mean_seperate <- unlist(meangen1, use.names = TRUE)
  mean_seperate_2 <- unlist(meangen2, use.names = TRUE)
    
  #pull only the A data from the vector
  data_a <- mean_seperate[c(1, 3, 5)]
  data2_a <- mean_seperate_2[c(1, 3, 5)]
  # pull only the st error from the vector
  data_ase <- mean_seperate[c(7, 9, 11)]
  data2_ase <- mean_seperate_2[c(7, 9, 11)]
  # combine data from each genotype
  A <- append(data_a, data2_a)
  allse_data <- append(data_ase, data2_ase)
  #pull only the gsw data from the vector
  data_gsw <- mean_seperate[c(2, 4, 6)]
  data2_gsw <- mean_seperate_2[c(2, 4, 6)]
  # pull only the st error from the vector
  data_gswse <- mean_seperate[c(8, 10, 12)]
  data2_gswse <- mean_seperate_2[c(8, 10, 12)]
  #combine gsw data from each genotype 
  gsw <- append(data_gsw, data2_gsw)
  allse_datagsw <- append(data_gswse, data2_gswse)
    
  #make a data frame to graph Anet for genotype 1 and 2
  ppm <- c("400", "800", "100", "400", "800", "100")
  id <- c("gen1_400ppm", "gen1_800ppm", "gen1_100ppm",
          "gen2_400ppm", "gen2_800ppm", "gen2_100ppm")
  #A mean values pulled from code above
  #A st error values pulled from code above
  #make data.frame to graph
  a_graph.data <- data.frame(ppm, id, A, allse_data)
    
  # produce graph of Anet values for genotype 1 and 2
  a_graph <- ggplot(data = a_graph.data, aes(x = id, y = A, fill = ppm)) +
    geom_bar(stat = "identity") +
    scale_x_discrete(limits = c("gen1_400ppm", "gen1_800ppm", "gen1_100ppm",
                                "gen2_400ppm", "gen2_800ppm", "gen2_100ppm")) +
    geom_errorbar(aes(ymin = A - allse_data, ymax = A + allse_data), width = .2,
                  position = position_dodge(.9)) +
    theme(axis.line.x = element_line(), axis.line.y = element_line(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.text = element_text(size = 12),axis.title = element_text(size = 12)) +
    theme(legend.position = ('none')) +
    scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"))
    
    
  #make a data frame to graph gsw for genotype 1 and 2
  ppm <- c("400", "800", "100", "400", "800", "100")
  id <- c("gen1_400ppm", "gen1_800ppm", "gen1_100ppm",
          "gen2_400ppm", "gen2_800ppm", "gen2_100ppm")
  #gsw mean values pulled from code above
  #gsw st error values pulled from code above
  #make data.frame to graph
  gsw_graph.data <- data.frame(ppm, id, A, allse_datagsw)
    
  # produce graph of gsw values for genotype 1 and 2
  gsw_graph <- ggplot(data = gsw_graph.data, aes(x = id, y = gsw, fill = ppm)) +
    geom_bar(stat = "identity") +
    scale_x_discrete(limits = c("gen1_400ppm", "gen1_800ppm", "gen1_100ppm",
                                "gen2_400ppm", "gen2_800ppm", "gen2_100ppm")) +
    geom_errorbar(aes(ymin = gsw - allse_datagsw, ymax = gsw + allse_datagsw),
                  width = .2, position = position_dodge(.9)) +
    theme(axis.line.x = element_line(), axis.line.y = element_line(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),panel.background = element_blank(),
          axis.line = element_blank(), axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)) +
    theme(legend.position = ('none')) +
    scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"))
    
  return(list(a_graph = a_graph, gsw_graph = gsw_graph))
  }
  
 graphs <- graphaverages(data_mean, data2_mean)
  
 return(list(ttestreturn, graphs ))
}


