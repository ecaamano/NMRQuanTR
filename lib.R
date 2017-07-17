
library(ggplot2)
library(gridExtra)

# Funtions for NMRQuantTR

parseSpinSysData <- function(data, verbose=F){
  
  X <- matrix("",nrow(data),3)
  
  Col1 <- as.character(data$Assign.F1)
  Col2 <- as.character(data$Assign.F2)
  
  for(i in 1:nrow(data)){
    temp1 <- strsplit(Col1[i], "[^[:digit:]]")[[1]]
    temp2 <- strsplit(Col2[i], "[^[:digit:]]")[[1]]
    
    temp1 <- temp1[temp1!=""]
    temp2 <- temp2[temp2!=""]
    
    if(length(temp1) == 2){
      X[i,1] <- temp1[1]
      X[i,2] <- temp1[2]
    } else if(length(temp1) == 1) {
      X[i,1] <- NA
      X[i,2] <- temp1
    } else {
      X[i,1] <- NA
      X[i,2] <- NA
    }
    
    if(length(temp2) == 2){
      X[i,3] <- temp2[2]
    } else if(length(temp2) == 1) {
      X[i,3] <- temp2
    } else {
      X[i,3] <- NA
    }
  }
  
  out_temp = matrix(as.numeric(X),nrow(data),3)
  colnames(out_temp) <- c("Sp.system","Assign.F1","Assign.F2")
  out = cbind(data[,1:4],out_temp,data[,7:8])
 
  if(verbose){
    cat("\nSpin systems in this file: ")
    cat(unique(out$Sp.system[!is.na(out$Sp.system)]))
    cat("\n\n")
  }
  
  return(out)
}

changeSpins2Names <- function(data, specNom){
  
  # A function that changes the spin system numbers to compound names in the data matrix
  # according to a number-name data frame.
  #
  # Input:  data      - the dataset after number parsing (removing brackets and splitting into 2 columns)
  #        specNom - a data frame of numbers and corresponding names. Created as shown in the example below.
  # Instead of a specNom we can use the file with all the assignments and call the second column, that has 
  # the spin systems and the first column.
  #
  # Ex.:specNom <- data.frame(nums=c(12,13,15,18,19),
  #                             names=c("qwe","asd","fgh","rty","iop"), 
  #                              stringsAsFactors=F)
  #
  # Output: dataTemp - a data frame where the numbers in the 5th column
  #                    are substituted according to the nameFrame.
  
  
  #dataTemp <- rep('', nrow(data))
  if(sum(is.na(data$Sp.system)) > 0) data <- data[-which(is.na(data$Sp.system)),]
  if(sum(is.na(data$Assign.F1)) > 0) data[is.na(data[,'Assign.F1']), 'Assign.F1'] <- 0
  if(sum(is.na(data$Assign.F2)) > 0) data[is.na(data[,'Assign.F2']), 'Assign.F2'] <- 0
  
  dataTemp <- as.character(data$Sp.system)
  
  for(i in 1:nrow(specNom)){
    dataTemp[data[,'Sp.system'] == specNom[i,'Spin_System']] <- specNom[i,'Metabolite']
  }
  data[,'Sp.system'] = dataTemp
  
  return(data)
}

normaliseByPeak <- function(data, normBy="none", exclude="none", concRef=2.5, concRefSC=2.5){
  
  # concRef - reference material concentration in the sample
  # concRefSC - reference material concentration in the samples collected for calibration curves
  
  #uniqueCarbons <- chooseUniqueSpinS(data)
  
  #exclude lines in the exclusion list
  if(length(exclude) > 1 | exclude != "none"){
    rowsToExclude <- which(data$Sp.system %in% exclude)
    if(length(rowsToExclude > 0)) data <- data[-rowsToExclude,]
    #rowsToExclude <- which(UniqueCarbons$Sp.system %in% exclude)
    #uniqueCarbons <- uniqueCarbons[-rowsToExclude,]
  }
  
  #Normalize by the line in the normBy variable
  if(normBy != "none" & concRef != "none"){
    normMat <- data[data$Sp.system == normBy, c("Height","Volume")]
    normMat$Height <- normMat$Height * (concRefSC / concRef)
    normMat$Volume <- normMat$Volume * (concRefSC / concRef)
    data <- data[-which(data$Sp.system == normBy),]
    
    data$Height <- data$Height / normMat$Height
    data$Volume <- data$Volume / normMat$Volume
    
    #rowsToExclude <- which(uniqueCarbons$Sp.system == normBy)
    #uniqueCarbons <- uniqueCarbons[-rowsToExclude,]
    
    return(data)
    
  } else {
    
    warning('normaliseByPeak() called without the normBy and/or concRef parameter. Data not normalised.')
    return(data)
  }
}

changeNomenclature <- function(data, specNom, excludeNotFound=T){
  #UniqueCarbons<-chooseUniqueSpinS(data3)
  
  data[is.na(data[,'Assign.F1']), 'Assign.F1'] = 0
  data[is.na(data[,'Assign.F2']), 'Assign.F2'] = 0
  #Create a data.frame for storage fo results
  finalmat <- data.frame(Carbon = rep('', nrow(data)), 
                         Height = data$Height, 
                         Volume = data$Volume, stringsAsFactors = F)
  
  # create a loop that for each row  it first checks if the numbers are
  # a given file and substitute their names for the value
  # and if not found, it takes the numbers and it says which numbers we have, 
  # the user gives an input with the carbon name and the latter with the sp system,
  # heigh and volume are stored in the matrix

  #databasepeaks <- read.csv(path, header=T, stringsAsFactors=F)
  
  for(i in 1:nrow(data)){
    temp <- specNom[(specNom[,3] == data$Assign.F1[i]) | (specNom[,4] == data$Assign.F2[i]),c(1,5)]

    if(nrow(temp) == 1){
      carbonName <- paste(temp[1], temp[2], sep="_")
      finalmat[i,1] <- carbonName
    } else if(data$Sp.system[i] == "0"){
      finalmat[i,1] <- "exclude"
    } else if(!excludeNotFound){
        carbon <- paste(data$Sp.system[i], data$Assign.F1[i], data$Assign.F2[i])
        cat(sprintf("\n Line %s. Enter name for the carbon and press <Enter>.\n Enter 'exclude' for cabon to be excluded. \n", carbon))
        carbonName <- readline(prompt="Name:")
        if(carbonName != 'exclude') carbonName <- paste(data$Sp.system[i],"_",carbonName,sep="")
        finalmat[i,1]<-carbonName
    } else{
      finalmat[i,1] <- 'exclude'
    }
  }
  
  linestoexclude <- which(finalmat[,1] == 'exclude')
  
  if(length(linestoexclude) != 0) finalmat <- finalmat[-linestoexclude,]
  
  return(finalmat)
}

removePeaksByName <- function(data, name){

    remRows <- grep(name, data[,1])
    if(length(remRows) > 0) return(data[-remRows,])
    else return(data) 
}

setNamesAndOrdering <- function(data){
  colnames(data) <- c("Carbon", "Height", "Volume")
  data <- data[order(data[,1]),]
  return(data)
}

calcConc <- function(data, coefficients){
  
  heightConc <- data.frame(data[,1], rep(0, nrow(data)), stringsAsFactors = F)
  volumeConc <- data.frame(data[,1], rep(0, nrow(data)), stringsAsFactors = F)
  
  #if(!is.null(path)) coefficients <- read.csv(path, header=T, stringsAsFactors = F)
  #else stop('Path not given to calcConc()')
  
  if(sum(duplicated(coefficients[,1])) > 0) stop('Duplicated entries in the coefficients file.')
  
  for(i in 1:nrow(data)){
    
    hCoeff <- coefficients[coefficients[,1] == data[i,1], 2]
    if(length(hCoeff != 0)) heightConc[i,2] <- data[i,2] / hCoeff
    else message(sprintf("Equation for height - %s not found..\n", data[i, 1]))
    
    vCoeff <- coefficients[coefficients[,1] == data[i,1], 3]
    if(length(vCoeff != 0)) volumeConc[i,2] <- data[i,3] / vCoeff
    else message(sprintf("Equation for volume - %s not found..\n", data[i,1]))
  }
  
  output <- heightConc
  output$Volume <- volumeConc[,2]
  colnames(output)<-c("Carbon","Conc_Height","Conc_Volume")
  return(output) 
}

adjustAliquot <- function(data, vNMRtube, vAliquot){
  # adjust the concentration of sample, (e.g. we take 800µL sample but freeze dry and resuspend in 300µL)
  data[,2:3] <- data[,2:3] * vNMRtube / vAliquot
  return(data)
}

pickOnePeak <- function(data, path=NULL){
  #In path, a file with two columns, one with the carbon we want to pick and other with how we want to rename it
  #In this case the file is in dropbox and called carbonToMetaboliteAllMetabs.csv
  if(!is.null(path)) carbons<-read.csv(path, header=T, stringsAsFactors=F)
  else stop('pickOnePeak() called without path parameter')
  
  for(i in 1:nrow(data)){
    for(j in 1:nrow(carbons)){
      if (data[i,1] == carbons[j,1]){
        data[i,1] <- carbons[j,2]
        break
      }
    }
  }
  remRows <- grep("_", data[,1])
  data <- data[-remRows,]
  return(data) 
} 

readFolder <- function(path){
   files <- list.files(path)
   X <- list()
   out <- data.frame()
   
   for(i in 1:length(files)){
      filePath <- paste(path,"/",files[i],sep="")
      dataTemp <- read.csv(filePath,header=T)
      dataTemp <- conv4curves(dataTemp)
      
      conc <- strsplit(files[i], "[^[:digit:] & ^.]")
      conc <- unlist(conc)
      conc <- conc[conc!=""]
      conc <- rep(conc,nrow(dataTemp))
      dataTemp$Conc <- as.factor(conc)
            
      X[[i]]<-dataTemp
   }
   
   for(i in 1:length(files)){
      out <- rbind(out,X[[i]])
   }
   cat("\nSpin systems in this folder: ")
   cat(unique(out$Sp.system[!is.na(out$Sp.system)]))
   cat("\n")
   return(out)
}

chooseUniqueSpins <- function(data){
   dat<-data[c("Sp.system","Assign.F1","Assign.F2")]
   unis<-unique(dat[,c("Sp.system","Assign.F1","Assign.F2")])
   unis<-unis[order(unis$Sp.system),]
   return(unis)
}

parseConfig <- function(configFile){
  conf = list(pathToNomenclFile = '',
              pathToCalCurveFile = '',
              normBy = '',
              remLst = '',
              vAliquot = 1,
              vNMRsample = 1,
              concRef = 1,
              concRefSC = 1,
              pathOut = getwd())
  
  conff <- file(configFile, 'r')
  while(TRUE){
    line = readLines(conff, n=1)
    if(length(line) != 0){
      line = strsplit(line, '=')[[1]]
      line = gsub(' ', '', line)
      conf[line[1]] <- line[2]
    } else {
      break
    }
  }
  close(conff)
  conf['vAliquot'] <- as.numeric(conf['vAliquot'])
  conf['vNMRsample'] <- as.numeric(conf['vNMRsample'])
  conf['concRef'] <- as.numeric(conf['concRef'])
  conf['concRefSC'] <- as.numeric(conf['concRefSC'])
  conf
}

# Plots

plotBarsWithErrors <- function(data, title, yheight, savePlots=T){
  
  names(data) <- c('Metabolite', 'Height', 'Volume')
  meansHeight <- tapply(data$Height,data$Metabolite, mean)
  Metabolite <- names(meansHeight)
  meansHeight <- as.numeric(as.character(meansHeight))
  meansVolume <- as.numeric(as.character(tapply(data$Volume,data$Metabolite, mean)))
  seHeight <- as.numeric(as.character(tapply(data$Height, data$Metabolite, function(x){sd(x)/sqrt(length(x))})))
  seVolume <- as.numeric(as.character(tapply(data$Volume, data$Metabolite, function(x){sd(x)/sqrt(length(x))})))
  
  #seHeight <- as.numeric(as.character(tapply(data$Height, data$Metabolite, function(x){sd(x)})))
  #seVolume <- as.numeric(as.character(tapply(data$Volume, data$Metabolite, function(x){sd(x)})))
  
  
  dataAll <- data.frame(Metabolite = Metabolite,
                        meansHeight = meansHeight,
                        meansVolume = meansVolume,
                        seHeight = seHeight,
                        seVolume = seVolume)
  
  errorsHeight <- ggplot(data = dataAll, aes(x=Metabolite, y=meansHeight, fill=Metabolite))+ 
    geom_bar(position=position_dodge(), stat="identity")+
    theme_bw()+
    geom_errorbar(aes(ymin=meansHeight-seHeight, ymax=meansHeight+seHeight),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))+
    ylim(0,yheight)+
    ylab("Concentration (mM)")+
    xlab("Metabolite")+
    theme(axis.title.y=element_text(face="bold",size=14),
          axis.title.x=element_text(face="bold",size=14),
          axis.text.x=element_text(angle=90,vjust=0,hjust=1,size=13,colour="black"),
          axis.text.y=element_text(size=14,colour="black"),legend.position="none")+
    ggtitle(title)+
    theme(plot.title=element_text(size=16,face="bold"))+
    geom_text(aes(y=meansHeight+3,label=round(meansHeight,1.5),angle=45,vjust=0.8,hjust=0.5),size=5)
  
  if(savePlots) ggsave(plot=errorsHeight, file=paste(title,"Height.pdf",sep=""))
  
  errorsVolume <- ggplot(data = dataAll, aes(x=Metabolite, y=meansVolume, fill=Metabolite))+ 
    geom_bar(position=position_dodge(), stat="identity")+
    theme_bw()+
    geom_errorbar(aes(ymin=meansVolume-seVolume, ymax=meansVolume+seVolume),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))+
    ylim(0,yheight)+
    ylab("Concentration (mM)")+ 
    xlab("Metabolite")+
    theme(axis.title.y=element_text(face="bold",size=14),axis.title.x=element_text(face="bold",size=14),
          axis.text.x=element_text(angle=90,vjust=0,hjust=1,size=13,colour="black"),axis.text.y=element_text(size=14,colour="black"),legend.position="none")+
    ggtitle(title)+
    theme(plot.title=element_text(size=16,face="bold"))+
    geom_text(aes(y=meansVolume+3,label=round(meansVolume,1.5),angle=45,vjust=0.8,hjust=0.5),size=5)
  
  if(savePlots) ggsave(plot=errorsVolume, file=paste(title,"Vol.pdf",sep=""))
}

plotNMRCalCurves <- function(data, normBy="none", exclude="none", saveAsPdf=F){
   
   #============================================================================================
   #The main plotting function. Takes a data frame outputted from readFolder and plots calibration curves
   #for each carbon peak height and volume with line parameters and 95% confidense interval.
   # Input:    data1       - (data.frame) Data extracted with readFolder()
   #           normBy      - (string) Compound name to normalize by
   #           exclude     - (string) a vector of compound names to be excluded from plotting
   #           saveAsPdf   - (boolean) spec. wether to save plots as PDFs(True) or show on the screen(False)
   #
   # Output:   Plots, either in the screen or files in the current directory.
   #============================================================================================
  
   UniqueCarbons <- chooseUniqueSpinS(data)
      
   #exclude row in the exclusion list
   if( length(exclude) > 1 | exclude!="none"){
      linesToExclude <- which(data$Sp.system %in% exclude)
      data <- data[-linesToExclude,]
      linesToExclude <- which(UniqueCarbons$Sp.system %in% exclude)
      UniqueCarbons <- UniqueCarbons[-linesToExclude,]
   }
   
   #Normalize by the row in the normBy variable
   if(normBy != "none"){
      normMat <- data[data$Sp.system == normBy, c("Height","Volume")]
      data <- data[-which(data$Sp.system == normBy),]
      
      uniqueConc <- unique(data$Conc)
      
      data[,"Height"]<-data[,"Height"]/mean(normMat[,"Height"])
      data[,"Volume"]<-data[,"Volume"]/mean(normMat[,"Volume"])
      
      linesToExclude<-which(UniqueCarbons$Sp.system == normBy)
      UniqueCarbons <- UniqueCarbons[-linesToExclude,]
   }
   
   outputCoeffs <- matrix(0, nrow(UniqueCarbons),2)
   outputNames <- rep("", nrow(UniqueCarbons)) 
   
   for(i in 1:nrow(UniqueCarbons)){
      dataTemp <- data[data$Sp.system==UniqueCarbons[i,1] & data$Assign.F1==UniqueCarbons[i,2] & data$Assign.F2==UniqueCarbons[i,3],c("Height","Volume","Conc")]
      dataTemp$Conc <- as.numeric(as.character(dataTemp$Conc))
      dataTemp <- dataTemp[order(dataTemp$Conc),]
      
      #Making the variables to be plotted:
      carbon <- paste(UniqueCarbons[i,1], UniqueCarbons[i,2], UniqueCarbons[i,3])
      
      print(dataTemp)
      mod1 <- lm(Height~Conc-1, data=dataTemp)
      coeffsHeight <- mod1$coefficients
      coeffsHeight <- paste(as.character(round(coeffsHeight, 4)))
      
      mod2 <- lm(Volume~Conc-1, data=dataTemp)
      coeffsVol <- mod2$coefficients
      coeffsVol <- paste(as.character(round(coeffsVol, 4)))   
      
      outputCoeffs[i, 1:2] <- c(mod1$coefficients, mod2$coefficients)
      
      cat(sprintf("\n Now plotting %s. Enter name for the carbon.\n",carbon))
      carbonName <- readline(prompt="Name: ")
      outputNames[i] <- paste(UniqueCarbons[i,1],carbonName, sep="_") # this is for the output
      carbonName <- paste(UniqueCarbons[i,1],carbonName, sep=" ")
      
      #Building the plot
   
      Eq1 <- paste(" y = ", coeffsHeight,"*x ", sep="")
      Rsq1 <- paste("R^2: ", round(summary(mod1)$adj.r.squared,4), sep="")
   
      Eq2 <- paste("y = ", coeffsVol,"*x ", sep="")
      Rsq2 <- paste("R^2: ", round(summary(mod2)$adj.r.squared,4), sep="")
   
      GG1<-ggplot(data=dataTemp, aes(x=Conc,y=Height))+
      geom_point(size=3)+  
      stat_smooth(method="lm")+
      annotate("rect", xmin = 0, xmax = 6, ymin = max(dataTemp$Height)-0.095*max(dataTemp$Height), ymax = max(dataTemp$Height)+0.095*max(dataTemp$Height), alpha = .2)+
      annotate("text", x = 3, y = max(dataTemp$Height)+0.04*max(dataTemp$Height), label=Eq1)+
      annotate("text", x = 3, y = max(dataTemp$Height)-0.03*max(dataTemp$Height), label=Rsq1, parse=T)+
      ggtitle(carbonName)+
      xlab("Concentration (mM)")
   
      GG2<-ggplot(data=dataTemp, aes(x=Conc,y=Volume))+
      geom_point(size=3)+
      stat_smooth(method="lm")+
      annotate("rect", xmin = 0, xmax = 6, ymin = max(dataTemp$Volume)-0.095*max(dataTemp$Volume), ymax = max(dataTemp$Volume)+0.095*max(dataTemp$Volume), alpha = .2)+
      annotate("text", x = 3, y = max(dataTemp$Volume)+0.04*max(dataTemp$Volume), label=Eq2)+
      annotate("text", x = 3, y = max(dataTemp$Volume)-0.03*max(dataTemp$Volume), label=Rsq2, parse=T)+
      ggtitle(carbonName)+
      xlab("Concentration (mM)")
      
   
      if(saveAsPdf){
         
         fileName1=paste(carbon[1],"_Height_normBy_",normBy,".pdf",sep="")
         ggsave(plot=GG1, file=fileName1)
                  
         fileName2=paste(carbon[1],"_Vol_normBy_",normBy,".pdf",sep="")
         ggsave(plot=GG2, file=fileName2)
         
      } else {
         grid.arrange(GG1,GG2, ncol=2)
      }
      
   }
   out <- data.frame(names = outputNames,
                   coeffHeight = outputCoeffs[,1],
                   coeffVolume = outputCoeffs[,2],
                   stringsAsFactors = FALSE)
   #out<-as.data.frame(UniqueCarbons)
   #out$coeffHeight <- outputCoeffs[,1]
   #out$coeffVolume <- outputCoeffs[,2]
   return(out)
}

