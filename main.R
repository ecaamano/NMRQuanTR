
# Main pipeline file
library(ggplot2)
library(gridExtra)

MainPipeline <- function(pathToData, pathToConfig, pipeline='quant', verbose=F){
  
  config = parseConfig(pathToConfig)
  specNom = read.csv(config[['pathToNomencFile']], header=T, stringsAsFactors = F)
  calCurveFile = read.csv(config[['pathToCalCurveFile']], header=T, stringsAsFactors = F)
  if(config[['remLst']] != '') {
    remLst <- strsplit(config[['remLst']], split = ',')[[1]]
    remLst <- gsub(' ','', remLst)
  }
  
  if(conf[['pathOut']] != ''){
    if(!dir.exists(conf[['pathOut']])) dir.create(conf[['pathOut']], showWarnings = T)
    outDir <- conf[['pathOut']]
  } else {
    outDir = getwd()
  }
  
  if(pipeline =='quant'){
    out <- runQuantPipeline(path=pathToData,
                     nomenclFile = specNom,
                     calCurveFile = calCurveFile,
                     normBy = conf[['normBy']],
                     remLst = remLst,
                     vAliquot = conf[['vAliquot']],
                     vNMRsample = conf[['vNMRsample']],
                     concRef = conf[['concRef']],
                     concRefSC = conf[['concRefSC']],
                     verbose=verbose)
    
  } else if(pipeline == 'makecc'){
    out <- runCalCurvePipeline(pathToData = pathToData,
                        nomenclFile = specNom,
                        normBy = conf[['normBy']],
                        saveAsPdf = T,
                        saveTo = conf[['pathOut']],
                        verbose = verbose)
    
  } else {
    stop(sprintf('Pipeline %s not found.', pipeline))
  }
  write.csv(x = out, file = file.path(outDir, 'Results.csv'), row.names = F, col.names = T)
}

quantPipeline = function(pathToDataFile, nomenclFile, calCurveFile, normBy = 'TSP', remLst = NULL, vAliquot, vNMRsample, concRef=2.5, concRefSC = 2.5, exclude = NULL, verbose = F){
  
  if(verbose) message(sprintf('Reading data from %s', pathToDataFile))
  data <- read.csv(pathToDataFile, header=T)
  data <- parseSpinSysData(data)
  data <- changeSpins2Names(data, nomenclFile)
  if(verbose) message(sprintf('Normalising by %s. Using concRef=%.2f concrefSC=%.2f', normBy, concRef, concRefSC))
  data <- normaliseByPeak(data, normBy = normBy, exclude = exclude, concRef = concRef, concRefSC = concRefSC)
  data <- changeNomenclature(data, nomenclatureFile, excludeNotFound=T)
  
  if(!is.null(remLst)) data <- RemovePeaksByName(data, remLst)
  
  data <- setNamesAndOrdering(data)
  if(verbose) message(sprintf('Calculating concentrations.'))
  data <- calcConc(data, calCurveFile)
  if(verbose) message(sprintf('Adjusting to aliquot volume: V_NMRsample=%d V_aliquot=%d', vNMRsample, vAliquot))
  data <- adjustAliquot(data, vNMRsample, vAliquot)
  data
}

runCalCurvePipeline <- function(pathToData, nomenclFile, normBy='none', exclude='none', saveAsPdf=F, saveTo=NULL, verbose=T){
  if(verbose) message(sprintf('Reading data from %s', pathToData))
  data <- readFolder(pathToData)
  data <- changeSpins2Names(data, nomenclFile)
  if(verbose) message(sprintf('Plotting'))
  out <- plotNMRCalCurves(data, normBy = normBy, exclude = exclude, saveAsPdf = saveAsPdf, saveTo = saveTo)
  out
}

runQuantPipeline <- function(path, nomenclFile, carb2metFile, calCurveFile, vAliquot, yHeight, remLst = NULL, exclude = "HEPES", title, plotBars = F, saveToFile=T, verbose=T){
  
  files <- list.files(path)
  processed <- lapply(1:length(files), function(i){
    filePath <- file.path(path, files[i])
    result <- quantPipeline(filePath, nomenclFile, calCurveFile, remLst, vAliquot, vNMRsample, exclude=NULL)
    result <- pickOnePeak(result, carb2MetFile)
    result$Condition <- i
    colnames(result) <- c("Metabolite","Height","Volume","Condition")
  })
  
  resultDF <- do.call(rbind, processed)
  resultDF$Condition <- as.factor(resultDF$Condition)
  if(!is.null(exclude)) resultDF <- resultDF[!(resultDF$Metabolite %in% exclude),]

  if(plotBars) plotBarsWithErrors(resultDF, title, yheight)
  if(saveToFile) write.csv(resultDF, paste(title,".csv"))
  resultDF
}
