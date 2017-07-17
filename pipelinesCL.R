

# A script for calling the pipelines from the command line through Rscript.
source('main.R')
source('lib.R')
# Parse command line parameters
args <- commandArgs(TRUE)

if (length(args) < 1) {
  args <- c("--help")
}

currDir < getwd()
if(!('main.R' %in% list.files(currDir)) | !('lib.R' %in% list.files(currDir)) ) {
  stop('Please run the script from the NMRquantRT folder.')}

if("--help" %in% args) {
  cat("
      pipelinesCL - a script for calling NMRquantTR pipelines

      Arguments:
      --data=path - path to data
      --config=path - pathto config file
      --pipeline=quant - which pipeline to call {quant/makeCC}
      --verbose=TRUE - (optional) show outputs in the command line? (TRUE/FALSE) (default=TRUE)

      Example:
      Rscript pipelinesCL.R --data=data/data_file_1 --config=config.txt --pipeline=quant --verbose=TRUE \n\n")
  q(save="no")
}

parsed <- function(x) strsplit(sub("^--", "", x),"=")
argsDF <- as.data.frame(do.call('rbind', parsed(args)))
args <- as.list(as.character(argsDF[,2]))
names(args) <- argsDF[,1]

if('data' %in% names(args)) {
  data <- args['data']
} else { stop('Path to data not specified.')}

if('config' %in% names(args)) {
  conf <- args['config']
} else { stop('Path to config not specified.')}

if('pipeline' %in% names(args)) {
  pipe <- args['pipeline']
} else { stop('Pipeline not specified.')}

if('verbose' %in% names(args)) {
  verbose <- ifelse(args['verbose']=='FALSE',FALSE, TRUE)
} else { verbose <- TRUE}

# run the pipeline
MainPipeline(pathToData = data, pathToConfig = config, pipeline = pipe, verbose = verbose)