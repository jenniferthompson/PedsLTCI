## -- Download TBI data by event -------------------------------------------------------------------
library(RCurl)

## Read in API token from a supersecret file
source('api_info.R')

## Set constant variables
use.uri <- 'https://redcap.vanderbilt.edu/api/'
use.format <- 'csv'
use.rawlabel <- 'label'
use.rawheader <- 'raw'

pscam.ltci <- read.csv(text = postForm(uri = use.uri,
                                       token = my.pedsltci.token,
                                       content = 'record',
                                       format = use.format,
                                       type = 'flat',
                                       rawOrLabel = use.rawlabel,
                                       rawOrLabelHeaders = use.rawheader,
                                       returnFormat = use.format),
                       stringsAsFactors = FALSE, na.strings = c('', 'NULL'))

names(pscam.ltci) <- gsub('_', '.', names(pscam.ltci))

save(pscam.ltci, file = 'pscam_ltci.Rdata')
