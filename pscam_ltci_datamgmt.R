library(dplyr)

## -- Download data from REDCap --------------------------------------------------------------------
library(RCurl)

## Read in API token from a supersecret file
source('api_info.R')

if(Sys.info()['sysname'] == 'Darwin'){
  load('/Volumes/thomps23/ICUDelirium/PsCAM/RiskFactors/pscam_risk_datasets.Rdata')
} else{
  load('/home/thomps23/ICUDelirium/PsCAM/RiskFactors/pscam_risk_datasets.Rdata')
}

pscam.ltci <- read.csv(text = postForm(uri = 'https://redcap.vanderbilt.edu/api/',
                                       token = my.pedsltci.token,
                                       content = 'record',
                                       format = 'csv',
                                       type = 'flat',
                                       rawOrLabel = 'label',
                                       rawOrLabelHeaders = 'raw',
                                       returnFormat = 'csv'),
                       stringsAsFactors = TRUE, na.strings = c('', 'NULL'))

names(pscam.ltci) <- gsub('_', '.', names(pscam.ltci))

## -- Create indicators for whether patient is considered delayed on each ASQ domain, ASQ:SE -------
## -- (cutoff differs depending on age)-------------------------------------------------------------

## Matrix of ages, cutoffs for each ASQ domain
asq.cutoffs <- matrix(NA, nrow = length(levels(pscam.ltci$asq)), ncol = 5)
rownames(asq.cutoffs) <- levels(pscam.ltci$asq)
colnames(asq.cutoffs) <- c('com', 'gm', 'fm', 'ps', 'social')

asq.cutoffs[,'com'] <- c(40, 35, rep(40, 11), 35)
asq.cutoffs[,'gm'] <- c(30, 35, 30, 40, 30, 40, 30, 35, 30, 40, 30, 35, 30, 35)
asq.cutoffs[,'fm'] <- c(30, 35, 30, 40, 30, 40, 30, 30, 30, 35, 30, 35, 30, 35)
asq.cutoffs[,'ps'] <- c(30, 30, 30, 35, 30, 35, 30, 30, 30, 40, 30, 40, 30, 35)
asq.cutoffs[,'social'] <- c(30, 30, 30, 40, 30, 40, 30, 40, 30, 40, 30, 25, 30, 40)

## Function to determine whether patient was delayed on a given test, depending on age
delayed.domain <- function(recnum, domain = colnames(asq.cutoffs)){
  domain.score <- paste0(domain, '.score')
  return(ifelse(is.na(pscam.ltci[recnum, domain.score]), NA,
                pscam.ltci[recnum, domain.score] <= asq.cutoffs[pscam.ltci[recnum, 'asq'], domain]))
}

## Wrapper function to create variables for each record
delayed.indicator <- function(domain = colnames(asq.cutoffs)){
  unlist(lapply(1:nrow(pscam.ltci), FUN = delayed.domain, domain = domain))
}

pscam.ltci$com.delayed <- delayed.indicator('com')
pscam.ltci$gm.delayed <- delayed.indicator('gm')
pscam.ltci$fm.delayed <- delayed.indicator('fm')
pscam.ltci$ps.delayed <- delayed.indicator('ps')
pscam.ltci$social.delayed <- delayed.indicator('social')

## Indicator for any developmental delay (delay on any domain)
pscam.ltci$any.asq.delayed <-
  rowSums(pscam.ltci[,paste0(colnames(asq.cutoffs), '.delayed')], na.rm = TRUE) > 0

## Determine whether patient was delayed on ASQ:SE
pscam.ltci$asqse.delayed <- with(pscam.ltci, {
  ifelse(is.na(asqse.total.score), NA,
         (asqse == '12 month' & asqse.total.score > 48) |
           (asqse == '18 month' & asqse.total.score > 50) |
           (asqse == '24 month' & asqse.total.score > 50) |
           (asqse == '30 month' & asqse.total.score > 57) |
           (asqse == '36 month' & asqse.total.score > 59) |
           (asqse == '48 month' & asqse.total.score > 70) |
           (asqse == '60 month' & asqse.total.score > 70)) })

## Indicator for any developmental delay (delay on any domain)
pscam.ltci$any.delayed <-
  rowSums(pscam.ltci[,c('any.asq.delayed', 'asqse.delayed')], na.rm = TRUE) > 0

## -- Calculate in-hospital summary variables ------------------------------------------------------
inhosp.summary <- analysis.data %>%
  group_by(study.id) %>%
  summarise(mean.benzo.wt = mean(benzo.wt, na.rm = TRUE))

## -- Calculate time between discharge and assessment ----------------------------------------------
pscam.ltci <- merge(pscam.ltci, subset(enroll.data, select = c(study.id, discharge.date)),
                    by.x = 'record.id', by.y = 'study.id', all.x = TRUE, all.y = FALSE)
pscam.ltci$asq.date <- as.Date(as.character(pscam.ltci$date), format = '%Y-%m-%d')
pscam.ltci$days.dc.asq <-
  with(pscam.ltci, as.numeric(difftime(asq.date, discharge.date, units = 'days')))

## -- Create final analysis data set ---------------------------------------------------------------
peds.ltci.data <- subset(demo.data.oneobs,
                         select = -days.enroll.died,
                         study.id %in% pscam.ltci$record.id) %>%
  left_join(inhosp.summary, by = 'study.id') %>%
  full_join(select(pscam.ltci, record.id, dev.delay, com.baseline, gm.baseline, fm.baseline,
                   ps.baseline, social.baseline, days.dc.asq, asq, com.score, com.delayed,
                   gm.score, gm.delayed, fm.score, fm.delayed, ps.score, ps.delayed, social.score,
                   social.delayed, any.asq.delayed, asqse, asqse.total.score, asqse.delayed,
                   any.delayed),
            by = c('study.id' = 'record.id'))

peds.ltci.data <- as.data.frame(peds.ltci.data)

label(peds.ltci.data$study.id) <- 'Hospitalization ID'
label(peds.ltci.data$mean.benzo.wt) <- 'Mean 24h benzodiazepines in hospital (mg/kg, midaz equiv.)'
label(peds.ltci.data$dev.delay) <- 'Developmental delay at baseline, per PCP records'
label(peds.ltci.data$com.baseline) <- 'Communication delay at baseline'
label(peds.ltci.data$gm.baseline) <- 'Gross motor delay at baseline'
label(peds.ltci.data$fm.baseline) <- 'Fine motor delay at baseline'
label(peds.ltci.data$ps.baseline) <- 'Problem solving delay at baseline'
label(peds.ltci.data$social.baseline) <- 'Personal/social delay at baseline'
label(peds.ltci.data$days.dc.asq) <- 'Days between hospital discharge and ASQ assessment'
label(peds.ltci.data$asq) <- 'ASQ version'
label(peds.ltci.data$com.score) <- 'ASQ communication score'
label(peds.ltci.data$com.delayed) <- 'Delayed on communication domain'
label(peds.ltci.data$gm.score) <- 'ASQ gross motor score'
label(peds.ltci.data$gm.delayed) <- 'Delayed on gross motor domain'
label(peds.ltci.data$fm.score) <- 'ASQ fine motor score'
label(peds.ltci.data$fm.delayed) <- 'Delayed on fine motor domain'
label(peds.ltci.data$ps.score) <- 'ASQ problem solving score'
label(peds.ltci.data$ps.delayed) <- 'Delayed on problem solving domain'
label(peds.ltci.data$social.score) <- 'ASQ personal-social score'
label(peds.ltci.data$social.delayed) <- 'Delayed on personal-social domain'
label(peds.ltci.data$any.asq.delayed) <- 'Delayed on any ASQ domain'
label(peds.ltci.data$asqse) <- 'ASQ:SE version'
label(peds.ltci.data$asqse.total.score) <- 'ASQ:SE score'
label(peds.ltci.data$asqse.delayed) <- 'Delayed on ASQ:SE'
label(peds.ltci.data$any.delayed) <- 'Delayed on any ASQ domain or ASQ:SE'

## Save date that analysis data sets were created
dataset.created.at <- Sys.time()

save(peds.ltci.data, dataset.created.at, file = 'pscam_ltci.Rdata')
