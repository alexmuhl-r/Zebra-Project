library(ez)
library(ggplot2)
library(plyr)
library(grid)
library(gridExtra)
library(rio)
library(ptinpoly)
library(sp)
library(zoo)
library(readxl)    
library(tidyr)
library(useful)
library(magick)
library(cowplot)
library(data.table)
library("viridis")     
library(lsr)

#Helper functions below - summary stats and trans.arcsine

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

## Norms the data within specified groups in a data frame; it normalizes each
## subject (identified by idvar) so that they have the same mean, within each group
## specified by betweenvars.
##   data: a data frame.
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
  library(plyr)
  
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}

## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same between-group mean),
##   standard deviation, standard error of the mean, and confidence interval.
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))
  
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}

trans.arcsine = function(x){
  asin(sqrt(x))
}

trans.arcsine.col = function(x){
  sapply(x,trans.arcsine)
}

prevCount = function(x) {ave(x==x,x,FUN=cumsum)}

#function for reading multi sheet excel workbooks
read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

#Helper function loading complete

################################################################################

#Load and merge data

allTrialsMousetrackingOrder1 = list.files(pattern= glob2rx("*r58b*mousetracking_trial-*.xlsx$*"), recursive = TRUE)
allTrialsMousetrackingOrder2 = list.files(pattern= glob2rx("*uitz*mousetracking_trial-*.xlsx$*"), recursive = TRUE)

allTrialsMousetrackingListOrder1 = list()
allTrialsMousetrackingListOrder2 = list()

for (n in 1:length(allTrialsMousetrackingOrder1))
{allTrialsMousetrackingListOrder1[[n]] = import(allTrialsMousetrackingOrder1[n])
allTrialsMousetrackingListOrder1[[n]]$mousedataFilename = allTrialsMousetrackingOrder1[n]}
for (n in 1:length(allTrialsMousetrackingOrder2))
{allTrialsMousetrackingListOrder2[[n]] = import(allTrialsMousetrackingOrder2[n])
allTrialsMousetrackingListOrder2[[n]]$mousedataFilename = allTrialsMousetrackingOrder2[n]}

gorillaDataOrder1= import("data_exp_34079-v13_task-r58b.csv")
gorillaDataOrder1$spreadsheet_row = gorillaDataOrder1$`Spreadsheet Row`

gorillaDataOrder2 = import("data_exp_34079-v13_task-uitz.csv")
gorillaDataOrder2$spreadsheet_row = gorillaDataOrder2$`Spreadsheet Row`

singleTrialDataOrder1 = gorillaDataOrder1[gorillaDataOrder1$display == "StimDisplay",]
singleTrialDataOrder1$participant_id = singleTrialDataOrder1$`Participant Private ID`
singleTrialDataOrder1 = singleTrialDataOrder1[!duplicated(singleTrialDataOrder1[c("spreadsheet_row","participant_id")]),]

singleTrialDataOrder2 = gorillaDataOrder2[gorillaDataOrder2$display == "StimDisplay",]
singleTrialDataOrder2$participant_id = singleTrialDataOrder2$`Participant Private ID`
singleTrialDataOrder2 = singleTrialDataOrder2[!duplicated(singleTrialDataOrder2[c("spreadsheet_row","participant_id")]),]

allTrialsMousetrackingDataOrder1 = ldply(allTrialsMousetrackingListOrder1, data.frame)
allTrialsMousetrackingDataOrder2 = ldply(allTrialsMousetrackingListOrder2, data.frame)

mergedDataOrder1 = merge(allTrialsMousetrackingDataOrder1,singleTrialDataOrder1, by = c("spreadsheet_row","participant_id"), all.y= T)
mergedDataOrder1 = mergedDataOrder1[order(mergedDataOrder1$participant_id,mergedDataOrder1$time_stamp),]
mergedDataOrder1$task = 'Order1'

mergedDataOrder2 = merge(allTrialsMousetrackingDataOrder2,singleTrialDataOrder2, by = c("spreadsheet_row","participant_id"), all.y = T)
mergedDataOrder2 = mergedDataOrder2[order(mergedDataOrder2$participant_id,mergedDataOrder2$time_stamp),]
mergedDataOrder2$task = 'Order2'

mergedData = rbind.fill(mergedDataOrder1,mergedDataOrder2)

count(mergedData,c("task","participant_id"))

mouseMovementsOnly = mergedData[mergedData$display == "StimDisplay",]
mouseMovementsOnly = mouseMovementsOnly[mouseMovementsOnly$type=="mouse",]
mouseMovementsOnly = mouseMovementsOnly[!is.na(mouseMovementsOnly$task),]

participant_ids = mergedData[!duplicated(mergedData$participant_id),]$participant_id

mouseMovementsOnly$`Trial Number` = as.numeric(mouseMovementsOnly$`Trial Number`)
mouseMovementsOnly = mouseMovementsOnly[!is.na(mouseMovementsOnly$`Trial Number`),]

#count(mouseMovementsOnly[!(duplicated(mouseMovementsOnly[c('participant_id','Trial Number')])),],c("task","participant_id"))
aggregate(`Trial Number` ~ task + participant_id, mouseMovementsOnly, FUN = min)

#optional - write data to files
#write.csv(mergedData,"mergedData.csv")
#write.csv(mouseMovementsOnly,"mouseMovementsOnly.csv")

participant_ids = mergedData[!duplicated(mergedData$participant_id),]$participant_id

#make a blank dataframe to add mouse data to
reconstructedMouseDataframe = data.frame(0)

mouseMovementsOnly$Trial.Number = mouseMovementsOnly$`Trial Number`
mergedData$Trial.Number = mergedData$`Trial Number`

screenPxData = data.frame(ppt = NA, screenX=NA,screenY=NA,imageN_xmin=NA,imageN_xmax=NA,imageN_ymin=NA,imageN_ymax=NA)

for (ppt in participant_ids){
  mouseMovementsOnlyPpt = mouseMovementsOnly[mouseMovementsOnly$participant_id==ppt,]
  
  mergedDataPpt = mergedData[mergedData$participant_id==ppt,]
  trialNumber = 1
  trial = mouseMovementsOnlyPpt[mouseMovementsOnlyPpt$Trial.Number==trialNumber,]
  mergedSingleTrial = mergedDataPpt[mergedDataPpt$Trial.Number==trialNumber,]
  
  screenX = mergedSingleTrial[mergedSingleTrial$zone_name=="screen",][!is.na(mergedSingleTrial[mergedSingleTrial$zone_name=="screen",]$zone_name),][1,]$zone_width
  screenY = mergedSingleTrial[mergedSingleTrial$zone_name=="screen",][!is.na(mergedSingleTrial[mergedSingleTrial$zone_name=="screen",]$zone_name),][1,]$zone_height
  
  while (nrow(trial)==0){
    trialNumber = trialNumber + 1
    trial = mouseMovementsOnlyPpt[mouseMovementsOnlyPpt$Trial.Number==trialNumber,]
    if (trialNumber > 192){
      break
    }
  }
  
  imageN = mergedDataPpt[mergedDataPpt$zone_name=="ImageN",][!is.na(mergedDataPpt[mergedDataPpt$zone_name=="ImageN",]$zone_name),][1,]
  
  imageN_xmin = imageN$zone_x
  imageN_xmax = imageN$zone_x + imageN$zone_width
  imageN_ymin = imageN$zone_y
  imageN_ymax = imageN$zone_y + imageN$zone_height
  
  pptScreenPxData = c(ppt,screenX,screenY,imageN_xmin,imageN_xmax,imageN_ymin,imageN_ymax)
  screenPxData = rbind(screenPxData,pptScreenPxData)
  
  for (t in 1:192){
    if (t %in% mouseMovementsOnlyPpt$Trial.Number){
      trialLoop = mouseMovementsOnlyPpt[mouseMovementsOnlyPpt$Trial.Number==t,]
      trialLoop$screenX = screenX
      trialLoop$screenY = screenY
      
      preTrialFixCrossLoop = trialLoop[trialLoop$filename=='mousetracking_pretrial_fixcross',]
      mousetrackingTrialLoop = trialLoop[trialLoop$filename=='mousetracking_trial',]
      
      reconstructedMouseDataframe = rbind.fill(reconstructedMouseDataframe,mousetrackingTrialLoop)
    }
  }
}

screenPxData = screenPxData[!is.na(screenPxData$ppt),]

#optional write data to file
#write.csv(screenPxData,"screenData.csv")

reconstructedMouseDataframe = reconstructedMouseDataframe[!is.na(reconstructedMouseDataframe$spreadsheet_row),]
reconstructedMouseDataframe$flipped_x = reconstructedMouseDataframe$x
reconstructedMouseDataframe$flipped_x_normalised = reconstructedMouseDataframe$x_normalised

#optional write data to file
#write.csv(reconstructedMouseDataframe,"reconstructedMouseData.csv")

#Select all left facing images and flip the x-coords to standardise direction
reconstructedMouseDataframe[reconstructedMouseDataframe$ImageN %like% "_L_", ]$flipped_x = 
  reconstructedMouseDataframe[reconstructedMouseDataframe$ImageN %like% "_L_", ]$screenX - 
  reconstructedMouseDataframe[reconstructedMouseDataframe$ImageN %like% "_L_", ]$x
reconstructedMouseDataframe[reconstructedMouseDataframe$ImageN %like% "_L_", ]$flipped_x_normalised = 
  1 - reconstructedMouseDataframe[reconstructedMouseDataframe$ImageN %like% "_L_", ]$x_normalised

#Take mouse sample with latest time stamp
reconstructedMouseDataframe = reconstructedMouseDataframe[order(reconstructedMouseDataframe$participant_id, reconstructedMouseDataframe$time_stamp),]
lastMouseCoordsDataFrame = reconstructedMouseDataframe[!duplicated(reconstructedMouseDataframe[, c("participant_id","Trial.Number")], fromLast=T),]

count(lastMouseCoordsDataFrame,c("participant_id"))

#optional write data to file
#write.csv(lastMouseCoordsDataFrame,"lastMouseCoordsData.csv")

#data loading complete

##################################################################################

#check position of stimuli on screen using mouse data (to account for variations in screen size

allPosCheck = list.files(pattern= glob2rx("*mousetracking_stimcheck*.xlsx$*"), recursive = TRUE)

allPosCheckList = list()

for (n in 1:length(allPosCheck))
{allPosCheckList[[n]] = import(allPosCheck[n])
allPosCheckList[[n]]$mousedataFilename = allPosCheck[n]}

allPosCheckDF = ldply(allPosCheckList, data.frame)

#exclude some participants for incomplete data
allPosCheckDF = subset(allPosCheckDF, !(allPosCheckDF$participant_id %in% c(2460441,2460450,2460483,2474755)))

allPosCheckDF = allPosCheckDF[allPosCheckDF$type=="mouse",]

allPosCheckDF = allPosCheckDF[order(allPosCheckDF$participant_id, allPosCheckDF$time_stamp),]
allPosCheckDF = allPosCheckDF[!duplicated(allPosCheckDF[, c("participant_id","spreadsheet_row")], fromLast=T),]

allPosCheckFirstDF = allPosCheckDF[allPosCheckDF$spreadsheet_row < 205,]
allPosCheckLastDF = allPosCheckDF[allPosCheckDF$spreadsheet_row > 205,]

allPosCheckFirstDF = allPosCheckFirstDF[c('participant_id','x','y','x_normalised','y_normalised')]
colnames(allPosCheckFirstDF) = c('participant_id','x_1','y_1','x_normalised_1','y_normalised_1')
allPosCheckLastDF = allPosCheckLastDF[c('participant_id','x','y','x_normalised','y_normalised')]
colnames(allPosCheckLastDF) = c('participant_id','x_2','y_2','x_normalised_2','y_normalised_2')
allPosCheckLastDF = rbind(allPosCheckLastDF,c('2460445',NA,NA,NA,NA))
allPosCheckFirstDF = allPosCheckFirstDF[order(allPosCheckFirstDF$participant_id),]
allPosCheckLastDF = allPosCheckLastDF[order(allPosCheckLastDF$participant_id),]

allPosCheckMeansDF = cbind(allPosCheckFirstDF,allPosCheckLastDF)
allPosCheckMeansDF$x_2 = as.numeric(allPosCheckMeansDF$x_2)
allPosCheckMeansDF$y_2 = as.numeric(allPosCheckMeansDF$y_2)
allPosCheckMeansDF$x_normalised_2 = as.numeric(allPosCheckMeansDF$x_normalised_2)
allPosCheckMeansDF$y_normalised_2 = as.numeric(allPosCheckMeansDF$y_normalised_2)

allPosCheckMeansDF$x_mean = rowMeans(subset(allPosCheckMeansDF, select = c('x_1', 'x_2')), na.rm = TRUE)
allPosCheckMeansDF$x_mean_norm = rowMeans(subset(allPosCheckMeansDF, select = c('x_normalised_1', 'x_normalised_2')), na.rm = TRUE)
allPosCheckMeansDF$y_mean = rowMeans(subset(allPosCheckMeansDF, select = c('y_1', 'y_2')), na.rm = TRUE)
allPosCheckMeansDF$y_mean_norm = rowMeans(subset(allPosCheckMeansDF, select = c('y_normalised_1', 'y_normalised_2')), na.rm = TRUE)

allPosCheckMeansDF = allPosCheckMeansDF[c('participant_id','x_mean','x_mean_norm','y_mean','y_mean_norm')]
allPosCheckMeansDF$coordNum = rep(1:4,24)
allPosCheckMeansDF = reshape(allPosCheckMeansDF,timevar = "coordNum",idvar = c("participant_id"),direction = "wide")

fullAveragedPosCheckDF = allPosCheckMeansDF
fullAveragedPosCheckDF$x_mean.1 = rowMeans(subset(fullAveragedPosCheckDF, select = c('x_mean.1', 'x_mean.4')), na.rm = TRUE)
fullAveragedPosCheckDF$x_mean_norm.1 = rowMeans(subset(allPosCheckMeansDF, select = c('x_mean_norm.1', 'x_mean_norm.4')), na.rm = TRUE)
fullAveragedPosCheckDF$x_mean.4 = fullAveragedPosCheckDF$x_mean.1
fullAveragedPosCheckDF$x_mean_norm.4 = fullAveragedPosCheckDF$x_mean_norm.1
fullAveragedPosCheckDF$x_mean.2 = rowMeans(subset(fullAveragedPosCheckDF, select = c('x_mean.2', 'x_mean.3')), na.rm = TRUE)
fullAveragedPosCheckDF$x_mean_norm.2 = rowMeans(subset(allPosCheckMeansDF, select = c('x_mean_norm.2', 'x_mean_norm.3')), na.rm = TRUE)
fullAveragedPosCheckDF$x_mean.3 = fullAveragedPosCheckDF$x_mean.2
fullAveragedPosCheckDF$x_mean_norm.3 = fullAveragedPosCheckDF$x_mean_norm.2
fullAveragedPosCheckDF$y_mean.1 = rowMeans(subset(fullAveragedPosCheckDF, select = c('y_mean.1', 'y_mean.2')), na.rm = TRUE)
fullAveragedPosCheckDF$y_mean_norm.1 = rowMeans(subset(allPosCheckMeansDF, select = c('y_mean_norm.1', 'y_mean_norm.2')), na.rm = TRUE)
fullAveragedPosCheckDF$y_mean.2 = fullAveragedPosCheckDF$y_mean.1
fullAveragedPosCheckDF$y_mean_norm.2 = fullAveragedPosCheckDF$y_mean_norm.1
fullAveragedPosCheckDF$y_mean.3 = rowMeans(subset(fullAveragedPosCheckDF, select = c('y_mean.3', 'y_mean.4')), na.rm = TRUE)
fullAveragedPosCheckDF$y_mean_norm.3 = rowMeans(subset(allPosCheckMeansDF, select = c('y_mean_norm.3', 'y_mean_norm.4')), na.rm = TRUE)
fullAveragedPosCheckDF$y_mean.3 = fullAveragedPosCheckDF$y_mean.4
fullAveragedPosCheckDF$y_mean_norm.3 = fullAveragedPosCheckDF$y_mean_norm.4

x1Mean = mean(fullAveragedPosCheckDF$x_mean_norm.1)
y1Mean = mean(fullAveragedPosCheckDF$y_mean_norm.1)
x2Mean = mean(fullAveragedPosCheckDF$x_mean_norm.2)
y2Mean = mean(fullAveragedPosCheckDF$y_mean_norm.2)
x3Mean = mean(fullAveragedPosCheckDF$x_mean_norm.3)
y3Mean = mean(fullAveragedPosCheckDF$y_mean_norm.3)
x4Mean = mean(fullAveragedPosCheckDF$x_mean_norm.4)
y4Mean = mean(fullAveragedPosCheckDF$y_mean_norm.4)

#position checking complete

##############################################################################################

#data analysis

# optional checking for left right bias - make sure no one over 74%
# lastMouseCoordsDataFrame$lookRearUnflipped = ifelse(lastMouseCoordsDataFrame$x_normalised < 0.5,1,0)
# freqTablePPT = count(lastMouseCoordsDataFrame,c("participant_id","lookRearUnflipped"))
# freqTablePPT$prop = ldply(by(freqTablePPT['freq'],freqTablePPT$participant_id,FUN=prop.table,simplify=F))$freq

lastMouseCoordsDataFrame$blurred = "N"
lastMouseCoordsDataFrame[lastMouseCoordsDataFrame$ImageN %like% "_B_", ]$blurred = "B"
lastMouseCoordsDataFrame$facing = "R"
lastMouseCoordsDataFrame[lastMouseCoordsDataFrame$ImageN %like% "_L_", ]$facing = "L"
lastMouseCoordsDataFrame$image_directionless = substr(lastMouseCoordsDataFrame$ImageN, 1, 6) 
lastMouseCoordsDataFrame$image_directionless_B <- with(lastMouseCoordsDataFrame, paste0(image_directionless, blurred))
lastMouseCoordsDataFrame$ImageN_flipped = gsub("L", "R", lastMouseCoordsDataFrame$ImageN)

ImageN_flipped = lastMouseCoordsDataFrame[!duplicated(lastMouseCoordsDataFrame$ImageN_flipped),]$ImageN_flipped

# #optional, create black and white plots for correction  with computational salience estimates
# for (stim in ImageN_flipped){
#   
#   stimPlotImgOnly = ggplot(lastMouseCoordsDataFrame[lastMouseCoordsDataFrame$ImageN_flipped == stim,], aes(x = flipped_x_normalised, y = y_normalised)) + 
#     stat_density2d(geom = "tile", aes(fill = ..density.., alpha=..density..), contour = FALSE) + 
#     xlim(x1Mean,x2Mean) +
#     ylim(y3Mean,y2Mean) + 
#     scale_fill_gradient(low="black", high="white") + #theme_nothing() + labs(x = NULL, y = NULL) 
#     theme(axis.line=element_blank(),axis.text.x=element_blank(),
#           axis.text.y=element_blank(),axis.ticks=element_blank(),
#           axis.title.x=element_blank(),
#           axis.title.y=element_blank(),legend.position="none",
#           panel.background=element_rect(fill="black"), plot.background = element_rect(fill = "black"), panel.border=element_blank(),panel.grid.major=element_blank(),
#           panel.grid.minor=element_blank())
#   
#   png(paste("blackAndWhite_stim_",stim,".png", sep=''),width=512,height=512)
#   #scale up the plot to cover white border
#   print(plot_grid(stimPlotImgOnly, scale=1.15))
#   dev.off()
# }

lastMouseCoordsDataFrame$species = substr(lastMouseCoordsDataFrame$ImageN, 4, 4) 
lastMouseCoordsDataFrame$blurred = as.factor(lastMouseCoordsDataFrame$blurred)
lastMouseCoordsDataFrame$species = as.factor(lastMouseCoordsDataFrame$species)
lastMouseCoordsDataFrame$facing = as.factor(lastMouseCoordsDataFrame$facing)

lastMouseCoordsDataFrame$lookRear = ifelse(lastMouseCoordsDataFrame$flipped_x_normalised < 0.5,1,0)

lastMouseCoordsDataFrame$lookRear = as.factor(lastMouseCoordsDataFrame$lookRear)
lastMouseCoordsDataFrame$participant_id = as.factor(lastMouseCoordsDataFrame$participant_id)

lastMouseCoordsDataFrame = subset(lastMouseCoordsDataFrame, !(lastMouseCoordsDataFrame$participant_id %in% c(2460441,2460450,2460483,2474755)))

freqTable = count(lastMouseCoordsDataFrame,c("participant_id","species","lookRear","blurred","facing"))

#optional write data 
# write.csv(freqTable,"temporaryFreqTable.csv")

#make blank data
speciesVector = c(rep("G", 8),rep("M",8),rep("S",8),rep("N",8))
blurredVector = rep(c(rep("B",2),rep("N",2)),8)
facingVector = rep(c("L","R"),16)
lookRearVector = rep(c(0,0,0,0,1,1,1,1),4)

participantVector = sort(rep(as.numeric(levels(freqTable[!duplicated(lastMouseCoordsDataFrame$participant_id),]$participant_id)),32))

placeholderDF = data.frame("participant_id" = participantVector, "species" = rep(speciesVector,27), 
                           "blurred" = rep(blurredVector,27), "facing" = rep(facingVector,27), 
                           "lookRear" = rep(lookRearVector,27), "freq" = NA)
placeholderDF = subset(placeholderDF, !(placeholderDF$participant_id %in% c(2460441,2460450,2460483,2474755)))

filledDF = merge(freqTable,placeholderDF,by = c("participant_id","species","lookRear","blurred","facing"),all.y=T)
filledDF = filledDF[,1:6]
colnames(filledDF) = c("participant_id","species","lookRear","blurred","facing","freq")
filledDF[is.na(filledDF$freq),]$freq = 0

library(dplyr)
freqTableDF = data.frame(filledDF %>%
                           group_by(participant_id, species, facing, blurred) %>%
                           mutate(prop = freq/sum(freq)))
detach("package:useful", unload=TRUE)
detach("package:dplyr", unload=TRUE)

freqTableDFRear = freqTableDF[freqTableDF$lookRear==1,]
freqTableDFRear$prop = freqTableDFRear$prop - 0.5

freqTableDFRear$blurred = as.factor(freqTableDFRear$blurred)
freqTableDFRear$species = as.factor(freqTableDFRear$species)
freqTableDFRear$facing = as.factor(freqTableDFRear$facing)
freqTableDFRear$participant_id = as.factor(freqTableDFRear$participant_id)

propSummaryLookRear = summarySEwithin(freqTableDFRear, measurevar="prop", withinvars=c("blurred","facing","species"),
                                      idvar="participant_id", na.rm=FALSE, conf.interval=.95)

summarySEwithin(freqTableDFRear, measurevar="prop", withinvars=c("facing", "species"),
                idvar="participant_id", na.rm=FALSE, conf.interval=.95)

#statistical analysis below

ezANOVA(freqTableDFRear, prop, participant_id, within=.(blurred,facing,species), detailed=T, type = 3)

freqTableDFRearG = freqTableDFRear[freqTableDFRear$species=="G",]
freqTableDFRearM = freqTableDFRear[freqTableDFRear$species=="M",]
freqTableDFRearS = freqTableDFRear[freqTableDFRear$species=="S",]
freqTableDFRearN = freqTableDFRear[freqTableDFRear$species=="N",]

ezANOVA(freqTableDFRearG, prop, participant_id, within=.(blurred,facing), detailed=T, type = 3)
ezANOVA(freqTableDFRearM, prop, participant_id, within=.(blurred,facing), detailed=T, type = 3)
ezANOVA(freqTableDFRearS, prop, participant_id, within=.(blurred,facing), detailed=T, type = 3)
ezANOVA(freqTableDFRearN, prop, participant_id, within=.(blurred,facing), detailed=T, type = 3)

freqTableDFRearGvM = freqTableDFRear[freqTableDFRear$species %in% c("G","M"),]
freqTableDFRearGvN = freqTableDFRear[freqTableDFRear$species %in% c("G","N"),]

ezANOVA(freqTableDFRearGvM, prop, participant_id, within=.(blurred,facing,species), detailed=T, type = 3)
ezANOVA(freqTableDFRearGvN, prop, participant_id, within=.(blurred,facing,species), detailed=T, type = 3)

freqTableDFRearG_B = freqTableDFRearG[freqTableDFRearG$blurred == "N",]
freqTableDFRearG_B_NOFACE = aggregate(prop ~ participant_id + species + lookRear + blurred, freqTableDFRearG_B, mean)
freqTableDFRearN_B = freqTableDFRearN[freqTableDFRearN$blurred == "N",]
freqTableDFRearN_B_NOFACE = aggregate(prop ~ participant_id + species + lookRear + blurred, freqTableDFRearN_B, mean)
freqTableDFRearM_B = freqTableDFRearM[freqTableDFRearM$blurred == "N",]
freqTableDFRearM_B_NOFACE = aggregate(prop ~ participant_id + species + lookRear + blurred, freqTableDFRearM_B, mean)

t.test(freqTableDFRearG_B_NOFACE$prop, mu = 0, alternative = "greater")
cohensD(freqTableDFRearG_B_NOFACE$prop, mu = 0)
t.test(freqTableDFRearN_B_NOFACE$prop, mu = 0, alternative = "greater")
cohensD(freqTableDFRearN_B_NOFACE$prop, mu = 0)
t.test(freqTableDFRearM_B_NOFACE$prop, mu = 0, alternative = "greater")
cohensD(freqTableDFRearM_B_NOFACE$prop, mu = 0)
