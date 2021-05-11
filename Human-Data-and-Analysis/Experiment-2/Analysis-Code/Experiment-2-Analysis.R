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
library(ggplot2)
library(tidyr)
library(useful)
library(magick)
library(cowplot)
library(data.table)
library("viridis")     
library(stringi)
library(stringr)
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

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
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

gorillaDataOrder1= import("data_exp_36126-v13_task-r58b.csv")
gorillaDataOrder1$spreadsheet_row = gorillaDataOrder1$`Spreadsheet Row`

gorillaDataOrder2 = import("data_exp_36126-v13_task-uitz.csv")
gorillaDataOrder2$spreadsheet_row = gorillaDataOrder2$`Spreadsheet Row`

singleTrialDataOrder1 = gorillaDataOrder1[gorillaDataOrder1$display == "mediumStimDisplay",]
singleTrialDataOrder1$participant_id = singleTrialDataOrder1$`Participant Private ID`
singleTrialDataOrder1 = singleTrialDataOrder1[!duplicated(singleTrialDataOrder1[c("spreadsheet_row","participant_id")]),]

singleTrialDataOrder2 = gorillaDataOrder2[gorillaDataOrder2$display == "mediumStimDisplay",]
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

mouseMovementsOnly = mergedData[mergedData$display == "mediumStimDisplay",]
mouseMovementsOnly = mouseMovementsOnly[mouseMovementsOnly$type=="mouse",]
mouseMovementsOnly = mouseMovementsOnly[!is.na(mouseMovementsOnly$task),]

participant_ids = mergedData[!duplicated(mergedData$participant_id),]$participant_id

mouseMovementsOnly$`Trial Number` = as.numeric(mouseMovementsOnly$`Trial Number`)
mouseMovementsOnly = mouseMovementsOnly[!is.na(mouseMovementsOnly$`Trial Number`),]

count(mouseMovementsOnly,c("task","participant_id"))
aggregate(`Trial Number` ~ task + participant_id, mouseMovementsOnly, FUN = min)

#optional write data to files
#write.csv(mergedData,"mergedData.csv")
#write.csv(mouseMovementsOnly,"mouseMovementsOnly.csv")

participant_ids = mergedData[!duplicated(mergedData$participant_id),]$participant_id
#crossImg <- as.raster(image_fill(image_read('cross.png'), 'none'))

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
    if (trialNumber > 288){
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
  
  for (t in 1:288){
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
allPosCheckDF = subset(allPosCheckDF, !(allPosCheckDF$participant_id %in% c(2639221,2639155)))

allPosCheckDF = allPosCheckDF[allPosCheckDF$type=="mouse",]

allPosCheckDF = allPosCheckDF[order(allPosCheckDF$participant_id, allPosCheckDF$time_stamp),]
allPosCheckDF = allPosCheckDF[!duplicated(allPosCheckDF[, c("participant_id","spreadsheet_row")], fromLast=T),]

allPosCheckFirstDF = allPosCheckDF[allPosCheckDF$spreadsheet_row < 205,]
allPosCheckLastDF = allPosCheckDF[allPosCheckDF$spreadsheet_row > 205,]

allPosCheckFirstDF = allPosCheckFirstDF[c('participant_id','x','y','x_normalised','y_normalised')]
colnames(allPosCheckFirstDF) = c('participant_id','x_1','y_1','x_normalised_1','y_normalised_1')
allPosCheckLastDF = allPosCheckLastDF[c('participant_id','x','y','x_normalised','y_normalised')]
colnames(allPosCheckLastDF) = c('participant_id','x_2','y_2','x_normalised_2','y_normalised_2')
#allPosCheckLastDF = rbind(allPosCheckLastDF,c('2460445',NA,NA,NA,NA))
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
allPosCheckMeansDF$coordNum = rep(1:4,28)
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

# Checking for left right bias - no one over 74%
# lastMouseCoordsDataFrame$lookRearUnflipped = ifelse(lastMouseCoordsDataFrame$x_normalised < 0.5,1,0)
# freqTablePPT = count(lastMouseCoordsDataFrame,c("participant_id","lookRearUnflipped"))
# freqTablePPT$prop = ldply(by(freqTablePPT['freq'],freqTablePPT$participant_id,FUN=prop.table,simplify=F))$freq

lastMouseCoordsDataFrame$blurred = "AC"
lastMouseCoordsDataFrame[lastMouseCoordsDataFrame$ImageN %like% "_MED_", ]$blurred = "MED"
lastMouseCoordsDataFrame[lastMouseCoordsDataFrame$ImageN %like% "_LOW_", ]$blurred = "LOW"
lastMouseCoordsDataFrame$facing = "R"
lastMouseCoordsDataFrame[lastMouseCoordsDataFrame$ImageN %like% "_L_", ]$facing = "L"

lastMouseCoordsDataFrame$image_directionless = sub("\\L_.*", "\\1", lastMouseCoordsDataFrame$ImageN)
lastMouseCoordsDataFrame$image_directionless = sub("\\R_.*", "\\1", lastMouseCoordsDataFrame$image_directionless)

lastMouseCoordsDataFrame$stimulus = str_sub(lastMouseCoordsDataFrame$image_directionless,-4,-1)
lastMouseCoordsDataFrame$stimulus = sub("_", "", lastMouseCoordsDataFrame$stimulus)

lastMouseCoordsDataFrame$size = NA
lastMouseCoordsDataFrame[lastMouseCoordsDataFrame$ImageN %like% "165", ]$size = 165
lastMouseCoordsDataFrame[lastMouseCoordsDataFrame$ImageN %like% "256", ]$size = 256

lastMouseCoordsDataFrame$ImageN_flipped = gsub("_L_", "_R_", lastMouseCoordsDataFrame$ImageN)

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
#   png(paste("blackAndWhite_stim_",stim,".png", sep=''),width=256,height=256)
#   #scale up the plot to cover white border
#   print(plot_grid(stimPlotImgOnly, scale=1.15))
#   dev.off()
# }

lastMouseCoordsDataFrame$species = substr(lastMouseCoordsDataFrame$stimulus, 1, 1) 
lastMouseCoordsDataFrame$blurred = as.factor(lastMouseCoordsDataFrame$blurred)
lastMouseCoordsDataFrame$species = as.factor(lastMouseCoordsDataFrame$species)
lastMouseCoordsDataFrame$facing = as.factor(lastMouseCoordsDataFrame$facing)
lastMouseCoordsDataFrame$size = as.factor(lastMouseCoordsDataFrame$size)

lastMouseCoordsDataFrame$lookRear = ifelse(lastMouseCoordsDataFrame$flipped_x_normalised < 0.5,1,0)

lastMouseCoordsDataFrame$lookRear = as.factor(lastMouseCoordsDataFrame$lookRear)
lastMouseCoordsDataFrame$participant_id = as.factor(lastMouseCoordsDataFrame$participant_id)

#exclusions
lastMouseCoordsDataFrame = subset(lastMouseCoordsDataFrame, !(lastMouseCoordsDataFrame$participant_id %in% c(2639221,2639155)))

##############################################################################

#analysis collapsed across facing direction

freqTableNF = count(lastMouseCoordsDataFrame,c("participant_id","species","lookRear","blurred","size"))

#make blank data
speciesVectorNF = c(rep("G", 12),rep("M",12),rep("S",12),rep("N",12))
blurredVectorNF = rep(c(rep("AC",4),rep("LOW",4),rep("MED",4)),4)
sizeVectorNF = rep(c(165,256),24)
lookRearVectorNF = rep(c(0,0,1,1,0,0,1,1,0,0,1,1),4)

participantVectorNF = sort(rep(as.numeric(levels(freqTableNF[!duplicated(lastMouseCoordsDataFrame$participant_id),]$participant_id)),48))

#each number here is numebr of ppts
placeholderDFNF = data.frame("participant_id" = participantVectorNF, "species" = rep(speciesVectorNF,30), 
                             "blurred" = rep(blurredVectorNF,30), 
                             "lookRear" = rep(lookRearVectorNF,30),"size" = rep(sizeVectorNF,30), "freq" = NA)
placeholderDFNF = subset(placeholderDFNF, !(placeholderDFNF$participant_id %in% c(2639221,2639155)))

filledDFNF = merge(freqTableNF,placeholderDFNF,by = c("participant_id","species","lookRear","blurred","size"),all.y=T)
filledDFNF = filledDFNF[,1:6]
colnames(filledDFNF) = c("participant_id","species","lookRear","blurred","size","freq")
filledDFNF[is.na(filledDFNF$freq),]$freq = 0

library(dplyr)
freqTableDFNF = data.frame(filledDFNF %>%
                             group_by(participant_id, species, blurred, size) %>%
                             mutate(prop = freq/sum(freq)))
detach("package:useful", unload=TRUE)
detach("package:dplyr", unload=TRUE)

freqTableDFRearNF = freqTableDFNF[freqTableDFNF$lookRear==1,]
freqTableDFRearNF$prop = freqTableDFRearNF$prop - 0.5

freqTableDFRearNF$blurred = as.factor(freqTableDFRearNF$blurred)
freqTableDFRearNF$species = as.factor(freqTableDFRearNF$species)
freqTableDFRearNF$size = as.factor(freqTableDFRearNF$size)
freqTableDFRearNF$participant_id = as.factor(freqTableDFRearNF$participant_id)

freqTableDFRearNF[freqTableDFRearNF$freq==0,]$prop = 0

propSummaryLookRearNF = summarySEwithin(freqTableDFRearNF, measurevar="prop", withinvars=c("blurred","species","size"),
                                        idvar="participant_id", na.rm=FALSE, conf.interval=.95)

#statistical analysis below

ezANOVA(freqTableDFRearNF, prop, participant_id, within=.(blurred,species,size), detailed=T, type = 3)

freqTableDFRearGNF = freqTableDFRearNF[freqTableDFRearNF$species=="G",]
freqTableDFRearMNF = freqTableDFRearNF[freqTableDFRearNF$species=="M",]
freqTableDFRearSNF = freqTableDFRearNF[freqTableDFRearNF$species=="S",]
freqTableDFRearNNF = freqTableDFRearNF[freqTableDFRearNF$species=="N",]

ezANOVA(freqTableDFRearGNF, prop, participant_id, within=.(blurred,size), detailed=T, type = 3)
ezANOVA(freqTableDFRearMNF, prop, participant_id, within=.(blurred,size), detailed=T, type = 3)
ezANOVA(freqTableDFRearSNF, prop, participant_id, within=.(blurred,size), detailed=T, type = 3)
ezANOVA(freqTableDFRearNNF, prop, participant_id, within=.(blurred,size), detailed=T, type = 3)

propSummaryLookRearNF$blurred <- factor(propSummaryLookRearNF$blurred, levels = c("LOW","MED","AC"))

#Grevy
t.test(prop ~ size, data = freqTableDFRearGNF[freqTableDFRearGNF$blurred == "LOW",], paired = TRUE)
cohensD(prop ~ size, data   = freqTableDFRearGNF[freqTableDFRearGNF$blurred == "LOW",], method = "paired")
t.test(prop ~ size, data = freqTableDFRearGNF[freqTableDFRearGNF$blurred == "MED",], paired = TRUE)
cohensD(prop ~ size, data   = freqTableDFRearGNF[freqTableDFRearGNF$blurred == "MED",], method = "paired")
t.test(prop ~ size, data = freqTableDFRearGNF[freqTableDFRearGNF$blurred == "AC",], paired = TRUE)
cohensD(prop ~ size, data   = freqTableDFRearGNF[freqTableDFRearGNF$blurred == "AC",], method = "paired")

#Mountain
t.test(prop ~ size, data = freqTableDFRearMNF[freqTableDFRearMNF$blurred == "LOW",], paired = TRUE)
cohensD(prop ~ size, data   = freqTableDFRearMNF[freqTableDFRearMNF$blurred == "LOW",], method = "paired")
t.test(prop ~ size, data = freqTableDFRearMNF[freqTableDFRearMNF$blurred == "MED",], paired = TRUE)
cohensD(prop ~ size, data   = freqTableDFRearMNF[freqTableDFRearMNF$blurred == "MED",], method = "paired")
t.test(prop ~ size, data = freqTableDFRearMNF[freqTableDFRearMNF$blurred == "AC",], paired = TRUE)
cohensD(prop ~ size, data   = freqTableDFRearMNF[freqTableDFRearMNF$blurred == "AC",], method = "paired")

#PlainsS
t.test(prop ~ size, data = freqTableDFRearSNF[freqTableDFRearSNF$blurred == "LOW",], paired = TRUE)
cohensD(prop ~ size, data   = freqTableDFRearSNF[freqTableDFRearSNF$blurred == "LOW",], method = "paired")
t.test(prop ~ size, data = freqTableDFRearSNF[freqTableDFRearSNF$blurred == "MED",], paired = TRUE)
cohensD(prop ~ size, data   = freqTableDFRearSNF[freqTableDFRearSNF$blurred == "MED",], method = "paired")
t.test(prop ~ size, data = freqTableDFRearSNF[freqTableDFRearSNF$blurred == "AC",], paired = TRUE)
cohensD(prop ~ size, data   = freqTableDFRearSNF[freqTableDFRearSNF$blurred == "AC",], method = "paired")

#PlainsN
t.test(prop ~ size, data = freqTableDFRearNNF[freqTableDFRearNNF$blurred == "LOW",], paired = TRUE)
cohensD(prop ~ size, data   = freqTableDFRearNNF[freqTableDFRearNNF$blurred == "LOW",], method = "paired")
t.test(prop ~ size, data = freqTableDFRearNNF[freqTableDFRearNNF$blurred == "MED",], paired = TRUE)
cohensD(prop ~ size, data   = freqTableDFRearNNF[freqTableDFRearNNF$blurred == "MED",], method = "paired")
t.test(prop ~ size, data = freqTableDFRearNNF[freqTableDFRearNNF$blurred == "AC",], paired = TRUE)
cohensD(prop ~ size, data   = freqTableDFRearNNF[freqTableDFRearNNF$blurred == "AC",], method = "paired")


#Mountain
tempMLow = freqTableDFRearMNF[freqTableDFRearMNF$blurred == "LOW",]
t.test(tempMLow[tempMLow$size==256,]$prop, mu = 0, alternative = "greater")
cohensD(tempMLow[tempMLow$size==256,]$prop, mu = 0)

tempMMed = freqTableDFRearMNF[freqTableDFRearMNF$blurred == "MED",]
t.test(tempMMed[tempMMed$size==165,]$prop, mu = 0, alternative = "greater")
cohensD(tempMMed[tempMMed$size==165,]$prop, mu = 0)
t.test(tempMMed[tempMMed$size==256,]$prop, mu = 0, alternative = "greater")
cohensD(tempMMed[tempMMed$size==256,]$prop, mu = 0)

tempMAc = freqTableDFRearMNF[freqTableDFRearMNF$blurred == "AC",]
t.test(tempMAc[tempMAc$size==165,]$prop, mu = 0, alternative = "greater")
cohensD(tempMAc[tempMAc$size==165,]$prop, mu = 0)
t.test(tempMAc[tempMAc$size==256,]$prop, mu = 0, alternative = "greater")
cohensD(tempMAc[tempMAc$size==256,]$prop, mu = 0)

#PlainsN
tempNLow = freqTableDFRearNNF[freqTableDFRearNNF$blurred == "LOW",]
t.test(tempNLow[tempNLow$size==256,]$prop, mu = 0, alternative = "greater")
cohensD(tempNLow[tempNLow$size==256,]$prop, mu = 0)

tempNMed = freqTableDFRearNNF[freqTableDFRearNNF$blurred == "MED",]
t.test(tempNMed[tempNMed$size==165,]$prop, mu = 0, alternative = "greater")
cohensD(tempNMed[tempNMed$size==165,]$prop, mu = 0)
t.test(tempNMed[tempNMed$size==256,]$prop, mu = 0, alternative = "greater")
cohensD(tempNMed[tempNMed$size==256,]$prop, mu = 0)

tempNAc = freqTableDFRearNNF[freqTableDFRearNNF$blurred == "AC",]
t.test(tempNAc[tempNAc$size==165,]$prop, mu = 0, alternative = "greater")
cohensD(tempNAc[tempNAc$size==165,]$prop, mu = 0)
t.test(tempNAc[tempNAc$size==256,]$prop, mu = 0, alternative = "greater")
cohensD(tempNAc[tempNAc$size==256,]$prop, mu = 0)