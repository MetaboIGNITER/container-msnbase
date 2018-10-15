#!/usr/bin/env Rscript

options(stringAsfactors = FALSE, useFancyQuotes = FALSE)

# Taking the command line arguments
args <- commandArgs(trailingOnly = TRUE)

if(length(args)==0)stop("No files have been specified!")


################# the following packages have been adapted from msnbasea and maldiquant
#################


numberOfCommonPeaks <- function(x, y, tolerance=25e-6, relative=TRUE) {
  sum(commonPeaks(x, y, tolerance=tolerance, relative=relative))
}

commonPeaks <- function(x, y, method=c("highest", "closest"),
                        tolerance=25e-6, relative=TRUE) {
  m <- matchPeaks(x, y, method=match.arg(method), tolerance=tolerance,
                  relative=relative)
  
  m[which(is.na(m))] <- 0L
  
  as.logical(m)
}


matchPeaks <- function(x, y, method=c("highest", "closest", "all"),
                       tolerance=25e-6, relative=TRUE) {
  method <- match.arg(method)
  
  y<-y[,"mz"]
  
  if (nrow(x) == 0 || length(y) == 0) {
    return(integer(peaksCount(x)))
  }
  
  m <- relaxedMatch((x[,"mz"]), y, nomatch=NA, tolerance=tolerance,
                    relative=relative)
  
  if (anyDuplicated(m)) {
    o <- order((x[,"intensity"]), decreasing=TRUE)
    sortedMatches <- m[o]
    sortedMatches[which(duplicated(sortedMatches))] <- NA
    m[o] <- sortedMatches
  }
  
  as.integer(m)
}


relaxedMatch <- function(x, table, nomatch=NA_integer_, tolerance=25e-6,
                         relative=TRUE) {
  
  if (relative) {
    if (tolerance > 1L) {
      stop(sQuote("tolerance"),
           " must be smaller than 1 for relative deviations.")
    }
    tolerance <- table*tolerance
  }
  
  MALDIquant:::match.closest(x, table, tolerance=tolerance, nomatch=nomatch)
}

bin_Spectra <- function(object1, object2, binSize = 1L,
                        breaks = seq(floor(min(c((object1[,"mz"]), (object2[,"mz"])))),
                                     ceiling(max(c((object1[,"mz"]), (object2[,"mz"])))),
                                     by = binSize)) {
  breaks <- .fix_breaks(breaks, range((object1[,"mz"]), (object2[,"mz"])))
  list(bin_Spectrum(object1, breaks = breaks),
       bin_Spectrum(object2, breaks = breaks))
}

bin_Spectrum <- function(object, binSize = 1L,
                         breaks = seq(floor(min((object[,"mz"]))),
                                      ceiling(max((object[,"mz"]))),
                                      by = binSize),
                         fun = sum,
                         msLevel.) {
  ## If msLevel. not missing, perform the trimming only if the msLevel
  ## of the spectrum matches (any of) the specified msLevels.
  # print(length(object[,"intensity"]))
  # print(length(object[,"mz"]))
  bins <- .bin_values(object[,"intensity"], object[,"mz"], binSize = binSize,
                      breaks = breaks, fun = fun)
  
  
  #data.frame(bins$mids,bins$x)
  object<-matrix(c(bins$mids,bins$x), nrow = length(bins$x), dimnames = list(1:length(bins$x), c("mz","intensity")))
  # object<-data.frame("mz" = bins$mids,"intensity" = bins$x)
  return(object)
}

dotproduct<-function (x, y) 
{
  as.vector(x %*% y)/(sqrt(sum(x * x)) * sqrt(sum(y * y)))
}
compare_Spectra2 <- function(x, y,
                             fun=c("common", "cor", "dotproduct"),binSize=0,tolerance =0,relative = F) {
  {
    
    if (fun == "cor" || fun == "dotproduct") {
      #binnedSpectra <- bin_Spectra(x, y, ...)
      
      binSize<-binSize
      
      breaks = seq(floor(min(c((x[,"mz"]), (y[,"mz"])))),
                   ceiling(max(c((x[,"mz"]), (y[,"mz"])))),
                   by = binSize)
      
      # breaks <- .fix_breaks(brks = breaks, rng = range((x[,"mz"]), (y[,"mz"]))) # fix!
      
      brks = breaks
      rng = range((x[,"mz"]), (y[,"mz"]))
      if (brks[length(brks)] <= rng[2])
        breaks <- c(brks, max((rng[2] + 1e-6),
                              brks[length(brks)] + mean(diff(brks))))
      
      
      breaks1 = breaks
      rng = range(x[,"mz"])
      if (brks[length(brks)] <= rng[2])
        breaks1 <- c(brks, max((rng[2] + 1e-6),
                               brks[length(brks)] + mean(diff(brks))))
      
      nbrks <- length(breaks1)
      idx <- findInterval(x[,"mz"], breaks1)
      ## Ensure that indices are within breaks.
      idx[which(idx < 1L)] <- 1L
      idx[which(idx >= nbrks)] <- nbrks - 1L
      
      ints <- double(nbrks - 1L)
      ints[unique(idx)] <- unlist(lapply(base::split(x[,"intensity"], idx), sum),
                                  use.names = FALSE)
      binsx=ints
      binsmids = (breaks1[-nbrks] + breaks1[-1L]) / 2L  
      
      list1<-binsx#matrix(c(binsmids,binsx), nrow = length(binsx), dimnames = list(1:length(binsx), c("mz","intensity")))
      
      
      # breaks2 <- .fix_breaks(breaks, range(y[,"mz"]))
      breaks2 = breaks
      brks = breaks
      rng = range(y[,"mz"])
      if (brks[length(brks)] <= rng[2])
        breaks2 <- c(brks, max((rng[2] + 1e-6),
                               brks[length(brks)] + mean(diff(brks))))
      nbrks <- length(breaks2)
      idx <- findInterval(y[,"mz"], breaks2)
      ## Ensure that indices are within breaks.
      idx[which(idx < 1L)] <- 1L
      idx[which(idx >= nbrks)] <- nbrks - 1L
      
      ints <- double(nbrks - 1L)
      ints[unique(idx)] <- unlist(lapply(base::split(y[,"intensity"], idx), sum),
                                  use.names = FALSE)
      bins<-list(x = ints, mids = (breaks2[-nbrks] + breaks[-1L]) / 2L)
      binsmids = (breaks2[-nbrks] + breaks2[-1L]) / 2L  
      binsx=ints
      
      list2<-binsx#matrix(c(binsmids,binsx), nrow = length(binsx), dimnames = list(1:length(binsx), c("mz","intensity")))
      
      
      #inten <- lapply(list(list1,list2), function(x){x[,"intensity"]})
      # ifelse(fun == "dotproduct",dotproduct(list1,list2),cor(list1,list2)
      return(c(dotproduct(list1,list2),cor(list1,list2)))
    } else if (fun == "common") {
      return(numberOfCommonPeaks(x, y, tolerance =tolerance,relative = relative))
    }
  } 
}

.fix_breaks <- function(brks, rng) {
  ## Assuming breaks being sorted.
  if (brks[length(brks)] <= rng[2])
    brks <- c(brks, max((rng[2] + 1e-6),
                        brks[length(brks)] + mean(diff(brks))))
  brks
}
.bin_values <- function(x, toBin, binSize = 1, breaks = seq(floor(min(toBin)),
                                                            ceiling(max(toBin)),
                                                            by = binSize),
                        fun = max) {
  if (length(x) != length(toBin))
    stop("lengths of 'x' and 'toBin' have to match.")
  fun <- match.fun(fun)
  breaks <- .fix_breaks(breaks, range(toBin))
  nbrks <- length(breaks)
  idx <- findInterval(toBin, breaks)
  ## Ensure that indices are within breaks.
  idx[which(idx < 1L)] <- 1L
  idx[which(idx >= nbrks)] <- nbrks - 1L
  
  ints <- double(nbrks - 1L)
  ints[unique(idx)] <- unlist(lapply(base::split(x, idx), fun),
                              use.names = FALSE)
  list(x = ints, mids = (breaks[-nbrks] + breaks[-1L]) / 2L)
  
}



#################
#################
#################



inputLibrary<-NA
inputMS2<-NA
outputCSV<-NA
readNameOftheMS2<-NA
# MS1 PPM tol
precursorPPMTol<-10

# MS2 abs tol
fragmentabsTol<-0.07

# MS2 ppm tol
fragmentPPMTol<-10

# MS1 RT tol
precursorRTTol<-20

# search based on feature RT or parent MS2
searchRange<-T

# preprocess MS2 ?
preprocess<-F

# estimate decoy ?
outputSemiDecoy<-T

# how many we should peak 0 is top -1 is all 
topHits<--1

# set ionization mode pos or neg
ionMode<-"pos"
# which score to choose
topScore<-"Scoredotproduct"

# resample
resample=1000

# read the parameters
for(arg in args)
{
  argCase<-strsplit(x = arg,split = "=")[[1]][1]
  value<-strsplit(x = arg,split = "=")[[1]][2]
  
  if(argCase=="inputMS2")
  {
    inputMS2=as.character(value)
  }
    if(argCase=="inputLibrary")
  {
    inputLibrary=as.character(value)
  }
    if(argCase=="readNameOftheMS2")
  {
    readNameOftheMS2=as.character(value)
  }
    if(argCase=="outputCSV")
  {
    outputCSV=as.character(value)
  }
    if(argCase=="precursorPPMTol")
  {
    precursorPPMTol=as.numeric(value)
  }
    if(argCase=="fragmentabsTol")
  {
    fragmentabsTol=as.numeric(value)
  }
      if(argCase=="fragmentPPMTol")
  {
    fragmentPPMTol=as.numeric(value)
  }
        if(argCase=="precursorRTTol")
  {
    precursorRTTol=as.numeric(value)
  }
        if(argCase=="searchRange")
  {
    searchRange=as.logical(value)
  }
  
          if(argCase=="outputSemiDecoy")
  {
    outputSemiDecoy=as.logical(value)
  }
  
    
          if(argCase=="topHits")
  {
    topHits=as.numeric(value)
  }
           if(argCase=="ionMode")
  {
    ionMode=as.character(value)
  }
  
           if(argCase=="topScore")
  {
    topScore=as.character(value)
  }
  
           if(argCase=="resample")
  {
    resample=as.numeric(value)
  }
 

  
}
# load MSnbase package for comparing spectra




require(MSnbase)
require(intervals)
if(is.na(inputLibrary) | is.na(inputMS2)) stop("Both inputs (library & MS2) are required")
# read library file
######################### UNCOMMENT #############################
MSlibrary<-read.csv(inputLibrary,stringsAsFactors = F)
if(nrow(MSlibrary)<1) stop("Library is empty!")


checkAllCol<-all(sapply(c("startRT","endRT","startMZ","endMZ","centermz","centerrt","intensity","MS2fileName" ,"ID","Name", "nmass","MS2mz","MS2rt","MS2intensity","MS2mzs","MS2intensities","featureGroup" ),
                        function(x){x%in%colnames(MSlibrary)}))

if(!checkAllCol)stop("Check the library file! it has to include the following columns:
                               startRT,endRT,startMZ,endMZ,centermz,centerrt, intensity,fileName ,ID,Name, nmass,MS2mz,MS2rt,MS2intensity,MS2mzs,featureGroup and MS2intensities" )

# limit the library to those with MS2
MSlibrary<-MSlibrary[MSlibrary[,"MS2mzs"]!="",]
# extract name of the file
MS2NameFileName<-readNameOftheMS2
if(is.na(MS2NameFileName))MS2NameFileName<-basename(inputMS2)

# extract precursor RT and mz from name of the file
precursorRT<-as.numeric(strsplit(x = MS2NameFileName,split = "_",fixed = T)[[1]][2])
precursorMZ<-as.numeric(strsplit(x = MS2NameFileName,split = "_",fixed = T)[[1]][3])
if(is.na(precursorRT) | is.na(precursorMZ)) stop("File name does not contain RT or mz. Check the file name!")

# extract MS2 information from MS2 file
MS2Information<-readLines(inputMS2)

# split the data into a dataframe for easy access
MS2DataFrame<-sapply(X = strsplit(x = MS2Information,split = " ",fixed = T)[[1]],FUN = function(x){strsplit(x = x,split = "=",fixed = T)[[1]]})
colnames(MS2DataFrame)<-MS2DataFrame[1,]

# set parameters from input file
# mz ppm
if(!"DatabaseSearchRelativeMassDeviation"%in%colnames(MS2DataFrame))stop("DatabaseSearchRelativeMassDeviation is not in the parameter file!")
if(is.na(precursorPPMTol))precursorPPMTol<-as.numeric(MS2DataFrame[2,"DatabaseSearchRelativeMassDeviation"])
if(is.na(precursorPPMTol)) stop("precursorPPMTol has not been provided and is not in the parameter file!")

# fragment absolute deviation
if(!"FragmentPeakMatchAbsoluteMassDeviation"%in%colnames(MS2DataFrame))stop("FragmentPeakMatchAbsoluteMassDeviation is not in the parameter file!")
if(is.na(fragmentabsTol))fragmentabsTol<-as.numeric(MS2DataFrame[2,"FragmentPeakMatchAbsoluteMassDeviation"])
if(is.na(fragmentabsTol)) stop("fragmentabsTol has not been provided and is not in the parameter file!")

# fragment PPM deviation

if(is.na(fragmentPPMTol))fragmentPPMTol<-as.numeric(MS2DataFrame[2,"FragmentPeakMatchRelativeMassDeviation"])
if(is.na(fragmentPPMTol)) stop("fragmentPPMTol has not been provided and is not in the parameter file!")
# fix if RT tol is missing!
if(is.na(precursorRTTol)) {cat("WARNNING: precursorRTTol has not been provided using largest RT region:",.Machine$double.xmax);precursorRTTol<-.Machine$double.xmax}

parentFile<-"NotFound"
if(!"SampleName"%in%colnames(MS2DataFrame))stop("SampleName is not in the parameter file!")
parentFile<-basename(MS2DataFrame[2,"SampleName"])
if(parentFile=="NotFound") {warning("SampleName was not found in the parameter file setting to NotFound")}

# extract MS2 peaks
if(!"PeakListString"%in%colnames(MS2DataFrame))stop("PeakListString is not in the parameter file!")
MS2TMP<-MS2DataFrame[2,"PeakListString"]
# if spectrum is empty do not continue
isMS2Emtpy<-MS2TMP==""

# define a function for calculating ppm!

ppmCal<-function(run,ppm)
{
  return((run*ppm)/1000000)
}
temp<-NA
if(!isMS2Emtpy)
{
  
  temp<-t(sapply(X=(strsplit(x = strsplit(x = MS2TMP,split = ";",fixed = T)[[1]],split = "_",fixed = T)),FUN = function(x){c(mz=as.numeric(x[1]),
                                                                                                                             int=as.numeric(x[2]))}))
  temp<-temp[temp[,1]!=0,]
  if(nrow(temp)<1)isMS2Emtpy<-F
}

if(!isMS2Emtpy){
  # read MS2 and convert to dataframe
  
  # create a msnbase file!
  targetMS2<-new("Spectrum2", mz=temp[,1], intensity=temp[,2])
  targetMS2DataFrame<-data.frame(mz=temp[,1],intensity=temp[,2])
  targetMS2DataFrame<- matrix(c(temp[,1],temp[,2]), nrow = length(temp[,2]), dimnames = list(1:length(temp[,2]), c("mz","intensity")))
  #names(targetMS2DataFrame)<-c("mz","intensity")
  # set search interval for MS2
  mzTarget<-Intervals_full(cbind(precursorMZ,precursorMZ))
  rtTarget<-Intervals_full(cbind(precursorRT,precursorRT))
  
  # set search interval for library that is either as range or centroid
  mzLib<-NA
  rtLib<-NA
  if(searchRange)
  {
    if(!"startMZ"%in%colnames(MSlibrary))stop("startMZ is not in library file!")
    if(!"endMZ"%in%colnames(MSlibrary))stop("endMZ is not in library file!")
    mzLib<-Intervals_full(cbind(MSlibrary$startMZ-
                                  ppmCal(MSlibrary$startMZ,precursorPPMTol),
                                MSlibrary$endMZ+
                                  ppmCal(MSlibrary$endMZ,precursorPPMTol)))
    
    if(!"startRT"%in%colnames(MSlibrary))stop("startRT is not in library file!")
    if(!"endRT"%in%colnames(MSlibrary))stop("endRT is not in library file!")
    rtLib<-Intervals_full(cbind(MSlibrary$startRT-
                                  precursorRTTol,
                                MSlibrary$endRT+
                                  precursorRTTol))
  }else{
    
    if(!"centermz"%in%colnames(MSlibrary))stop("centermz is not in library file!")
    mzLib<-Intervals_full(cbind(MSlibrary$centermz-
                                  ppmCal(MSlibrary$centermz,precursorPPMTol),
                                MSlibrary$centermz+
                                  ppmCal(MSlibrary$centermz,precursorPPMTol)))
    if(!"centerrt"%in%colnames(MSlibrary))stop("centerrt is not in library file!")
    rtLib<-Intervals_full(cbind(MSlibrary$centerrt-
                                  precursorRTTol,
                                MSlibrary$centerrt+
                                  precursorRTTol))
  }
  
  # do precursor mass search
  Mass_iii <- interval_overlap(mzTarget,mzLib)
  
  # do precursor mass search
  Time_ii <- interval_overlap(rtTarget,rtLib)
  
  # check if there is any hit ?!
  imatch = mapply(intersect,Time_ii,Mass_iii)
  foundHit<-length(imatch[[1]])>0
  
  # create an empty array for the results
  results<-c()
  
 # compareSpectra needs relative deviation in fractions
fragmentPPMTol<-fragmentPPMTol/1000000
  
  if(foundHit)
  {
    cat("Number of hits: ",length(imatch),"\n")
    for(i in imatch)
    {
      hitTMP<-MSlibrary[i,]

      if(hitTMP[,"MS2mzs"]!="")
     { 

     
       # tmpResults<-c()
       # for(j in 1:length(parentmzs))
        {
          #MS2sTMPLib<-parentMS2s[[j]]
          tempmz<-as.numeric(strsplit(hitTMP[,"MS2mzs"],split = ";",fixed = T)[[1]])
          tempint<-as.numeric(strsplit(hitTMP[,"MS2intensities"],split = ";",fixed = T)[[1]])
          
          tempmz<-tempmz[tempint!=0]
          tempint<-tempint[tempint!=0]
    
          
          # if all the peaks were zero, skip rest of the loop!
          if(length(tempint)<1)next
          
          # extract name of the metabolite

          hitName<-hitTMP[,"Name"]
          hitChI<-NA
          Identifier<-ifelse(test = is.null(results),yes = 1,no = (nrow(results)+1))
          fileName<-parentFile
          parentMZ<-precursorMZ
          parentRT<-precursorRT
        
          restOfLibInformation<-hitTMP[,c("startRT","endRT","startMZ","endMZ","centermz","centerrt","intensity","MS2fileName" ,"ID","Name", "nmass","MS2mz","MS2rt","MS2intensity","featureGroup" )]
          featureGroup<-hitTMP[,"featureGroup"]
          nmass<-as.numeric(hitTMP[,"nmass"])
          MS1mzTolTh<-NA
         # if(ionMode=="pos")
          MS1mzTolTh<-((parentMZ-(nmass))/(nmass))*1000000
         # if(ionMode=="neg")
         # MS1mzTolTh<-((parentMZ-(nmass-1.00727647))/(nmass-1.00727647))*1000000
          
          MS1RTTol<-NA
          centerRT<-as.numeric(hitTMP[,"centerrt"])
          MS1RTTol<-parentRT-centerRT
          # create MS2 object
          libMS2Obj<-new("Spectrum2", mz=tempmz, intensity=tempint)
          

          
          ## output three types of scores: dotproduct, common peaks and correlation (dotproduct will be our main score)
          dotPeaks<-NA
          tryCatch({
            dotPeaks<-compareSpectra(targetMS2, libMS2Obj, fun="dotproduct",binSize =fragmentabsTol)
          }, warning = function(w) {
          }, error = function(e) {
           # cat(e,"\n")
          }, finally = {
           # dotPeaks<-"Error"
          })
          
          nPeakCommon<-NA
          tryCatch({
            nPeakCommon<-compareSpectra(targetMS2, libMS2Obj, fun="common",tolerance =fragmentPPMTol,relative = TRUE)
          }, warning = function(w) {
          }, error = function(e) {
           # cat(e,"\n")
          }, finally = {
           # nPeakCommon<-"Error"
          })
          
          corPeaks<-NA
          tryCatch({
            corPeaks<-compareSpectra(targetMS2, libMS2Obj, fun="cor",binSize =fragmentabsTol)
          }, warning = function(w) {
          }, error = function(e) {
          #  cat(e,"\n")
          }, finally = {
           # corPeaks<-"Error"
          })
          ##
          # if requrested create a decoyDatabase for this specific MS2 and estimate "e-value"
          # this will repeat the whole process but for rest of the peaks
          # it will take LONG time!
          decoyScore<-list()
          decoyScore[["dotproduct"]]<-c()
          decoyScore[["common"]]<-c()
          decoyScore[["cor"]]<-c()
          if(outputSemiDecoy)
          {
            decoyLib<-MSlibrary[-as.vector(imatch),]
            decoyScoresTMP<-c()
            #nrow(decoyLib)
            start_time <- Sys.time()
            allMzs<-as.character(decoyLib[,"MS2mzs"])
            allInts<-as.character(decoyLib[,"MS2intensities"])
            
            rowNumbersDecoy<-1:nrow(decoyLib)
            if(resample>0)
            {
              set.seed(resample)
              rowNumbersDecoy<-sample(x = rowNumbersDecoy,size = resample,replace = F)
            }
            for(k in rowNumbersDecoy)
            {
             tempmzDecoy<-as.numeric(strsplit(allMzs[k],split = ";",fixed = T)[[1]])
             tempintDecoy<-as.numeric(strsplit(allInts[k],split = ";",fixed = T)[[1]])
             tempmzDecoy<-tempmzDecoy[tempintDecoy!=0]
             tempintDecoy<-tempintDecoy[tempintDecoy!=0]
  
             # if(length(tempintDecoy)<1)next
             tempDecoy<- matrix(c(tempmzDecoy,tempintDecoy), nrow = length(tempintDecoy), dimnames = list(1:length(tempintDecoy), c("mz","intensity")))
            #  tempDecoy<-list(mz=tempmzDecoy,intensity=tempintDecoy)
         #   tempDecoy<-data.frame(mz=tempmzDecoy,intensity=tempintDecoy)
              dotPeaksDecoy<-NA
              nPeakCommonDecoy<-NA
              tryCatch({
                tmp<-c(NA,NA)
                tmp<-compare_Spectra2(targetMS2DataFrame, tempDecoy, fun="dotproduct",binSize =fragmentabsTol)
                dotPeaksDecoy<-tmp[1]
                corPeaksDecoy<-tmp[2]
              }, warning = function(w) {
              }, error = function(e) {
                #cat(e,"\n")
              }, finally = {
                #  dotPeaksDecoy<-"Error"
              })
              
              nPeakCommonDecoy<-NA
              tryCatch({
                nPeakCommonDecoy<-compare_Spectra2(targetMS2DataFrame, tempDecoy, fun="common",tolerance =fragmentPPMTol,relative = TRUE)
              }, warning = function(w) {
              }, error = function(e) {
              #  cat(e,"\n")
              }, finally = {
                # nPeakCommonDecoy<-"Error"
              })
              # 
              # corPeaksDecoy<-NA
              # tryCatch({
              #   corPeaksDecoy<-compare_Spectra2(targetMS2, tempDecoy, fun="cor",binSize =fragmentabsTol)
              # }, warning = function(w) {
              # }, error = function(e) {
              #  # cat(e,"\n")
              # }, finally = {
              #   #corPeaks<-"Error"
              # })
              
             decoyScore[["dotproduct"]]<-c(decoyScore[["dotproduct"]],dotPeaksDecoy)
             decoyScore[["common"]]<-c(decoyScore[["common"]],nPeakCommonDecoy)
            decoyScore[["cor"]]<-c(decoyScore[["cor"]],corPeaksDecoy)
              
            }
           
            # for(k in 1:nrow(decoyLib))
            # {
            #   hitTMPDecoy<-MSlibrary[k,]
            #   
            #   parentmzsDecoy<-strsplit(x = hitTMPDecoy[,"parentmzs"],split = ";",fixed = T)[[1]]
            #   parentMS2sDecoy<-sapply(strsplit(x = hitTMPDecoy[,"MS2s"],split = ";",fixed = T)[[1]],func)
            #   for(p in 1:length(parentmzsDecoy))
            #   {
            #      MS2sTMPLibDecoy<-parentMS2sDecoy[[p]]
            #      asd<-(strsplit(x = strsplit(x = MS2sTMPLibDecoy,split = ":",fixed = T)[[1]],split = "_",fixed = T))
            #     # df <- data.frame(matrix(unlist(asd,use.names=FALSE), nrow=length(asd), byrow=T))
            #   #    n <- length(asd[[1]])
            #   #    tempDecoy <- t(structure(asd, row.names = c(NA, -n), class = "data.frame"))
            #   #    tempDecoy<-data.frame(tempDecoy)
            #   #    names(tempDecoy)<-c("mz","intensity")
            #   #    tempDecoy$mz<-as.numeric(tempDecoy$mz)
            #   #    tempDecoy$intensity<-as.numeric(tempDecoy$intensity)
            #   #   # tempDecoy<-t(sapply(X=(strsplit(x = strsplit(x = MS2sTMPLibDecoy,split = ":",fixed = T)[[1]],split = "_",fixed = T)),FUN = function(x){c(mz=as.numeric(x[1]),
            #   #   #                                                                                                                                         int=as.numeric(x[2]))}))
            #   #    tempDecoy<-data.frame(tempDecoy[tempDecoy[,1]!=0,])
            #   # #   names(tempDecoy)<-c("mz","intensity")
            #   # #   # if all the peaks were zero, skip rest of the loop!
            #   #   if(nrow(tempDecoy)<1)next
            #   # aaa<- compare_Spectra(tempDecoy,targetMS2DataFrame,fun="common",tolerance =0.001,relative = TRUE)
            #     #libMS2ObjDecoy<-new("Spectrum2", mz=tempDecoy[,1], intensity=tempDecoy[,2])
            # 
            # 
            #     ## output three types of scores: dotproduct, common peaks and correlation (dotproduct will be our main score)
            #     dotPeaksDecoy<-NA
            #     tryCatch({
            #       dotPeaksDecoy<-compare_Spectra(targetMS2DataFrame, tempDecoy, fun="dotproduct",binSize =fragmentabsTol)
            #     }, warning = function(w) {
            #     }, error = function(e) {
            #       cat(e,"\n")
            #     }, finally = {
            #     #  dotPeaksDecoy<-"Error"
            #     })
            # 
            #     nPeakCommonDecoy<-NA
            #     tryCatch({
            #       nPeakCommonDecoy<-compareSpectra(targetMS2, libMS2ObjDecoy, fun="common",tolerance =fragmentPPMTol,relative = TRUE)
            #     }, warning = function(w) {
            #     }, error = function(e) {
            #       cat(e,"\n")
            #     }, finally = {
            #      # nPeakCommonDecoy<-"Error"
            #     })
            # 
            #     corPeaksDecoy<-NA
            #     tryCatch({
            #       corPeaksDecoy<-compareSpectra(targetMS2, libMS2ObjDecoy, fun="cor",binSize =fragmentabsTol)
            #     }, warning = function(w) {
            #     }, error = function(e) {
            #       cat(e,"\n")
            #     }, finally = {
            #       #corPeaks<-"Error"
            #     })
            # 
            #    decoyScore[["dotproduct"]]<-c(decoyScore[["dotproduct"]],dotPeaksDecoy)
            #    decoyScore[["common"]]<-c(decoyScore[["common"]],nPeakCommonDecoy)
            #    decoyScore[["cor"]]<-c(decoyScore[["cor"]],corPeaksDecoy)
            #   }
            #   
            # }
            # end_time <- Sys.time()
            # print(end_time - start_time)
            decoyScore[["dotproduct"]]<-na.omit(decoyScore[["dotproduct"]])
            decoyScore[["common"]]<-na.omit(decoyScore[["common"]])
            decoyScore[["cor"]]<-na.omit(decoyScore[["cor"]])
            dotPeaksDecoy<-NA
            nPeakCommonDecoy<-NA
            corPeaksDecoy<-NA
            if(!is.na(dotPeaks))
            dotPeaksDecoy<-sum(decoyScore[["dotproduct"]]>dotPeaks)/length(decoyScore[["dotproduct"]]   )
            if(!is.na(nPeakCommon))
            nPeakCommonDecoy<-sum(decoyScore[["common"]]>nPeakCommon)/length(decoyScore[["common"]])
            if(!is.na(corPeaks))
            corPeaksDecoy<-sum(decoyScore[["cor"]]>corPeaks)/length(decoyScore[["cor"]])
            
          }
          results<-rbind(results,  data.frame(fileName=parentFile,parentMZ=parentMZ,parentRT=parentRT,Name=fileName,Identifier=Identifier,InChI=NA,
                                                    MS1mzTolTh,MS1RTTol,
                     Scoredotproduct=dotPeaks,Scorecommon=nPeakCommon,ScoreCorrelation=corPeaks,
                     ScoredotproductEValue=dotPeaksDecoy,ScorecommonEValue=nPeakCommonDecoy,ScoreCorrelationEValue=corPeaksDecoy,
                     score=dotPeaks,scoreEValue=dotPeaksDecoy,restOfLibInformation,MS1RTTol=MS1RTTol,featureGroup=featureGroup))
          
 
        }
 
        
      }
    }
    # limit the results as requrested by user: tophit:-1 = all, tophit:0 = top, tophit:>0 = top tophits score higher the better (for now)
    resTMP<-c()
    if(topHits!=-1 & nrow(results)>1)
    {
      if(topHits==0)
      {
        for(groupNumber in unique(results[,"featureGroup"]))
        {
          tmpResults<-data.frame(results[results[,"featureGroup"]==groupNumber,])
          resTMP<-rbind(resTMP,tmpResults[which.max(tmpResults[,topScore]),])
        }
        
      }else{
        for(groupNumber in unique(results[,"featureGroup"]))
        {
          tmpResults<-data.frame(results[results[,"featureGroup"]==groupNumber,])
          resTMP<-rbind(resTMP,tmpResults[order(tmpResults[,topScore],decreasing = T,na.last = T),][seq(1,topHits),])
        }
        
      }
      results<-resTMP
    }
    
    write.csv(x =results, outputCSV)
  }else{file.create(outputCSV)}
}else{
  file.create(outputCSV)
  }


