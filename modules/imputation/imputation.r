library(imputeLCMD)
library(missForest)
library(dplyr)
library(optparse)
library(parallel)
library(caret)




option_list <- list(
  make_option(c("--input_data"), type="character", help="Path to input CSV"),
  make_option(c("--file_name"), type="character", help="Base file name"),
  make_option(c("--output"), type="character", help="Output CSV file name")
 
)

opt_parser <- OptionParser(option_list=option_list)
args <- parse_args(opt_parser)







removeDataMM_Zero_safe = function(data,percentMiss,percentBelowThresh_I,LowAbundThresh_II,percentMisslowAbund_III, label = TRUE){

  #"""
  #altered Missed Missingness algorithm to exclude 0 values for impsing missing values in the complete data set
  #"""
  data <- as.matrix(data)
  zero_mask <- data == 0 
  zero_safe_data <- data
  zero_safe_data[zero_mask] <- NA # convert all real Zeros to NA values

  percentMisslowAbund_III <- percentMisslowAbund_III/0.9
  percentMVmedAbund <- 0.5*percentMisslowAbund_III
  TotalNumValRemove <- round(percentMiss/100*ncol(zero_safe_data)*nrow(zero_safe_data))

  #sort metablites based on mean (zeros get ignored due to beeing NA)
  temp = sort(apply(data, 1, mean), decreasing = FALSE, index = TRUE)
  avgMetabolite <- temp[["x"]] 
  sortAvgIdx <- temp[["ix"]] 
  
  
  numMetBelowThresh <- round(percentBelowThresh_I/100*length(avgMetabolite))
  belowThreshMet <- zero_safe_data[sortAvgIdx[seq_len(numMetBelowThresh)],]
  belowThreshIdx <- sortAvgIdx[seq_len(numMetBelowThresh)]
  aboveThreshMet <- zero_safe_data[sortAvgIdx[(numMetBelowThresh+1):length(avgMetabolite)],] #High AAG-mtb
  aboveThreshIdx <- sortAvgIdx[(numMetBelowThresh+1):length(avgMetabolite)]

  numMetLowAbund <- round(LowAbundThresh_II/100*nrow(belowThreshMet))
  lowAbundMet <- as.matrix(zero_safe_data[sortAvgIdx[seq_len(numMetLowAbund)],]) #low AAG-mtb
  medAbundmet <- as.matrix(zero_safe_data[sortAvgIdx[(numMetLowAbund+1):nrow(belowThreshMet)],])#Medium AAG-mtb

  
  #numbers of samples per metabolite that will be missing inside ther low AAG-mtb
  numSampBelowAbundThreshLowAbund <- round(percentMisslowAbund_III/100*ncol(lowAbundMet))



  if (numSampBelowAbundThreshLowAbund == 0){dw <- 0}else{dw <- 1} #marker if no missing in low AAG-mtb
  
  #sortLowAbundmet <- t(apply(lowAbundMet, 1, sort)) # on default NA values get put at the end 
  sortLowAbundmet <- t(apply(lowAbundMet, 1, function(x) {sort(x, na.last = TRUE)}))
  #print(sortLowAbundmet)

  sortAbundLowAbundidx <- t(apply(lowAbundMet, 1, order))


  #Matrix of high values of low mean metabolites 
  LowLowAbund = as.matrix(sortLowAbundmet[,seq_len(numSampBelowAbundThreshLowAbund)])
  #print(dim(LowLowAbund))

  #Matrix of low values of low mean metabolites
  HighLowAbund = as.matrix(sortLowAbundmet[,(numSampBelowAbundThreshLowAbund+1):ncol(sortLowAbundmet)]) #include NA / Zero values 
  #print(dim(HighLowAbund))

  alwaysMV = round(0.8*ncol(LowLowAbund)) #allways introduce 80% MNAR values in low values of low mean metabolites
  LowLowAbund[,alwaysMV] = "MNAR" 
  
  #introduce randome MNAR 
  if (alwaysMV < ncol(LowLowAbund)){
    LowLowAbundRandomMV = as.matrix(LowLowAbund[,(alwaysMV+1):ncol(LowLowAbund)])
    LowLowAbund = as.matrix(LowLowAbund[,-((alwaysMV+1):ncol(LowLowAbund))])
    numVal = ncol(LowLowAbundRandomMV)*nrow(LowLowAbundRandomMV)
    LowLowAbundRandomMV_vector = as.vector(t(LowLowAbundRandomMV))
    #zero safe indx
    zero_safe_idx = which(!is.na(LowLowAbundRandomMV_vector))
    #if(length(zero_safe_idx) < numVal){zero_safe_idx <- seq_along(LowLowAbundRandomMV_vector)} #only allows zero to e missing if not enough values are aviable other wise #maybe let it fail?
    #print(numVal)
    numValRemove = round(0.5*numVal)
    #removeIdx = sample(numVal,numValRemove)
    removeIdx = sample(zero_safe_idx,numValRemove)
    temp = as.vector(t(LowLowAbundRandomMV))
    temp[removeIdx] = "MNAR"
    LowLowAbundRandomMV = matrix(temp,nrow = nrow(LowLowAbundRandomMV),ncol = ncol(LowLowAbundRandomMV), byrow = TRUE)
    LowLowAbund = cbind(LowLowAbund, LowLowAbundRandomMV)
    
  }

  lowMNARidx = which(LowLowAbund == "MNAR", arr.ind = TRUE) #safe Index of MNAR values 

  #numbers of samples per metabolite that will be missing inside ther medium AAG-mtb
  numSampBelowAbundThreshMedAbund = round(percentMVmedAbund/100*ncol(medAbundmet))

  if (numSampBelowAbundThreshMedAbund == 0){dz = 0}else{dz=1} #marker if no missing in medium AAG-mtb

  #sortmedAbundmet = t(apply(medAbundmet,1, sort)) # on default NA values get put at the end
  sortmedAbundmet <- t(apply(medAbundmet, 1, function(x) {sort(x, na.last = TRUE)}))

  sortAbundMedAbundidx = t(apply(medAbundmet,1, order))

  LowMedAbund = as.matrix(sortmedAbundmet[,seq_len(numSampBelowAbundThreshMedAbund)])
  
  HighMedAbund = as.matrix(sortmedAbundmet[,(numSampBelowAbundThreshMedAbund+1):ncol(sortmedAbundmet)]) #include NA values 
  alwaysMV = round(0.8*ncol(LowMedAbund))
  LowMedAbund[,seq_len(alwaysMV)] = "MNAR"

  if (alwaysMV < ncol(LowMedAbund)){LowMedAbundRandomMV = as.matrix(LowMedAbund[,(alwaysMV+1):ncol(LowMedAbund)])
    LowMedAbund = as.matrix(LowMedAbund[,-((alwaysMV+1):ncol(LowMedAbund))])
    numVal = ncol(LowMedAbundRandomMV)*nrow(LowMedAbundRandomMV)
    LowMedAbundRandomMV_vector = as.vector(t(LowMedAbundRandomMV))
    #zero safe indx
    zero_safe_idx = which(!is.na(LowMedAbundRandomMV_vector))
    #if(length(zero_safe_idx) < numVal){zero_safe_idx <- seq_along(LowLowAbundRandomMV_vector)} #only allows zero to e missing if not enough values are aviable other wise #maybe let it fail?
    numValRemove = round(0.5*numVal)
    #removeIdx = sample(numVal,numValRemove)
    removeIdx = sample(zero_safe_idx,numValRemove)
    temp = as.vector(t(LowMedAbundRandomMV))
    temp[removeIdx] = "MNAR"
    LowMedAbundRandomMV = matrix(temp,nrow = nrow(LowMedAbundRandomMV),ncol = ncol(LowMedAbundRandomMV),byrow = TRUE)
    LowMedAbund = cbind(LowMedAbund, LowMedAbundRandomMV)
  }
  LMedMNARidx = which(LowMedAbund== "MNAR", arr.ind = TRUE)#safe Index of MNAR values 


  numMVbelowThreshLOD <- (sum(LowMedAbund == "MNAR", na.rm = TRUE) + sum(LowLowAbund == "MNAR", na.rm = TRUE)) * dz
  

  numMVrandom = TotalNumValRemove - numMVbelowThreshLOD #rest of missing values filled up with MCAR
  

  
  if (numMVrandom < 0){numMVrandom = 0}

  # with = without 
  MVrandomAvailable = ncol(zero_safe_data)*nrow(zero_safe_data)-dz*ncol(LowMedAbund)*nrow(LowMedAbund)-dw*ncol(LowLowAbund)*nrow(LowLowAbund)
  aboveThreshWithMV = aboveThreshMet
  HighLowAbundWithMV = HighLowAbund
  HighMedAbundWithMV = HighMedAbund
  size_aboveThreshWithMV = dim(aboveThreshWithMV)
  size_HighLowAbundWithMV = dim(HighLowAbundWithMV)
  size_HighMedAbundWithMV = dim(HighMedAbundWithMV)
  

  remainingValuesVector = c(as.vector(t(aboveThreshWithMV)),as.vector(t(HighLowAbundWithMV)),as.vector(t(HighMedAbundWithMV)))

  zero_safe_Idx <- which(!is.na(remainingValuesVector))   # only sample from zero safe positions
  removeIdx <- sample(zero_safe_Idx, numMVrandom)         
  remainingValuesVector[removeIdx] <- "MCAR"


  #rebuilding matrix

  aboveThreshWithMV = remainingValuesVector[seq_len((ncol(aboveThreshWithMV)*nrow(aboveThreshWithMV)))]
  HighLowAbundWithMV = remainingValuesVector[(length(aboveThreshWithMV)+1):(length(aboveThreshWithMV)+ncol(HighLowAbundWithMV)*nrow(HighLowAbundWithMV))]
  HighMedAbundWithMV = remainingValuesVector[(length(aboveThreshWithMV)+(length(HighLowAbundWithMV)+1)):length(remainingValuesVector)]

  aboveThreshWithMV = matrix(aboveThreshWithMV,
                             nrow = size_aboveThreshWithMV[1],
                             ncol = size_aboveThreshWithMV[2],
                             byrow = TRUE)
  HighLowAbundWithMV = matrix(HighLowAbundWithMV,
                              nrow = size_HighLowAbundWithMV[1],
                              ncol = size_HighLowAbundWithMV[2],
                              byrow = TRUE)
  HighMedAbundWithMV = matrix(HighMedAbundWithMV,
                              nrow = size_HighMedAbundWithMV[1],
                              ncol = size_HighMedAbundWithMV[2],
                              byrow = TRUE)

  if(dz==0){  #no missing in medium AAG-mtb
    MedAbundWithMV = HighMedAbundWithMV
  }else{
    MedAbundWithMV = cbind(LowMedAbund, HighMedAbundWithMV)
  }

  if (ncol(MedAbundWithMV) > 1){
    for (i in seq_len(nrow(MedAbundWithMV))){
      MedAbundWithMV[i,] = MedAbundWithMV[i,][order(sortAbundMedAbundidx[i,])]
    }}

  if(dw==0){ #no missing in low AAG-mtb
    LowAbundWithMV = HighLowAbundWithMV
  }else{
    
    

    LowAbundWithMV = cbind(LowLowAbund, HighLowAbundWithMV)
  }
  if (ncol(LowAbundWithMV)>1){
    for (i in seq_len(nrow(LowAbundWithMV))){
      LowAbundWithMV[i,] = LowAbundWithMV[i,][order(sortAbundLowAbundidx[i,])]
    }}

  if (dz==0 & dw == 0){
    belowThreshWithMV = cbind(LowAbundWithMV, MedAbundWithMV)
  } else if (dz==0 & dw ==1){
    belowThreshWithMV = rbind(LowAbundWithMV, t(MedAbundWithMV))
  }else{
    belowThreshWithMV = rbind(LowAbundWithMV, MedAbundWithMV)
  }
  if (dz==0 & dw==0){
    dataMVtemp = rbind(t(belowThreshWithMV), aboveThreshWithMV)
  }else{
    dataMVtemp = rbind(belowThreshWithMV, aboveThreshWithMV)
  }

  dataMV = dataMVtemp[order(sortAvgIdx),]
  dataMV[is.na(dataMV)] = 0

  if(!label){
    dataMV[dataMV == "MCAR"] = NA
    dataMV[dataMV == "MNAR"] = NA
    dataMV <- matrix(as.numeric(dataMV), nrow = nrow(dataMV), ncol = ncol(dataMV))
  }

  return(dataMV)
}

# Altered NSKNN mixed missingness function
removeDataMM_altered = function(data,
                                percentMiss,
                                percentBelowThresh_I,
                                LowAbundThresh_II,
                                percentMisslowAbund_III){

  #"""
 # https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04659-1

 # """


  data = as.matrix(data)

  #safe zero pstions 
  zero_mask <- data == 0 
  zero_safe_data <- data
  zero_safe_data[zero_mask] <- NA 

  


  percentMisslowAbund_III = percentMisslowAbund_III/0.9
  percentMVmedAbund = 0.5*percentMisslowAbund_III
  TotalNumValRemove = round(percentMiss/100*ncol(data)*nrow(data))
 

  #zero safe mean sort
  zero_safe_row_mean <- rowMeans(zero_safe_data, na.rm = TRUE )
  temp = sort(zero_safe_row_mean, decreasing = FALSE, index = TRUE)
  avgMetabolite = temp[["x"]]
  sortAvgIdx = temp[["ix"]]

  numMetBelowThresh = round(percentBelowThresh_I/100*length(avgMetabolite))
  belowThreshMet = data[sortAvgIdx[seq_len(numMetBelowThresh)],, drop = FALSE]
  belowThreshIdx = sortAvgIdx[seq_len(numMetBelowThresh)]
  aboveThreshMet = data[sortAvgIdx[(numMetBelowThresh+1):length(avgMetabolite)],, drop = FALSE]
  aboveThreshIdx = sortAvgIdx[(numMetBelowThresh+1):length(avgMetabolite)]

  numMetLowAbund = round(LowAbundThresh_II/100*nrow(belowThreshMet))
  lowAbundMet = as.matrix(data[sortAvgIdx[seq_len(numMetLowAbund)],])
  medAbundmet = as.matrix(data[sortAvgIdx[(numMetLowAbund+1):nrow(belowThreshMet)],])
  
  #zero safe lowAbund values
  zero_safe_low <- lowAbundMet#alterd
  zero_safe_low[zero_mask[sortAvgIdx[seq_len(numMetLowAbund)],]] <- NA 

  #zero safe medAbund values 
  zero_safe_med <- medAbundmet#alterd
  zero_safe_med[zero_mask[sortAvgIdx[(numMetLowAbund+1):nrow(belowThreshMet)],]] <- NA

  # # of samples for each metabolite that belong to the low abundance group
  # of metabolites that will be missing.
  numSampBelowAbundThreshLowAbund = round(percentMisslowAbund_III/100*ncol(zero_safe_low))
  if (numSampBelowAbundThreshLowAbund == 0){dw = 0}else{dw=1}
  # Sort matrix with low mean abundance metabolites from low to high abundance
  # values for each metabolite.

  sort_only_nonzero <- function(x) {  #function to only sort non zero values 
    nz <- x
    nz[is.na(nz)] <- Inf  # so NAs go to end
    ord <- order(nz)
    x[ord]
  }

  sort_idx_only_nonzero <- function(x) {#function to only sort non zero indexes values 
    nz <- x
    nz[is.na(nz)] <- Inf
    order(nz)
  }

 
  
  #zero safe sort
  sortLowAbundmet = t(apply(zero_safe_low, 1, sort_only_nonzero))
  sortAbundLowAbundidx = t(apply(zero_safe_low, 1, sort_idx_only_nonzero))

  # Create matrix with low abundance values for each metabolite with a low mean
  # abundance.
 

  LowLowAbund = as.matrix(sortLowAbundmet[,seq_len(numSampBelowAbundThreshLowAbund), drop = FALSE])
  # Create a matrix with high abundance values for each metabolite with a low
  # mean abundance
  
  HighLowAbund = as.matrix(sortLowAbundmet[,(numSampBelowAbundThreshLowAbund+1):ncol(sortLowAbundmet), drop = FALSE])

  alwaysMV = round(0.8*ncol(LowLowAbund))
  LowLowAbund[,seq_len(alwaysMV)] = "MNAR"

  #low abundance MNAR values
  if (alwaysMV < ncol(LowLowAbund)){
    LowLowAbundRandomMV = as.matrix(
      LowLowAbund[,(alwaysMV+1):ncol(LowLowAbund), drop = FALSE]
      )
    LowLowAbund = as.matrix(LowLowAbund[,-((alwaysMV+1):ncol(LowLowAbund)), drop = FALSE])
    numVal = ncol(LowLowAbundRandomMV)*nrow(LowLowAbundRandomMV)
    
    numValRemove = round(0.5*numVal)

    safe_select <- which(!is.na(LowLowAbundRandomMV))#alterd
    removeIdx = sample(safe_select,numValRemove)#alterd

    temp = as.vector(t(LowLowAbundRandomMV))
    temp[removeIdx] = "MNAR"
    LowLowAbundRandomMV = matrix(temp, nrow = nrow(LowLowAbundRandomMV),
                                 ncol = ncol(LowLowAbundRandomMV), byrow = TRUE)
    LowLowAbund = cbind(LowLowAbund, LowLowAbundRandomMV)
  }

  # # of samples for each metabolite that belong to the medium abundance group
  # of metabolites that will be missing.
  numSampBelowAbundThreshMedAbund = round(percentMVmedAbund/100*ncol(zero_safe_med))
  if (numSampBelowAbundThreshMedAbund == 0){dz = 0}else{dz=1}
  # Sort matrix with high mean abundance metabolites from low to high abundance
  # values for each metabolite.

  #zero safe sort
  sortmedAbundmet = t(apply(zero_safe_med,1, sort_only_nonzero))#alterd
  sortAbundMedAbundidx = t(apply(zero_safe_med,1, sort_idx_only_nonzero))#alterd


  # Create matrix with low abundance values for each metabolite with a high mean
  # abundance.
  LowMedAbund = as.matrix(sortmedAbundmet[,seq_len(numSampBelowAbundThreshMedAbund), drop = FALSE])
  # Create a matrix with high abundance values for each metabolite with a high
  # mean abundance.
  HighMedAbund = as.matrix(sortmedAbundmet[,(numSampBelowAbundThreshMedAbund+1):ncol(sortmedAbundmet), drop = FALSE])
  alwaysMV = round(0.8*ncol(LowMedAbund))
  LowMedAbund[,seq_len(alwaysMV)] = "MNAR"

  #med abundance MNAR values
  if (alwaysMV < ncol(LowMedAbund)){
    LowMedAbundRandomMV = as.matrix(
      LowMedAbund[,(alwaysMV+1):ncol(LowMedAbund), drop = FALSE]
      )
    LowMedAbund = as.matrix(LowMedAbund[,-((alwaysMV+1):ncol(LowMedAbund)), drop = FALSE])
    numVal = ncol(LowMedAbundRandomMV)*nrow(LowMedAbundRandomMV)
    numValRemove = round(0.5*numVal)

    safe_select <- which(!is.na(LowMedAbundRandomMV))#alterd
    removeIdx = sample(safe_select,numValRemove)#alterd

    
    temp = as.vector(t(LowMedAbundRandomMV))
    temp[removeIdx] = "MNAR"
    LowMedAbundRandomMV = matrix(temp,
                                 nrow = nrow(LowMedAbundRandomMV),
                                 ncol = ncol(LowMedAbundRandomMV),
                                 byrow = TRUE)
    LowMedAbund = cbind(LowMedAbund, LowMedAbundRandomMV)
  }

  # Add MCAR values
  numMVbelowThreshLOD = (sum(sum(is.na(LowMedAbund))) +
                           sum(sum(is.na(LowLowAbund))))*dz
  numMVrandom = TotalNumValRemove - numMVbelowThreshLOD
  if (numMVrandom < 0){
    numMVrandom = 0
  }

  MVrandomAvailable = ncol(data)*nrow(data)-dz*ncol(LowMedAbund)*nrow(LowMedAbund)-dw*ncol(LowLowAbund)*nrow(LowLowAbund)
  aboveThreshWithMV = aboveThreshMet
  HighLowAbundWithMV = HighLowAbund
  HighMedAbundWithMV = HighMedAbund


  size_aboveThreshWithMV = dim(aboveThreshWithMV)
  size_HighLowAbundWithMV = dim(HighLowAbundWithMV)
  size_HighMedAbundWithMV = dim(HighMedAbundWithMV)


  #zero safe remaining value 
  remainingValuesVector = c(as.vector(t(aboveThreshWithMV)),as.vector(t(HighLowAbundWithMV)),as.vector(t(HighMedAbundWithMV)))

  print(remainingValuesVector)
  nonZeroIdx <- which(remainingValuesVector != 0)
  
  #MCAR valeues in Remaining values
  removeIdx = sample(nonZeroIdx, numMVrandom)
  remainingValuesVector[removeIdx] = "MCAR"


  
  aboveThreshWithMV = remainingValuesVector[seq_len((ncol(aboveThreshWithMV)*nrow(aboveThreshWithMV)))]
  HighLowAbundWithMV = remainingValuesVector[(length(aboveThreshWithMV)+1):(length(aboveThreshWithMV)+ncol(HighLowAbundWithMV)*nrow(HighLowAbundWithMV))]
  HighMedAbundWithMV = remainingValuesVector[(length(aboveThreshWithMV)+(length(HighLowAbundWithMV)+1)):length(remainingValuesVector)]

 


  
  dataMV[is.na(dataMV)] = "MCAR"

  
  return(dataMV)
}

removeDataMM = function(data,
                        percentMiss,
                        percentBelowThresh_I,
                        LowAbundThresh_II,
                        percentMisslowAbund_III){

 # """
 # https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04659-1
  
 # """
  data = as.matrix(data)

  #safe zero pstions 
  zero_mask <- data == 0 
  zero_safe_data <- data
  zero_safe_data[zero_mask] <- NA 

  


  percentMisslowAbund_III = percentMisslowAbund_III/0.9
  percentMVmedAbund = 0.5*percentMisslowAbund_III
  TotalNumValRemove = round(percentMiss/100*ncol(data)*nrow(data))
 

  #zero safe mean sort
  zero_safe_row_mean <- rowMeans(zero_safe_data, na.rm = TRUE )
  temp = sort(zero_safe_row_mean, decreasing = FALSE, index = TRUE)
  avgMetabolite = temp[["x"]]
  sortAvgIdx = temp[["ix"]]

  numMetBelowThresh = round(percentBelowThresh_I/100*length(avgMetabolite))
  belowThreshMet = data[sortAvgIdx[seq_len(numMetBelowThresh)],, drop = FALSE]
  belowThreshIdx = sortAvgIdx[seq_len(numMetBelowThresh)]
  aboveThreshMet = data[sortAvgIdx[(numMetBelowThresh+1):length(avgMetabolite)],, drop = FALSE]
  aboveThreshIdx = sortAvgIdx[(numMetBelowThresh+1):length(avgMetabolite)]

  numMetLowAbund = round(LowAbundThresh_II/100*nrow(belowThreshMet))
  lowAbundMet = as.matrix(data[sortAvgIdx[seq_len(numMetLowAbund)],])
  medAbundmet = as.matrix(data[sortAvgIdx[(numMetLowAbund+1):nrow(belowThreshMet)],])
  
  #zero safe lowAbund values
  zero_safe_low <- lowAbundMet#alterd
  zero_safe_low[zero_mask[sortAvgIdx[seq_len(numMetLowAbund)],]] <- NA #alterd

  #zero safe medAbund values 
  zero_safe_med <- medAbundmet#alterd
  zero_safe_med[zero_mask[sortAvgIdx[(numMetLowAbund+1):nrow(belowThreshMet)],]] <- NA

  # # of samples for each metabolite that belong to the low abundance group
  # of metabolites that will be missing.
  numSampBelowAbundThreshLowAbund = round(percentMisslowAbund_III/100*ncol(zero_safe_low))
  if (numSampBelowAbundThreshLowAbund == 0){dw = 0}else{dw=1}
  # Sort matrix with low mean abundance metabolites from low to high abundance
  # values for each metabolite.

  sort_only_nonzero <- function(x) {  #function to only sort non zero values 
    nz <- x
    nz[is.na(nz)] <- Inf  # so NAs go to end
    ord <- order(nz)
    x[ord]
  }

  sort_idx_only_nonzero <- function(x) {#function to only sort non zero indexes values 
    nz <- x
    nz[is.na(nz)] <- Inf
    order(nz)
  }

 
  
  #zero safe sort
  sortLowAbundmet = t(apply(zero_safe_low, 1, sort_only_nonzero))
  sortAbundLowAbundidx = t(apply(zero_safe_low, 1, sort_idx_only_nonzero))

  # Create matrix with low abundance values for each metabolite with a low mean
  # abundance.
 

  LowLowAbund = as.matrix(sortLowAbundmet[,seq_len(numSampBelowAbundThreshLowAbund), drop = FALSE])
  # Create a matrix with high abundance values for each metabolite with a low
  # mean abundance
  
  HighLowAbund = as.matrix(sortLowAbundmet[,(numSampBelowAbundThreshLowAbund+1):ncol(sortLowAbundmet), drop = FALSE])

  alwaysMV = round(0.8*ncol(LowLowAbund))
  LowLowAbund[,seq_len(alwaysMV)] = "MNAR"

  #low abundance MNAR values
  if (alwaysMV < ncol(LowLowAbund)){
    LowLowAbundRandomMV = as.matrix(
      LowLowAbund[,(alwaysMV+1):ncol(LowLowAbund), drop = FALSE]
      )
    LowLowAbund = as.matrix(LowLowAbund[,-((alwaysMV+1):ncol(LowLowAbund)), drop = FALSE])
    numVal = ncol(LowLowAbundRandomMV)*nrow(LowLowAbundRandomMV)
    
    numValRemove = round(0.5*numVal)

    safe_select <- which(!is.na(LowLowAbundRandomMV))#alterd
    removeIdx = sample(safe_select,numValRemove)#alterd

    temp = as.vector(t(LowLowAbundRandomMV))
    temp[removeIdx] = "MNAR"
    LowLowAbundRandomMV = matrix(temp, nrow = nrow(LowLowAbundRandomMV),
                                 ncol = ncol(LowLowAbundRandomMV), byrow = TRUE)
    LowLowAbund = cbind(LowLowAbund, LowLowAbundRandomMV)
  }

  # # of samples for each metabolite that belong to the medium abundance group
  # of metabolites that will be missing.
  numSampBelowAbundThreshMedAbund = round(percentMVmedAbund/100*ncol(zero_safe_med))
  if (numSampBelowAbundThreshMedAbund == 0){dz = 0}else{dz=1}
  # Sort matrix with high mean abundance metabolites from low to high abundance
  # values for each metabolite.

  #zero safe sort
  sortmedAbundmet = t(apply(zero_safe_med,1, sort_only_nonzero))#alterd
  sortAbundMedAbundidx = t(apply(zero_safe_med,1, sort_idx_only_nonzero))#alterd


  # Create matrix with low abundance values for each metabolite with a high mean
  # abundance.
  LowMedAbund = as.matrix(sortmedAbundmet[,seq_len(numSampBelowAbundThreshMedAbund), drop = FALSE])
  # Create a matrix with high abundance values for each metabolite with a high
  # mean abundance.
  HighMedAbund = as.matrix(sortmedAbundmet[,(numSampBelowAbundThreshMedAbund+1):ncol(sortmedAbundmet), drop = FALSE])
  alwaysMV = round(0.8*ncol(LowMedAbund))
  LowMedAbund[,seq_len(alwaysMV)] = "MNAR"

  #med abundance MNAR values
  if (alwaysMV < ncol(LowMedAbund)){
    LowMedAbundRandomMV = as.matrix(
      LowMedAbund[,(alwaysMV+1):ncol(LowMedAbund), drop = FALSE]
      )
    LowMedAbund = as.matrix(LowMedAbund[,-((alwaysMV+1):ncol(LowMedAbund)), drop = FALSE])
    numVal = ncol(LowMedAbundRandomMV)*nrow(LowMedAbundRandomMV)
    numValRemove = round(0.5*numVal)

    safe_select <- which(!is.na(LowMedAbundRandomMV))#alterd
    removeIdx = sample(safe_select,numValRemove)#alterd

    
    temp = as.vector(t(LowMedAbundRandomMV))
    temp[removeIdx] = "MNAR"
    LowMedAbundRandomMV = matrix(temp,
                                 nrow = nrow(LowMedAbundRandomMV),
                                 ncol = ncol(LowMedAbundRandomMV),
                                 byrow = TRUE)
    LowMedAbund = cbind(LowMedAbund, LowMedAbundRandomMV)
  }

  # Add MCAR values
  numMVbelowThreshLOD = (sum(sum(is.na(LowMedAbund))) +
                           sum(sum(is.na(LowLowAbund))))*dz
  numMVrandom = TotalNumValRemove - numMVbelowThreshLOD
  if (numMVrandom < 0){
    numMVrandom = 0
  }

  MVrandomAvailable = ncol(data)*nrow(data)-dz*ncol(LowMedAbund)*nrow(LowMedAbund)-dw*ncol(LowLowAbund)*nrow(LowLowAbund)
  aboveThreshWithMV = aboveThreshMet
  HighLowAbundWithMV = HighLowAbund
  HighMedAbundWithMV = HighMedAbund
  size_aboveThreshWithMV = dim(aboveThreshWithMV)
  size_HighLowAbundWithMV = dim(HighLowAbundWithMV)
  size_HighMedAbundWithMV = dim(HighMedAbundWithMV)


  #zero safe remaining value 
  remainingValuesVector = c(as.vector(t(aboveThreshWithMV)),as.vector(t(HighLowAbundWithMV)),as.vector(t(HighMedAbundWithMV)))
  ramaining_Values_Zero_mask <- as.vector(t(zero_mask[c(rownames(HighLowAbund),rownames(HighMedAbund),aboveThreshIdx), ]))
  zero_save_ramaining_Values <- which(!ramaining_Values_Zero_mask)

  #MCAR valeues in Remaining values
  removeIdx = sample(zero_save_ramaining_Values, numMVrandom)
  remainingValuesVector[removeIdx] = "MCAR"


  #rejoin matrix 
  medWithMV <- cbind(LowMedAbund, HighMedAbund)
  lowWithMV <- cbind(LowLowAbund, HighLowAbund)
  belowThreshWithMV <- rbind(lowWithMV, medWithMV)
  num_cols <- ncol(data)
  num_rows_above <- nrow(aboveThreshMet)
  above <- matrix(remainingValuesVector,
                  nrow=num_rows_above,
                  ncol = num_cols,
                  byrow=TRUE)

  

  dataMVtemp <- rbind(belowThreshWithMV, above)
  dataMV <- dataMVtemp[order(sortAvgIdx), ]
  
  #reinsert missing values
  dataMV[is.na(dataMV)] <- 0  #rebuild zeros
  dataMV[dataMV == "MNAR"] <- NA
   dataMV[dataMV == "MCAR"] <- NA
  
  return(dataMV)
}

nsKNN = function(data_miss,
                 k,
                 iters=1,
                 weighted = FALSE,
                 scale = FALSE,
                 shuffle = TRUE){

  #"""
 # https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04659-1
  
  #"""
  if (scale){
    data_miss_scale = scale(data_miss)
    data_miss = data_miss_scale
  }


  findKneighbors = function(distances, k, col){
    distances = distances[distances[,1]==col | distances[,3]==col, ]
    sortedDistances = distances[order(distances[,2]),]
    neighbors = sortedDistances[seq_len(k), ]
    n = list()
    for(i in seq_len(nrow(neighbors))){n[i] = ifelse(neighbors[i,1] == col,
                                              list(neighbors[i, c(3,2)]),
                                              list(neighbors[i, c(1,2)]))}
    n = matrix(unlist(n), nrow = k, ncol = 2, byrow = TRUE)
    return(n)
  }



  distances = as.matrix(dist(t(data_miss))/sqrt(nrow(data_miss)))
  distances = distances[lower.tri(distances)]
  distances = data.frame(col1 = combn(paste(seq_len(ncol(data_miss)),
                                            sep = ""), 2)[1,],
                         distance = distances,
                         col2 = combn(paste(seq_len(ncol(data_miss)),
                                            sep = ""), 2)[2,])
  distances = apply(distances, 2, as.numeric)

  minVals = apply(data_miss, 1, min, na.rm=TRUE)

  data_imputed = data_miss

  for (col in seq_len(ncol(data_miss))){
    neighbors = findKneighbors(distances, k, col)
    missingNum = which(is.na(data_miss[,col]))

    if(length(missingNum)==0){
      next
    }

    neighborVals = data_miss[missingNum, c(col, neighbors[,1])]
    mins = minVals[missingNum]
    if (is.null(dim(neighborVals))){
      listns = list(c(neighborVals[-1], mins))
      listns = lapply(listns, function(x){ifelse(is.na(x), x[k+1], x)})

      filledNeighborVals = unlist(listns)[-length(neighborVals)]
    } else {
      listns = asplit(cbind(neighborVals[,-1], mins), 1)
      listns = lapply(listns, function(x){ifelse(is.na(x), x[k+1], x)})


    filledNeighborVals = cbind(matrix(unlist(listns),
                                      ncol = ncol(neighborVals),
                                      nrow = nrow(neighborVals),
                                      byrow = TRUE))[,-ncol(neighborVals)]
    }

    weightMultiplier = (1/neighbors[,2])/sum(1/neighbors[,2])

    if (weighted){
      if (is.null(dim(neighborVals))){
        data_imputed[missingNum, col] = sum(
          filledNeighborVals%*%diag(weightMultiplier),
          na.rm=TRUE)
      } else {
        data_imputed[missingNum, col] = apply(
          filledNeighborVals%*%diag(weightMultiplier),
          1, sum, na.rm=TRUE)
      }

    } else {
      if (is.null(dim(neighborVals))){
        data_imputed[missingNum, col] = mean(
          filledNeighborVals%*%diag(weightMultiplier),
          na.rm=TRUE)
      }else{
        data_imputed[missingNum, col] = apply(filledNeighborVals,
                                              1, mean, na.rm=TRUE)
      }

    }

  }

  if (iters < 2){
    if (scale){
      data_imputed = t(apply(data_imputed, 1, function(r){
        r*attr(data_miss_scale,'scaled:scale') +
          attr(data_miss_scale, 'scaled:center')
        }))
      return(data_imputed)
    } else {
      return(data_imputed)
    }

  } else {
    for (j in seq_len((iters-1))){
      if (shuffle){
        numCols = sample(seq_len(ncol(data_imputed)),
                         ncol(data_imputed), replace = FALSE)
      } else {
        numCols = seq_len(ncol(data_imputed))
      }

      for (numCol in numCols){
        neighborCols = numCols[numCols!=numCol]
        distances = as.matrix(dist(t(data_imputed[,c(numCol, neighborCols)])))
        distances = cbind(neighborCols, distances[1,][-1])
        sortDistances = distances[order(distances[,2]),]
        NN = sortDistances[seq_len(k)]
        missingNum = which(is.na(data_miss[,numCol]))
        neighborVals = data_imputed[missingNum, c(NN)]
        weightMultiplier = (1/sortDistances[,2][seq_len(k)])/sum(1/sortDistances[,2][seq_len(k)])
        if (weighted){
          if (is.null(dim(neighborVals))){
            data_imputed[missingNum, numCol] = sum(
              as.matrix(neighborVals)%*%diag(weightMultiplier),
              na.rm=TRUE)
          } else {
            data_imputed[missingNum, numCol] = apply(
              as.matrix(neighborVals)%*%diag(weightMultiplier),
              1, sum, na.rm=TRUE)
          }

        } else {
          if (is.null(dim(neighborVals))){
            data_imputed[missingNum, numCol] = mean(neighborVals,
                                                    na.rm=TRUE)
          } else {
            data_imputed[missingNum, numCol] = apply(neighborVals,
                                                     1, mean, na.rm=TRUE)
          }

        }

      }

    }
   if (scale){
     data_imputed = t(apply(data_imputed,
                            1, function(r){
                              r*attr(data_miss_scale,'scaled:scale') +
                                attr(data_miss_scale, 'scaled:center')
                              }))
     return(data_imputed)
   } else {
     return(data_imputed)
   }

  }

}

## SINGLE IMPUTATION functions
##
## input: missData - a data matrix with NAs for missingData
## output: output - the data matrix with imputed missing values
##X|censored or X|censored, Z
##Kap impute 2 uses density function

kapImpute2 <- function(missData){
  #"""
  #https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04659-1
  
  #"""
  #missData <- as.matrix(missData)
  # Count missing values per variable
  ROWna <- rowSums(is.na(missData))
  # Sort row indices by missingness
  sortRowNA <- sort(ROWna, index = TRUE)
  # Find row with minimum missingness
  MissVecStart <- (min(which(sortRowNA$x != 0)))

  nSubj <- dim(missData)[2]
  nVar <- dim(missData)[1]

  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

  output <- missData
  for (i in MissVecStart:nVar){
    # NOTE: index of ordered minimum count missing is i

    # to get row with respective missing count
    index <- sortRowNA$ix[i]
    rowMissingCount <- sortRowNA$x[i]

    # column indices for NA in row of consideration
    missSubj <- which(is.na(output[index,]))

    # estimating censoring point as min of metabolite
    c <- min(output[index,], na.rm = TRUE)

    for (j in seq_len(length(missSubj))){
      S.est <- 1 - (rowMissingCount/(nVar -rowMissingCount))
      #print(S.est)
      meanEst <- rowMeans(output, na.rm = TRUE)[index]
      sdEst <- sd(output[index,], na.rm = TRUE)
      #print("new line mean,sd,c")
      #print(meanEst)
      #print(sdEst)
      #print(c)

      # E[x_j | x_j < c]
      numerator <- (meanEst/2)*erf((c-meanEst)/(sqrt(2)*sdEst)) + (meanEst/2) -
        (sdEst/sqrt(2*pi))*exp(-((c-meanEst)^2)/(2*sdEst^2))
      denominator <- pnorm(c, meanEst, sdEst)
      #print(denominator)
      output[index, missSubj[j]] <- numerator/denominator

      rowMissingCount = rowMissingCount - 1
    }
  }
  return(output)
}

# Function to estimate thresholds
check_distance = function(data_miss, data_sub, PercentMiss){

  #"""
  #https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04659-1
  #
  #"""

  euc.dist = function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
  ## Approximate thresholds
  # Find highest to low average metabolite
  avgMetabols = sort(apply(data_miss, 1, mean, na.rm=TRUE),
                     decreasing = TRUE, index = TRUE)
  # Sort vertically high to low metabolite
  rowSortedData = data_miss[avgMetabols$ix,]
  fullSortedData = split(rowSortedData, seq(nrow(rowSortedData)))

  fullSortedData = lapply(fullSortedData, function(x){
    sum(is.na(x)/ncol(data_miss))# Count missing vals per metabolite
  })
  fullSortedData = as.data.frame(do.call(rbind, fullSortedData))

  if (PercentMiss < 5){
    thresh_I = c(5)
  } else {
    thresh_I = seq(5,PercentMiss,5)
  }

  thresh_II = seq(60,80,5)
  thresh_III = seq(5, 60, 5)
  grid = base::expand.grid(thresh_I, thresh_II, thresh_III)
  threshs = list()
  distances = numeric()
  try({
    for (i in seq_len(nrow(grid))) {
      alpha = grid$Var1[i]
      beta = grid$Var2[i]
      gamma = grid$Var3[i]
      ## Approximate thresholds
      data_test = suppressWarnings({
        tryCatch({
          removeDataMM_Zero_safe(data_sub, percentMiss = PercentMiss, alpha, beta, gamma,label = FALSE)
        }, error = function(e) e)
      })

      if(inherits(data_test, "error")) next
      # Find highest to low average metabolite
      data_avgMetabols = sort(apply(data_test, 1, mean, na.rm=TRUE),
                              decreasing = TRUE, index = TRUE)
      # Sort vertically high to low metabolite
      data_rowSortedData = data_test[data_avgMetabols$ix,]
      data_fullSortedData = split(data_rowSortedData,
                                  seq(nrow(data_rowSortedData)))

      data_fullSortedData = lapply(data_fullSortedData, function(x){
        sum(is.na(x)/ncol(data_test))# Count missing vals per metabolite
      })
      data_fullSortedData = as.data.frame(do.call(rbind, data_fullSortedData))
      distance = euc.dist(data_fullSortedData, fullSortedData)
      threshs[[i]] = cbind(alpha, beta, gamma) # Store current threshold I
      distances[i] = distance # Store distance computed

    }
  }, silent = TRUE)

  return(list(thresh = threshs[which.min(distances)],
              distance = distances[which.min(distances)]))
}

# function to get largest complete data subset
largest_complete_subset = function(data_miss){

  #"""
  #https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04659-1
  
  #"""
  shuffle_rows = function(x) {
    x = x[, sample(ncol(x))]
    return(x)
  }
  moveNAs = function(x) {
    # count NA
    num.na = sum(is.na(x))
    # remove NA
    x = x[!is.na(x)]
    # glue the number of NAs at the end
    x = c(x, rep(NA, num.na))
    return(x)
  }
  # Get largest possible subset of complete data
  data_original = shuffle_rows(data_miss)
  data_original = split(data_original, seq(nrow(data_miss)))
  data_original = lapply(data_original, moveNAs)

  SortedData = lapply(data_original, function(x){
    item = unlist(x)
    item = item[order(item)]
    return(item)

  })
  SortedData = do.call(rbind, SortedData)
  mostNAs = lapply(data_original, function(x){
    sum(is.na(x))
  })
  mostNAs.ix = which(unlist(mostNAs)==max(unlist(mostNAs)), arr.ind = TRUE)[1]
  data_original = do.call(rbind, data_original)
  remove.ix = ncol(data_original)-mostNAs[[mostNAs.ix]]
  data = data_original[,seq_len(remove.ix)]
  return(data)
}

# function to generate predictors/labels
generate_predictors = function(missingData, labels=TRUE){

 # """
 # https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04659-1
  
 # """
  # Normalization function
  normaliz = function(x, eps=0.001){
    return((x-(min(x, na.rm = TRUE)-eps))/((max(x, na.rm = TRUE)-min(x, na.rm = TRUE))+eps))
    }

  numcols = ncol(missingData)
  numrows = nrow(missingData)
  # calc max_metabolite predictor
  max_metabolite = apply(apply(missingData, 2, as.numeric), 1, max, na.rm=TRUE)
  # calc min_metabolite predictor
  min_metabolite = apply(apply(missingData, 2, as.numeric), 1, min, na.rm=TRUE)
  # calc median_metabolite predictor
  median_metabolite = apply(apply(missingData, 2, as.numeric),
                            1, median, na.rm=TRUE)
  # calc average_metabolite predictor
  avg_metabolite = apply(apply(missingData, 2, as.numeric),
                         1, mean, na.rm=TRUE)
  # calc metabolite quantiles
  metabolite_quantiles = apply(apply(missingData, 2, as.numeric),
                               1, quantile, na.rm=TRUE)
  # calc amount missing per metabolite
  num_missing = apply(apply(missingData, 2, as.numeric),
                      1, function(x){sum(is.na(x))})
  # calc percent missing per metabolite predictor
  num_missing = num_missing/ncol(missingData)

  levels = matrix(NA, ncol = numcols, nrow = numrows) # initalize matrix
  # Use normalized data
  missingData_norm = t(apply(apply(missingData, 2, as.numeric), 1, normaliz))
  for (i in seq_len(numrows)){ # calc metabolite quantile categories predictor
    for (j in seq_len(numcols)){
      levels[i, j] = ifelse(missingData_norm[i,j] > metabolite_quantiles[3, i], "high",
                            ifelse(missingData_norm[i,j]<metabolite_quantiles[2, i], "low",
                                   ifelse(is.na(missingData_norm[i,j]), missingData[i,j],"medium")))
    }
  }
  levels[is.na(levels)] = "none"
  frame.as.vector = as.vector(t(missingData)) # Vectorize matrix

  # append max_metabolite predictor
  frame.as.vector = cbind(frame.as.vector, rep(max_metabolite, each=numcols))
  # append min_metabolite predictor
  frame.as.vector = cbind(frame.as.vector, rep(min_metabolite, each=numcols))
  # append median_metabolite predictor
  frame.as.vector = cbind(frame.as.vector, rep(median_metabolite, each=numcols))
  # append average_metabolite predictor
  frame.as.vector = cbind(frame.as.vector, rep(avg_metabolite, each=numcols))
  # append percent missing per metabolite predictor
  frame.as.vector = cbind(frame.as.vector, rep(num_missing, each=numcols))
  # append metabolite quantile categories predictor
  frame.as.vector = cbind(frame.as.vector, as.vector(t(levels)))

  if (labels){
    target = ifelse(grepl("^[A-Za-z\\/]+$", frame.as.vector[,1], perl = TRUE),
                    frame.as.vector[,1], 'O') # Generate label vector
    frame.as.vector = cbind(frame.as.vector, target) # append labels
  }

  # normalize metabolite abundances
  missingvals = as.numeric(frame.as.vector[,1])
  # replace NAs with 0s
  missingvals[is.na(missingvals)] = 0
  # Append replacement to data frame
  frame.as.vector[,1] = missingvals
  # data structure change
  frame.as.vector = as.data.frame(frame.as.vector)
  # Ensure numeric predictors are numeric type
  frame.as.vector[,seq_len(6)] = apply(frame.as.vector[,seq_len(6)],
                                           2, as.numeric)
  # in case any rows are still missing omit
  frame.as.vector = na.omit(frame.as.vector)
  return(frame.as.vector)
}

# function to impute data usig specific imputation algorithms
imputation_algorithms = function(missingData,
                                 predictions,
                                 MCAR_algorithm,
                                 MNAR_algorithm,
                                 n_cores){

  #"""
  #https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04659-1
  
  #"""
  # Imputation step
  numcols = ncol(missingData)
  numrows = nrow(missingData)
  predicted_missingData = matrix(predictions,
                                 ncol = numcols, nrow = numrows, byrow = TRUE)
  MCAR_index = which(predicted_missingData == "MCAR", arr.ind = TRUE)
  MNAR_index = which(predicted_missingData == "MNAR", arr.ind = TRUE)
  non_missing = which(predicted_missingData == "O", arr.ind = TRUE)
  Imputed_data = matrix(NA, ncol = numcols, nrow = numrows)

  if (n_cores != 1){
    plan(multisession, workers = 2)

  SingleAlgorithmImputations = future_lapply(list(MNAR_algorithm,
                                                  MCAR_algorithm),function(l){

    if (l == "Single"){
      MNAR_imputation = kapImpute2(missingData)
    }

    if (l == 'nsKNN'){
      k = floor(sqrt(nrow(missingData)))
      if (k>numcols){
        k=numcols-1
      }
      MNAR_imputation = nsKNN(missingData, k, iters=1)
    }

    if (l == "BPCA"){
      MCAR_imputation = pca(missingData,
                            method = "bpca",
                            verbose=FALSE)@completeObs
    }
    if (l == "random_forest"){
      MCAR_imputation = missForest(missingData)[["ximp"]]
    }
    if (l == 'Multi_nsKNN'){
      k = floor(sqrt(nrow(missingData)))
      if (k>numcols){
        k=numcols-1
      }
      MCAR_imputation = nsKNN(missingData, k,
                              iters=4, weighted=TRUE, scale = TRUE,
                              shuffle = TRUE)
    }

    if (exists('MNAR_imputation')){
      return(MNAR_imputation = MNAR_imputation)
    } else {
      return(MCAR_imputation = MCAR_imputation)
    }
  }, future.globals = c("kapImpute2", "nsKNN"), future.seed=TRUE)

  plan(sequential)


  MNAR_imputation = SingleAlgorithmImputations[[1]]
  MCAR_imputation = SingleAlgorithmImputations[[2]]
  } else {
    if (MNAR_algorithm == "Single"){
      MNAR_imputation = kapImpute2(missingData)
    }

    if (MNAR_algorithm == 'nsKNN'){
      k = floor(sqrt(nrow(missingData)))
      if (k>numcols){
        k=numcols-1
      }
      MNAR_imputation = nsKNN(missingData, k, iters=1)
    }

    if (MCAR_algorithm == "BPCA"){
      MCAR_imputation = pca(missingData,
                            method = "bpca",
                            verbose=FALSE)@completeObs
    }
    if (MCAR_algorithm == "random_forest"){
      MCAR_imputation = missForest(missingData)[["ximp"]]
    }
    if (MCAR_algorithm == 'Multi_nsKNN'){
      k = floor(sqrt(nrow(missingData)))
      if (k>numcols){
        k=numcols-1
      }
      MCAR_imputation = nsKNN(missingData, k,
                              iters=4, weighted=TRUE, scale = TRUE,
                              shuffle = TRUE)
    }
  }

  Imputed_data[MCAR_index] = MCAR_imputation[MCAR_index]
  Imputed_data[MNAR_index] = MNAR_imputation[MNAR_index]
  Imputed_data[is.na(Imputed_data)] = missingData[non_missing]
  return(list(MAI = Imputed_data,
              MCAR_imps = as.matrix(MCAR_imputation),
              MNAR_imps = as.matrix(MNAR_imputation)))
}
modified_MAI = function(data_miss,
               MCAR_algorithm = c("BPCA", "Multi_nsKNN", "random_forest"),
               MNAR_algorithm = c('nsKNN', "Single"),
               n_cores = 1,
               assay_ix = 1,
               forest_list_args = list(
                   ntree = 300,
                   proximity = FALSE
               ),
               verbose = TRUE, SKIP_ZERO = FALSE){


  #"""
  #-altered MAI function to give out the matrix of classified Missing values / complete dataset / imposed Dataset
  #https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04659-1
  
  #"""
  
  check <- TRUE
  while(check){

    if(!SKIP_ZERO){
      check <-FALSE
    }


    MCAR_algorithm = match.arg(MCAR_algorithm)
    MNAR_algorithm = match.arg(MNAR_algorithm)

    #library(MAI)
    
    #set.seed(42) #Estimating the missingness patterns involves randomly imposing missing values in a complete dataset #setting seed for reproducibility 






    # Possible imputation algorithms
    MCAR_algorithms = c("BPCA", "random_forest", "Multi_nsKNN")
    MNAR_algorithms = c("nsKNN", "Single")

    # Check argument inputs
    if (!MCAR_algorithm %in% MCAR_algorithms){
      stop("Incorrect MCAR algorithm selected. Possible choices include:BPCA,
          random_forest, or Multi_nsKNN ")
    }

    if (!MNAR_algorithm %in% MNAR_algorithms){
      stop("Incorrect MNAR algorithm selected. Possible choices include: nsKNN
          or Single ")
    }

    # Check data
    if (is(data_miss, "SummarizedExperiment")){
      if (sum(is.na(assay(data_miss))) == 0){
        stop("No missing values detected to impute")
      }
      if (any(rowSums(is.na(assay(data_miss))) == ncol(data_miss))){
        stop("Detected a row that contains no observed data, omit this row.")
      }
      if (any(rowSums(is.na(assay(data_miss))) > 0.8*ncol(data_miss))){
        warning("Detected a row with more than 80% missing values.
              Consider omitting rows that are > 80% missing.")
      }
      if (ncol(data_miss) < 50 | nrow(data_miss) < 50){
        warning("The smallest data set MAI has been tested on was a 50x50 matrix.Accuracy can not be guaranteed for smaller data sets.")
      }
    }else{

    if (sum(is.na(data_miss)) == 0){
      stop("No missing values detected to impute")
    }
    if (any(rowSums(is.na(data_miss)) == ncol(data_miss))){
      stop("Detected a row that contains no observed data, omit this row.")
    }
    if (any(rowSums(is.na(data_miss)) > 0.8*ncol(data_miss))){
      warning("Detected a row with more than 80% missing values.
              Consider omitting rows that are > 80% missing.")
    }
    if (ncol(data_miss) < 50 | nrow(data_miss) < 50){
      #warning("The smallest data set MAI has been tested on was a 50x50 matrix.Accuracy can not be guaranteed for smaller data sets.")
      }
    }

    # Check n_cores
    if (n_cores%%1 != 0){
      stop("n_cores must be an integer greater than or equal to -1.")
    }
    if (n_cores > detectCores()){
      stop("You specified more cores than is available to your system.")
    }
    if (n_cores < -1){
      stop("n_cores must be an integer greater than or equal to -1.")
    }
    if (n_cores == 0){
      stop("n_cores must be a non-zero integer.")
    }

    # Check summarized experiment
    if (is(data_miss, "SummarizedExperiment")) {
      fldata = data_miss
      data_miss = assay(data_miss)
    }


    # Original data percent missing (total)
    PercentMiss = (sum(is.na(data_miss))/(nrow(data_miss)*ncol(data_miss)))*100
    data = largest_complete_subset(data_miss)
    # Estimate threshold I
    threshs = list()
    data_export <- data
    if (verbose){
      message("Estimating pattern of missingness")
    }

    for (i in seq_len(10)){
      try({threshs[[i]] = check_distance(data_miss, data, PercentMiss)},
          silent = TRUE)
    }

    threshs = threshs[!sapply(threshs,is.null)]


    tryCatch(
      {smallestDistance = threshs[[which.min(lapply(threshs, function(x){unlist(x[["distance"]])}))]][["thresh"]][[1]]},
      error=function(e) {
        message("Error message:")
        message(e)
        # Choose a return value in case of error
        return(stop("Not enough data to estimate the pattern of missingness!"))
        
      }
    )

    # Set Threshold params
    alpha = smallestDistance[1]
    beta = smallestDistance[2]
    gamma = smallestDistance[3]


    # generate missingness
    if (verbose){
      message('Imposing missingness')
    }

    missingData = removeDataMM_Zero_safe(data,
                                      percentMiss = PercentMiss,
                                      percentBelowThresh_I = alpha,
                                      LowAbundThresh_II = beta,
                                      percentMisslowAbund_III = gamma, label = TRUE)

    # generate predictors
    if (verbose){
      message('Generating features')
    }

    missing_data_export <- missingData


    #print(data_export)
    #print(missing_data_export)

    
    skip_zero_missing_data <- replace(missing_data_export, missing_data_export == "MCAR",NA)
    skip_zero_missing_data <- replace(skip_zero_missing_data, skip_zero_missing_data == "MNAR",NA)

    skip_zero_positions = which(is.na(skip_zero_missing_data), arr.ind = TRUE)

    skip_zero_values <- data_export[skip_zero_positions]
    #print(skip_zero_values)

    if(!(0 %in% skip_zero_values)){
      check <- FALSE
    }


    


  }






  full_data = suppressWarnings(generate_predictors(missingData, labels=TRUE))

  # Split into training and testing sets
  trainIndex = createDataPartition(full_data$target,
                                   p=0.95,
                                   list = FALSE,
                                   times = 1)
  trainingSet = full_data[trainIndex,]
  testingSet = full_data[-trainIndex,]
  # Train model
  fitControl = trainControl(
    method = "cv",
    number = 2)
  # Parallelize training
  if (n_cores > 1)
  {
    cores = n_cores
    cl = makePSOCKcluster(cores)
    registerDoParallel(cl)
  }
  if (n_cores == -1)
  {
    cores = detectCores()
    cl = makePSOCKcluster(cores)
    registerDoParallel(cl)
  }

  if (verbose){
    message("Training")
  }

  # Random forest
  rfFit = suppressMessages(
      do.call(
          train,
          c(
              list(
                  form = target ~ .,
                  data = trainingSet,
                  method = "rf",
                  trControl = fitControl,
                  verbose = FALSE
              ),
              forest_list_args
          )
      )
  )

  if (n_cores != 1){
    stopCluster(cl)
    registerDoSEQ()
  }

  # Use model to make predictions on initial data
  if (verbose){
    message("Predicting")
  }

  predictions = predict(rfFit,
                        newdata = suppressWarnings(
                          generate_predictors(data_miss, labels=FALSE))
                        )



    ############################Edite


  if (FALSE) {
  
  
    # Impute based on classification predictions
    if (verbose){
      message("Imputing")
    }
    Imputed_data = imputation_algorithms(data_miss,
                                        predictions,
                                        MCAR_algorithm = MCAR_algorithm,
                                        MNAR_algorithm = MNAR_algorithm,
                                        n_cores = n_cores)

    if(exists("fldata")){
      assay(fldata, assay_ix, withDimnames = FALSE) = as.matrix(Imputed_data[[1]])
      MCAR_imputations = Imputed_data[[2]]
      MNAR_imputations = Imputed_data[[3]]
      estimated_Params = list(
        Alpha = alpha,
        Beta = beta,
        Gamma = gamma
      )

      metadata(fldata) = c(metadata(fldata),
                          list(list(
                            Estimated_Params = estimated_Params
                          ))
                              )
      names(metadata(fldata))[[length(metadata(fldata))]] =
        paste("meta_assay", assay_ix, sep = "_")

      return(fldata)
    }else{
      return(list(Imputed_data = Imputed_data[["MAI"]],
                  Estimated_Params = list(
                    Alpha = alpha,
                    Beta = beta,
                    Gamma = gamma
                  )))
    }
  }

  return (list(predictions = predictions,complete_data = data_export,complete_missing_data = missing_data_export))
  
}



imputate = function(data){

   # """
   # -Imputating a given dataset with MAI classification and selcted methods for MCAR and MNAR
    #"""
    
    data_check <- TRUE
    while(data_check){

        data_ <- data
   
        #load data as datamatrix  "data/imputation_format_data.csv"
        result <- modified_MAI(data_miss = data_ , MCAR_algorithm = "BPCA", MNAR_algorithm = "Single",forest_list_args = list(ntree = 300,proximity = FALSE),verbose = FALSE,SKIP_ZERO = FALSE)
        #using modified MAI for detcting Missing value type

        predicted_missingData = matrix(result$predictions,ncol = ncol(data_), nrow = nrow(data_), byrow = TRUE)

        mask <- predicted_missingData != "O"
        data_before <- data_
        data_[mask] <- predicted_missingData[mask] # merging dataset with predicted missingness

        

        if(!any(is.na(data_))){
            data_check <- FALSE
            next
        }
        cat("NaN value in marked Dataset. Repeating ")
    }


    
    
    
    imputation_base_dataset <- replace(data_, data_ == "MCAR",NA)
    imputation_base_dataset <- replace(imputation_base_dataset, imputation_base_dataset == "MNAR",NA)
    imputation_base_dataset <- as.data.frame(imputation_base_dataset)
    imputation_base_dataset <- apply(imputation_base_dataset, 2, function(x) {as.numeric(as.character(x))})
    


   
    #print(imputation_base_dataset)
    modified_imputation_base_dataset <- replace(imputation_base_dataset, imputation_base_dataset == 0,NA)
    modified_imputation_base_dataset <- log(modified_imputation_base_dataset)
    MNAR_Imputation <- impute.QRILC(modified_imputation_base_dataset)[[1]] #QRILC imputation
  

    


    MNAR_Imputation <- exp(MNAR_Imputation)
    modified_imputation_base_dataset <- exp(modified_imputation_base_dataset)
    positions = which(is.na(imputation_base_dataset), arr.ind = TRUE) #missing value locations

    
   
    



 
    

   
    set.seed(42) #reprocubility
   
    MCAR_Imputation <- missForest(modified_imputation_base_dataset)$ximp #RF imputation


    


    #merging datasets

    

    
    for(row in 1:nrow(data_)){

        




        for(col in 1:ncol(data_)){

            if (is.na(data_[row, col])) {
                
                positions = which(is.na(data_), arr.ind = TRUE)
                print(data_before[positions])
                print(data_[positions])
                print(mask[positions])
                stop("Unexpected value in data")
                next
            }
            
            if (data_[row,col] != "MCAR" && data_[row,col] != "MNAR" ){
                next
            } 
            else if(data_[row,col] == "MCAR"){
                imputation_base_dataset[row,col] <- MCAR_Imputation[ row,col]
            }
            else if(data_[row,col] == "MNAR"){
                imputation_base_dataset[row,col] <- MNAR_Imputation[ row,col]
            }
            else {
                print(data_[row,col])
                stop("unecpected value in Imputation Results")
                
            }   

        }

        
    }

    imputated_values <- imputation_base_dataset[positions]

   



    complete_missing_data <- replace(result$complete_missing_data, result$complete_missing_data == "MCAR",NA)
    complete_missing_data <- replace(complete_missing_data, complete_missing_data == "MNAR",NA)
    return(list(imputation_result = imputation_base_dataset,complete_data = result$complete_data,complete_missing_data = complete_missing_data))
    
} 

imputate_Bio_Groups = function(dataset_list){

   # """
    #-takes a list of the batch seperated by biological group
    #-calls imputate function on each group seperatly
    #-merges groupdataframes to one complete dataframe
    #"""

  

    imputated_Dataset_list <- list()

    all_metabolites <- vector()
    skipped_groups <- c()

    for (i in seq_along(dataset_list)){

      skip <- FALSE

      dataset <- dataset_list[[i]]
      data <- unname(dataset)

      if (any(is.na(dataset))){
        tryCatch({
          imputated_data <- imputate(data)
          },error = function(e){
            skipped_groups <<- c(skipped_groups, i) 
            skip <<- TRUE
            })
        
        if (skip){
          next
        }



        result <- imputated_data$imputation_result
      }else {
        result <- data
      }

      

      row_names <- rownames(dataset)
      col_names <- colnames(dataset)

      data <- unname(dataset)


      

      
      

      rownames(result) <- row_names
      colnames(result) <- col_names

      imputated_Dataset_list[[i]] <- result

      if (i == 1){
          all_metabolites <- row_names 
      } else {
          all_metabolites <- union(all_metabolites, row_names)
      }

    }

    for (i in seq_along(imputated_Dataset_list)){ #add rows with only 0 values back in
        dataset <- imputated_Dataset_list[[i]]

        row_names <- rownames(dataset)

        for (metabolite in all_metabolites){

            if (!(metabolite %in% row_names)){

               new_row <- rep(0, ncol(dataset)) 
               new_row_named <- setNames(new_row, colnames(dataset))
               dataset <- rbind(dataset, new_row_named)
               rownames(dataset)[nrow(dataset)] <- metabolite

            }

        }
        
        imputated_Dataset_list[[i]] <- dataset
    }
    

    dfs <- lapply(imputated_Dataset_list, function(m){
        df <- as.data.frame(m)
        df$Row <- rownames(m) 
        df})

    merged_df <- Reduce(function(x, y) full_join(x, y, by="Row"), dfs)
    rownames(merged_df) <- merged_df$Row
    merged_df$Row <- NULL   
    merged_matrix <- as.matrix(merged_df)

    
    return(list(df = merged_matrix,groups = skipped_groups))
}






df_in <- read.csv(args$input_data,header = TRUE, sep = ",",row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

metadata_rows <- c("Group", "TECH_replica","qc")

bio_groups <- unique(as.character(df_in["Group",]))
dfs_by_group <- list()
group_dic <- list()

for (group in bio_groups){

    df_subset_columns <- which(df_in["Group", ] == group)
    df_subset <- df_in[!(rownames(df_in) %in% metadata_rows),df_subset_columns,drop = FALSE]
    df_subset_metadata <- df_in[(rownames(df_in) %in% metadata_rows),df_subset_columns,drop = FALSE]
    df_subset[df_subset == "-"] <- NA
    df_subset[] <- lapply(df_subset, function(x) as.numeric(x))
    group_dic[[group]] <- as.matrix(df_subset_metadata)
    dfs_by_group[[paste0("group_", group)]] <- df_subset
    dfs_by_group <- lapply(dfs_by_group, function(df) {as.matrix(df)})


}


imputation_result = imputate_Bio_Groups(dfs_by_group)
group_dic_filtered <- group_dic[!names(group_dic) %in% imputation_result$groups]
metadata_combined <- do.call(cbind, group_dic_filtered)
df_out = rbind(imputation_result$df,metadata_combined)



write.csv(df_out, file = args$output, row.names = TRUE)













#metadata_rows <- c("Group", "TECH_replica","qc")
#df_metadata <- df_in[rownames(df_in) %in% metadata_rows,,drop = FALSE]
#bio_groups <- unique(as.character(df_in["Group",]))
#dfs_by_group <- list()

#for (group in bio_groups){
#
#    df_subset_columns <- which(df_in["Group", ] == group)
#    df_subset <- df_in[!(rownames(df_in) %in% metadata_rows),df_subset_columns,drop = FALSE]
#    df_subset[df_subset == "-"] <- NA
#    df_subset[] <- lapply(df_subset, function(x) as.numeric(x))
#    dfs_by_group[[paste0("group_", group)]] <- df_subset
#    dfs_by_group <- lapply(dfs_by_group, function(df) {as.matrix(df)})
#
#
#}
#
#df_imputated = imputate_Bio_Groups(dfs_by_group)
#df_out = rbind(df_imputated,df_metadata)
#











