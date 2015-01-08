
bydf = function(df, groupColName, valColName, fun, newColName="temp_by_col_name"){
  byds = by(as.numeric(df[, valColName]), df[, groupColName], fun)
  #byMyDf = data.frame(names = names(byds), vals_9999 = unclass(byds))
  byMyDf = data.frame(names = names(byds), vals_9999 = as.vector(unclass(byds)), stringsAsFactors=FALSE)
  tempDf = renameColDf(byMyDf, "vals_9999", newColName)
  tempDf2 = renameColDf(tempDf, "names", groupColName)
  rownames(tempDf2) = NULL
  return(tempDf2)
}


bydfa = function(df, groupColName, valColName, fun, newColName="temp_by_col_name"){
  byds = by(df[, valColName], df[, groupColName], fun)
  #byMyDf = data.frame(names = names(byds), vals_9999 = unclass(byds))
  byMyDf = data.frame(names = names(byds), vals_9999 = as.vector(unclass(byds)), stringsAsFactors=FALSE)
  tempDf = merge(df, byMyDf, all.x=TRUE, by.x=groupColName, by.y="names")
  tempDf2 = renameColDf(tempDf, "vals_9999", newColName)
  return(tempDf2)
}



delColDf = function(df, colName){
  if(colName %in% names(df)) return( df[ , -which(names(df) %in% c(colName))] )
  else return(df)
}


delete_undefined_rows = function(df, colNames, quiet=FALSE){
  screen2=df
  ncr = nrow(df)
  for(i in 1:length(colNames)){
    if(!is.null(colNames[i])){
      screen2 = subset(screen2, !is.na(screen2[, colNames[i]]) & !is.null(screen2[, colNames[i]]) & screen2[, colNames[i]]!="")
      if(nrow(screen2)<ncr){
        ncr=nrow(screen2)
        if(!quiet) cat(paste('WARNING: Column ', colNames[i], ' contains undefined values. We have deleted the corresponding rows.', sep=""))
      }
    }
  }
  return(screen2)
}


intersectAll = function(...){
  vectors <- list(...)
  intersection = vectors[[1]]
  for(i in 1:length(vectors)){
    intersection = intersect(intersection, vectors[[i]])
  }
  return(intersection)
}


randomizeInner = function (df, baseColStr, sortColStr, reverse = FALSE) 
{
  currentVal = NULL
  baseVals = df[, baseColStr]
  sortVals = df[, sortColStr]
  maxSortVal = 0
  vectToSort = c()
  sortedV = c()
  for (i in 1:nrow(df)) {
    if (is.null(currentVal) || baseVals[i] != currentVal) {
      sortedVtemp = NULL
      if (!is.null(currentVal)) {
        sortedVtemp = sample(1:length(vectToSort))        
        sortedV = c(sortedV, (maxSortVal + sortedVtemp))
        maxSortVal = max((maxSortVal + sortedVtemp))
      }
      vectToSort = c()
      currentVal = baseVals[i]
    }
    vectToSort = c(vectToSort, sortVals[i])
  }
  sortedVtemp = NULL
  sortedVtemp = sample(1:length(vectToSort)) 
  sortedV = c(sortedV, (maxSortVal + sortedVtemp))
  maxSortVal = max((maxSortVal + sortedVtemp))
  return(df[sortedV, ])
}


renameColDf = function(df, colOldName, colNewName){
  if(! (colOldName %in% names(df)) ) print(paste("ERROR: We cannot find ", colOldName, " in the data frame.", sep=""))
  names(df)[ which( names( df ) == colOldName ) ] <- colNewName
  return(df)
}


replace_non_null_elements = function(inputVect, replacementVect){
  for(i in 1:length(inputVect)){
    if(!is.na(replacementVect[i]) && !is.null(replacementVect[i])) inputVect[i] = replacementVect[i]
  }
  return(inputVect)
}



sortInner = function(df, baseColStr, sortColStr, reverse = FALSE){
  currentVal = NULL
  baseVals = df[,baseColStr]
  sortVals = df[, sortColStr]
  maxSortVal = 0
  vectToSort = c()
  sortedV = c()
  for(i in 1:nrow(df)){
    if(is.null(currentVal) || baseVals[i]!=currentVal){
      sortedVtemp = NULL
      if(!is.null(currentVal)) {
        sortedVtemp = order(vectToSort)
        if(reverse) sortedVtemp = order(-vectToSort)
        sortedV = c(sortedV, (maxSortVal + sortedVtemp) )
        maxSortVal = max((maxSortVal + sortedVtemp))
      }
      vectToSort = c()
      currentVal = baseVals[i]
    }
    vectToSort = c(vectToSort, sortVals[i])
  }
  sortedVtemp = NULL
  sortedVtemp = order(vectToSort)
  if(reverse) sortedVtemp = order(-vectToSort)
  sortedV = c(sortedV, (maxSortVal + sortedVtemp) )
  maxSortVal = max((maxSortVal + sortedVtemp))
  return(df[sortedV,])
}


split_df = function(df, strIdCol, linesToGet){
  currentId=df[,strIdCol][1]
  counter=0
  resultV=c()
  for(i in 1:nrow(df)){
    if(df[,strIdCol][i]!=currentId){
      counter=0
      currentId=df[,strIdCol][i]
    }  
    counter=counter+1
    if(counter %in% linesToGet){
      resultV=c(resultV,i)
    }
  }
  return(df[resultV,])
}


randomSortOnVal = function(screen, strColVal){
  screen = arrange(screen, screen[,strColVal])
  randV = sample(1:nrow(screen))
  currentGene=NULL
  j=0
  newGeneID = rep(NA, nrow(screen))
  colVal = screen[, strColVal]
  for(i in 1:nrow(screen)){
    if(is.null(currentGene) || currentGene!=colVal[i]){
      currentGene=colVal[i]
      j=j+1
    }
    newGeneID[i] = randV[j]
  }
  screen = data.frame(screen, newGeneID=newGeneID)
  screen = arrange(screen, newGeneID)
  screen = delColDf(screen,  "newGeneID")
  return(screen)
}


transcribe_seqs = function(df, seqColName="siRNA_seq", toDNA=FALSE, progress_bar=FALSE){
    if( exists("progress_bar") && progress_bar){
        pb <- txtProgressBar(min = 0, max = nrow(df)/1000, style = 3, width=100)
    }
    seqV=as.character(df[, seqColName])
    seqVtranscribed = rep(NA, length(seqV))
    for(i in 1:length(seqV)){
        if(exists("progress_bar") && progress_bar && i%%1000==0 ){
            setTxtProgressBar(pb, i/1000)
        }
        if(!is.na(seqV[i]) & nchar(seqV[i])>0 ){
            if(grepl( "U", seqV[i]) | grepl("u", seqV[i])) {
                seqTemp = DNAString(RNAString(seqV[i]))
            }else{
                seqTemp = DNAString(seqV[i])
            }
            seqVtranscribed[i] =  as.character(RNAString(reverseComplement(DNAString(seqTemp))))
            if(toDNA) {
                seqVtranscribed[i] = as.character(DNAString(RNAString(seqVtranscribed[i])))
            }
        }
    }
    df[seqColName] = seqVtranscribed 
    return(df)
}


