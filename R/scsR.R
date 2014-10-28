add_rank_col = function(screen, reverse=FALSE,  scoreColName="score", geneColName="GeneID"){
  genesMedian = sqldf(paste('select ',geneColName,', median(',scoreColName,') as median from screen group by ',geneColName, sep=""), stringsAsFactors=FALSE)
  screen5 = merge(screen, genesMedian, by.x=geneColName, by.y=geneColName, x.all=TRUE, na.rm=TRUE)
  genesAvg = sqldf(paste('select ',geneColName,', avg(',scoreColName,') as average from screen group by ',geneColName, sep=""), stringsAsFactors=FALSE)
  screen5 = merge(screen5, genesAvg, by.x=geneColName, by.y=geneColName, x.all=TRUE, na.rm=TRUE)
  genesMax = sqldf(paste('select ',geneColName,', max(',scoreColName,') as max from screen group by ',geneColName, sep=""), stringsAsFactors=FALSE)
  screen5 = merge(screen5, genesMax, by.x=geneColName, by.y=geneColName, x.all=TRUE, na.rm=TRUE)
  genesMin = sqldf(paste('select ',geneColName,', min(',scoreColName,') as min from screen group by ',geneColName, sep=""), stringsAsFactors=FALSE)
  screen5 = merge(screen5, genesMin, by.x=geneColName, by.y=geneColName, x.all=TRUE, na.rm=TRUE)
  rsa_df = launch_RSA(screen5, LB=-10, UB=10, reverse=reverse, strScoreCol=scoreColName, strGeneCol=geneColName)
  return(rsa_df)
}


add_seed = function(df, seqColName="siRNA_seq", seedLength=7, startPosition=2){
  v=df[, seqColName]
  seedv = rep(NA, length(v))
  for(i in 1:length(v)){
    seed = substring(v[i],startPosition,(seedLength+startPosition-1))
    seed = toupper(seed)
    seed = gsub("T", "U", seed)
    seedv[i] = sub("T", "U", seed)
  }
  df[, paste("seed", seedLength, sep="")]=seedv
  return(df)
}


check_consistency = function(screen, scoreColName = "score", geneColName = "GeneID",
                             seqColName="siRNA_seq"){
  countGenes <- NULL #initialize to not get useless complains from R CMD CHECK 
  x<-NULL #initialize to not get useless complains from R CMD CHECK
  if(length(grep("[a-z]", screen[,seqColName]))){
    cat('WARNING: Your sequences contain lower case letters. For consistency we will convert them to upper case letters.')
    screen[,seqColName] = toupper(screen[,seqColName])
  }
  if(length(grep("[A,C,G,T,U]", screen[,seqColName], invert=TRUE))){
    cat('ERROR: Unsupported characters have been found in your sequences.\n')
    cat('        Your sequences must be RNA (or DNA) sequences. Please check again your sequences and try again.')
    stop()
  }
  if(length(unique(screen[,seqColName])) != length(screen[,seqColName]) ) 
    cat("WARNING: Your screen seems to contain replicates (i.e. some sequences are present more than once) ! 
        You should take the median of the scores of those replicates.
        You can use our median_replicates function in order to do that.\n")
  if(length(grep("T", screen[,seqColName])) >0) 
    cat("WARNING: Your sequences contain the letter T. 
        Are you sure to have provided the correct antisense sequences ? 
        Please remember that you must provide the antisense sequence (NOT the target sequence !) \n")  
  af=subset(count(substring(screen[,seqColName],1,1)), x=="A")$freq
  cf=subset(count(substring(screen[,seqColName],1,1)), x=="C")$freq
  gf=subset(count(substring(screen[,seqColName],1,1)), x=="G")$freq
  tf=subset(count(substring(screen[,seqColName],1,1)), x=="T")$freq
  uf=subset(count(substring(screen[,seqColName],1,1)), x=="U")$freq
  if(length(tf)==0){tf=0}
  if(length(uf)==0){uf=0}
  if(length(af)==0){af=0}
  if((uf+tf+af) < 0.8*nrow(screen)){
    cat("WARNING: Your sequences contain in the first nucleotide a low percentage of T, U and A. 
        Are you sure to have provided the correct antisense sequences ? 
        Please remember that you must provide the antisense sequence (NOT the target sequence !) \n")
  }
  screen=delete_undefined_rows(screen, colNames=c(seqColName, scoreColName, geneColName))
  screen2=sqldf(paste('select count(distinct ',geneColName,') as countGenes, ',seqColName,' from screen group by ', seqColName, sep=""))
  if(nrow(subset(screen2, countGenes!=1))>0) print("WARNING: Some of your siRNA sequence refers to more than one Gene Identifier ! Please check your input")
  return(screen)
}


compare_sorted_geneSets = function (genesetA1, genesetA2, genesetB, background, limA=NULL, limB=NULL) 
{
  genesetA1 = unique(genesetA1[genesetA1 %in% background])
  genesetA2 = unique(genesetA2[genesetA2 %in% background])
  genesetB = unique(genesetB[genesetB %in% background])
  
  if(is.null(limA)){
    minTemp = min(length(genesetA1), length(genesetA2))
    genesetA1 = genesetA1[1:minTemp]
    genesetA2 = genesetA2[1:minTemp]
  }else{
    if(length(genesetA1) <limA){cat("WARNING: The length of genesetA1 is less than limA")}
    if(length(genesetA1) <limA){cat("WARNING: The length of genesetA2 is less than limA")}
    genesetA1 = genesetA1[1:limA]
    genesetA2 = genesetA2[1:limA]
  }
  if(!is.null(limB)){
    if(length(genesetB) <limB){cat("WARNING: The length of genesetB is less than limB")}
    genesetB = genesetB[1:limB]
  }
  enrichment_geneSet(genesetA1, genesetB, background)
  enrichment_geneSet(genesetA2, genesetB, background)
}


create_sd_matrix = function(screen, seedColName="seed7", scoreColName="score"){
  nintervals = 8
  quntils = 21
  sd_matr = matrix(NA, nintervals+2, quntils)
  sa = subset(seeds_analysis(screen, seedColName=seedColName, scoreColName=scoreColName), !is.na(sd) )
  sd_matr[1,] = quantile(subset(sa, sa[,scoreColName]<=-1.5)$sd, probs = seq(0, 1, 0.05))
  sd_matr[10,] = quantile(subset(sa, sa[,scoreColName]>1.5)$sd, probs = seq(0, 1, 0.05))
  middleSa = subset(sa, sa[,scoreColName]>-1.5 && sa[,scoreColName]<=1.5)
  for(i in 1:nintervals){
    tempSa = subset(sa, sa[,scoreColName]>(-1.5 + 3/nintervals*(i-1)) & sa[,scoreColName]<=(-1.5 + 3/nintervals*i))
    sd_matr[i+1,] = quantile(tempSa$sd, probs = seq(0, 1, 0.05))
  }
  return(sd_matr)
}



enrichment_geneSet = function(genesetA, genesetB, background=NULL, quiet=FALSE){
  genesetA=unique(genesetA)
  genesetB=unique(genesetB)
  backgroundSize = 18000
  if(!is.null(background)){
    background=unique(background)
    backgroundSize = length(background)
    genesetA=intersect(genesetA, background)
    genesetB=intersect(genesetB, background)
  } 
  intersectionGeneSets = intersect(genesetA, genesetB)
  k=as.numeric(length(genesetA))
  n=as.numeric(backgroundSize) - as.numeric(length(genesetB))
  m=as.numeric(length(genesetB))
  x=as.numeric(length(intersectionGeneSets))
  if(!quiet) {
    cat(paste("sizeA: ", length(genesetA), "       sizeB: ", length(genesetB),
              "       sizeCommon: ", length(intersectionGeneSets), "\n",sep=""))
    cat(paste( "p-value (Hypergeometric test): ", phyper(x, m, n, k, FALSE), "\n",sep=""))
  }
  return(phyper(x, m, n, k, FALSE))
}


# enrichment_ppi_vector <- function(genesVectors, limit=400, STRINGversion="9_05", species_ncbi_taxonomy_id=9606){
#   
#   string_db <- STRINGdb$new( version=STRINGversion, species=species_ncbi_taxonomy_id, score_threshold=0, input_directory="" )
#   enrichList = c()
#   
#   genes=genesVectors[[1]]
#   for(i in 2:length(genesVectors)){ genes=intersect(genes, genesVectors[[i]])  }
#   for(i in 1:length(genesVectors)){genesVectors[[i]] = intersect(genesVectors[[i]], genes)}
#   cat("INFO: we are benchmarking ",length(genes), "genes")
#   
#   for(i in 1:length(genesVectors)){
#     stringGenes = string_db$mp(genesVectors[[i]][0:limit])
#     enrichList[i] = string_db$get_ppi_enrichment(stringGenes)$enrichment
#   }
#   return(enrichList)
# }


get_sd_quant = function(sdval, score, sd_matrix){
  perc = NA
  if(score<=-1.5) perc = table(findInterval(sd_matrix[1,], sdval))[1]
  else if(score>1.5) perc = table(findInterval(sd_matrix[10,], sdval))[1]
  else{
    perc = table(findInterval(sd_matrix[floor(((score+1.5)/(3/8))+2),], sdval))[1]
  }
  if(names(perc)=="1") perc = 1
  perc = min(perc, 20)
  return(perc)
}


get_seed_oligos_df=function(screen, seedColName="seed7",  scoreColName="score", geneColName="GeneID", gene_interval = c(1,100), 
                            min_oligos_x_gene=4, min_oligos_x_statistics=4, random=FALSE, kolmogorovSampleSize=5000, progress_bar=FALSE ){
  
  seeds7 = subset(seeds_analysis(screen, seedColName=seedColName, scoreColName=scoreColName), select=c(seedColName, "pvalue"))
  
  if(random) screen = randomSortOnVal(screen, geneColName)
  
  # extract the interesting oligos from the whole screen
  tempGenes3 = unique(screen[,geneColName])
  goodGenes = tempGenes3
  if(!is.null(min_oligos_x_gene)){
    tempGenes1 = as.numeric(sqldf(paste('select ', geneColName,' from screen group by GeneID having count(GeneID) >=',min_oligos_x_gene, paste=""))$GeneID)
    tempGenes2 = unique(screen[,geneColName])
    tempGenes3 = tempGenes2[tempGenes2 %in% tempGenes1]
  }
  if(!is.null(gene_interval)) goodGenes = tempGenes3[gene_interval[1]:gene_interval[2]]
  goodGenesDf = subset(screen, screen[,geneColName] %in% goodGenes)
  screent1 = subset(screen, screen[,seedColName] %in% goodGenesDf[,seedColName])
  
  screent2 = subset(screent1, select=c("GeneID", seedColName, scoreColName))
  hit_th_val=unname(quantile(screen[, scoreColName], c(0.1)))
  hit_th_val_enh=unname(quantile(screen[, scoreColName], c(0.9)))
  if(!is.null(kolmogorovSampleSize)) scoreSample = sample(screen[, scoreColName])[1:min(nrow(screen), kolmogorovSampleSize)]
  graphMatr = matrix(NA, nrow(goodGenesDf), 5)
  
  if(exists("progress_bar") && progress_bar) pb <- txtProgressBar(min = 0, max = nrow(goodGenesDf)/10, style = 3, width=100)
  currentGene = NULL
  gn=0
  for(i in 1:nrow(goodGenesDf)){
    if((exists("progress_bar") && progress_bar ) && i%%10==0) setTxtProgressBar(pb, i/10)
    scores = subset(screent2, screent2[,geneColName] != goodGenesDf[, geneColName][i] & screent2[,seedColName]==goodGenesDf[, seedColName][i])[, scoreColName]
    
    if(length(scores)>=min_oligos_x_statistics){
      countHits = length(scores[scores<=hit_th_val])
      countHitsEnh = length(scores[scores>=hit_th_val_enh])
      coeff=1
      if(mean(scores) >0){
        coeff = (-1)
        graphMatr[i,2] = coeff*max(log(phyper((countHitsEnh-1), (nrow(screen)*0.1), nrow(screen)-(nrow(screen)*0.1), length(scores), lower.tail=FALSE), base=2),-10)
      }else{
        graphMatr[i,2] = max(log(phyper((countHits-1), (nrow(screen)*0.1), nrow(screen)-(nrow(screen)*0.1), length(scores), lower.tail=FALSE), base=2), -10)
      }
      if(!is.null(kolmogorovSampleSize)) graphMatr[i,1] = coeff*max(log(ks.test(scores, scoreSample)$p.value, base=2), -10)
    }else{
      if(!is.null(kolmogorovSampleSize)) graphMatr[i,1] = 0
      graphMatr[i,2] = 0
    }
    graphMatr[i,3] = length(scores)
    if(length(scores)>0) graphMatr[i,4] = mean(scores) 
    
    if(is.null(currentGene) || goodGenesDf[,geneColName][i]!=currentGene){
      currentGene = goodGenesDf[,geneColName][i]
      gn=gn+1
    }
    graphMatr[i,5] = gn
  }
  
  resultDf = data.frame(goodGenesDf, seed_log_pval_ks = graphMatr[,1], seed_log_pval = graphMatr[,2], 
                        seed_oligos_count = graphMatr[,3], seed_oligos_mean = graphMatr[,4], geneRank = graphMatr[,5] )
  
  if(!is.null(kolmogorovSampleSize)) 
    resultDf2 = sortInner(resultDf, geneColName, "seed_log_pval_ks", reverse = FALSE)
  else resultDf2 = sortInner(resultDf, geneColName, "seed_log_pval", reverse = FALSE)
  
  return(resultDf2)
  
}


launch_RSA = function(df, LB=-100, UB=100, reverse=FALSE, strScoreCol="", strGeneCol="Gene_ID", keepAllRSAReturnFields=FALSE){
  rsaOptions = list(LB=LB,UB=UB,reverse=reverse);
  df_filtered = df[!is.na(df[,strGeneCol]) & df[,strGeneCol] != "" & !is.na(df[,strScoreCol]),]
  df_rsa = OPIrsa(df_filtered[,strGeneCol], df_filtered[,strScoreCol], rsaOptions, df_filtered)  
  
  df_rsa = renameColDf(df_rsa, "OPI_Rank", "rank_rsa")
  df_rsa = renameColDf(df_rsa, "LogP", "log_pval_rsa")
  df_rsa = renameColDf(df_rsa, "X.totalWell", "num_wells")
  
  if(! keepAllRSAReturnFields){
    df_rsa = delColDf(df_rsa, "EXP_Rank")
    df_rsa = delColDf(df_rsa, "Cutoff_Rank")
    df_rsa = delColDf(df_rsa, "rank")
    df_rsa = delColDf(df_rsa, "X.hitWell")
    df_rsa = delColDf(df_rsa, "OPI_Hit")
    df_rsa = delColDf(df_rsa, "Score")
  }
  
  return(df_rsa)
}


median_replicates = function(screen, seedColName = "seed7", scoreColName = "score",
                             geneColName = "GeneID", seqColName="siRNA_seq", spAvgColName = NULL
){
  spAvgCmd = ""
  if(!is.null(spAvgColName)) spAvgCmd=paste(', median(',spAvgColName,') as ', spAvgColName,  sep="")
  aggregateDf5 = sqldf(paste('select median(', geneColName,') as ',geneColName,', median(',scoreColName,') as ',scoreColName,' ,
                             ',seqColName,spAvgCmd,' from screen group by ',seqColName,' order by ',scoreColName, sep=""))
  return(aggregateDf5)
}


OPIrsa = function(Groups,Scores,opts,Data=NULL)
{
  t = data.frame(Gene_ID = Groups, Score = Scores)
  Data = data.frame(Data, Score = Scores)
  Sorted_Order = order(t$Score,decreasing=opts$reverse);
  Data = Data[Sorted_Order,]
  t = t[Sorted_Order,]
  t = do.call("rbind", tapply(seq(nrow(t)), list(t$Gene_ID), function(i,dataset, optsb) 
  {
    if(optsb$reverse)
    {
      i_max = sum(dataset$Score[i]>=optsb$LB)
      i_min = max(1,sum(dataset$Score[i]>=optsb$UB))
    }else
    {
      i_max = sum(dataset$Score[i]<=optsb$UB)
      i_min = max(1,sum(dataset$Score[i]<=optsb$LB))
    }
    r = OPIrsaScore(i,nrow(dataset),i_min,i_max)
    return ( cbind(
      LogP = r["logp"]
      ,OPI_Hit=as.numeric(seq(length(i))<=r["cutoff"])
      ,"#hitWell"=i_max
      ,"#totalWell"=length(i)
      ,rank = i))
  }, dataset = t, opts))
  t_sorted = t[order(t[,"rank"]),]
  t = data.frame(cbind(Data, t_sorted))
  
  # add OPI_Rank
  t = t[order(t$LogP,t$Score*ifelse(opts$reverse,-1,1)),]
  t$OPI_Rank = cumsum(t$OPI_Hit)
  t$OPI_Rank[t$OPI_Hit == 0] = 999999
  
  # add Cutoff_Rank
  t = t[order(t$Score*(ifelse(opts$reverse,-1,1)),t$LogP),]
  
  if(opts$reverse){tmp = t$Score>=opts$LB} else {tmp = t$Score<=opts$UB}
  t$Cutoff_Rank = cumsum(tmp)
  t$Cutoff_Rank[!tmp] = 999999
  
  # add EXP_Rank
  t$EXP_Rank = pmin(t$OPI_Rank,t$Cutoff_Rank)
  t$EXP_Rank = pmin(t$OPI_Rank,t$Cutoff_Rank)
  if(opts$reverse) {
    return(t[order(t$OPI_Rank, -t$Score),])
  } else {
    return(t[order(t$OPI_Rank, t$Score),])
  }
}


OPIrsaScore = function(I_rank, N, i_min=1, i_max=-1)
{
  n_drawn = length(I_rank) # number of black
  if(i_max == -1)
  {
    i_max=n_drawn
  }
  r1 = c(logp=1.0,cutoff=0)
  if( i_max < i_min) return (r1)
  # phyper(x, lower.tail = F), x = x-1, when lower.tail = F
  logp =  apply(cbind(seq(i_min,i_max),I_rank[i_min:i_max]),1,function(x) { phyper(x[1]-1,x[2] ,N-x[2], n_drawn,lower.tail = FALSE,log.p=TRUE)})
  
  logp = logp/log(10)
  logp[logp<(-100)] = -100
  if(all(is.na(logp))) {
    return  (r1)
  }else
    return   ( c(logp=min(logp),cutoff = i_min-1+which.min(logp)))
}


removeSharedOffTargets = function(screenA, screenB, seedColName="seed7",
                                  geneColName="GeneID",
                                  seqColName="siRNA_seq",
                                  removeGenes=FALSE){
  seeds = intersect(screenA[,seedColName], screenB[,seedColName])
  genes = intersect(screenA[,geneColName], screenB[,geneColName])
  screenA2 = subset(screenA, screenA[,seedColName] %in% seeds & screenA[,geneColName] %in% genes)
  screenB2 = subset(screenB, screenB[,seedColName] %in% seeds & screenB[,geneColName] %in% genes)
  gene_seeds=hash()
  for(i in 1:nrow(screenB2)){
    geneSeeds = c()
    if(has.key(as.character(screenB2[i,geneColName]) ,gene_seeds)) geneSeeds = gene_seeds[[as.character(screenB2[i,geneColName])]]
    geneSeeds=c(geneSeeds, as.character(screenB2[i,seedColName]))
    gene_seeds[[as.character(screenB2[i,geneColName])]] = geneSeeds
  }
  
  oligosToRemove = c()
  for(i in 1:nrow(screenA2)){
    if(has.key(as.character(screenA2[i,geneColName]) ,gene_seeds)){
      if(as.character(screenA2[i,seedColName]) %in% gene_seeds[[as.character(screenA2[i,geneColName])]]){
        oligosToRemove=c(oligosToRemove, screenA2[i,seqColName])
      }
    } 
  }
  if(removeGenes){
    genesToRemove = unique(subset(screenA2, screenA2[, seqColName] %in% oligosToRemove)[, geneColName])
    screenA3 = subset(screenA, !(screenA[, geneColName] %in% genesToRemove))
  }else{
    screenA3 = subset(screenA, !(screenA[, seqColName]  %in% oligosToRemove))
  }
  
  return(screenA3)
}


seeds_analysis = function(screen, seedColName="seed7",  scoreColName="score", hit_th_val=NULL, 
                          enhancer_analysis=FALSE, spAvgColName=NULL, 
                          minCount=NULL, ks_enabled=FALSE, miRBase=NULL){
  pvalue=NULL
  miRBase_df=NULL
  if(!is.null(miRBase)){
    #if(class(miRBase)!="RNAStringSet"){
     # cat("ERROR: the miRBase parameter should be an instance of RNAStringSet class of the package biostrings.")
    #  stop()
    #}
    #miRBase_df=data.frame(miRNA=names(miRBase), seq=as.character(miRBase))
    #rownames(miRBase_df)=NULL
    miRBase_df=miRBase
    if(!("miRNA" %in% colnames(miRBase_df) & "seq" %in% colnames(miRBase_df) )){
      cat("ERROR: miRBase should be a data frame containing the columns miRNA and seq.\n")
      cat("         Please check the name of the columns. ")
      stop()
    }
  }
  
  if(is.null(hit_th_val)){
    if(enhancer_analysis) hit_th_val=unname(quantile(screen[,scoreColName], c(0.9)))
    else hit_th_val=unname(quantile(screen[,scoreColName], c(0.1)))
  }
  spAvgCmd =""
  if(!is.null(spAvgColName)) spAvgCmd = paste(',avg(', spAvgColName, ') as ', spAvgColName, sep="")
  seeds = sqldf(paste('select ',seedColName,', avg(',scoreColName,') as score, count(',seedColName,') as count',spAvgCmd,' from screen group by ',seedColName, sep=""))
  screen6 = sqldf( paste('select * from screen where ',scoreColName,'<=', hit_th_val, sep="") )
  if(enhancer_analysis) screen6 = sqldf( paste('select * from screen where ',scoreColName,'>=', hit_th_val, sep="") )
  seedsHits = sqldf(paste('select ',seedColName,', count(',seedColName,') as countHits from screen6 group by ',seedColName, sep=""))
  seeds2 = merge(seeds, seedsHits, all.x=TRUE)
  seeds21=seeds2
  if(ks_enabled==TRUE){
    scoreSample = sample(screen[, scoreColName])[1:min(nrow(screen), 5000)]
    seedsPvaluesKs = bydf(screen, seedColName, scoreColName, function(x){return(ks.test(x, scoreSample)$p.value)}, newColName="pvalue_ks")
    seeds21 = merge(seeds2, seedsPvaluesKs, all.x=TRUE)
  }
  ratiosHitsVsCount = NULL
  for(i in 1:nrow(seeds21)){ ratiosHitsVsCount = c(ratiosHitsVsCount, seeds21$countHits[i]/seeds21$count[i])}
  seeds3 = data.frame(seeds21, ratiosHitsVsCount = ratiosHitsVsCount)
  seeds4 = data.frame(seeds3, pvalue=phyper((seeds3$countHits-1), nrow(screen6), nrow(screen)-nrow(screen6), seeds3$count, lower.tail=FALSE) )
  if(!is.null(minCount)) seeds4=seedsAnalysisDf = sqldf(paste('select * from seeds4 where count >= ', minCount, sep=""))
  sddf = bydf(screen, seedColName, "score", sd, "sd")
  seeds5 = merge(seeds4, sddf, all.x=TRUE, by.x=seedColName, by.y=seedColName)
  seeds6 = bydfa(seeds5, seedColName,  seedColName, function(seq){
    seq = toupper(seq)
    seq = strsplit(seq,"")[[1]]
    counter=0
    for(i in 1:length(seq)) 
      if(seq[i]=="C" || seq[i]=="G") counter=counter+1
    return(counter)
  }, "countGC")
  if(enhancer_analysis) {seeds7 = arrange(seeds6, -score)}else{seeds7 = arrange(seeds6, score)}
  seeds7[, seedColName] = as.character(seeds7[, seedColName])
  
  if(!is.null(miRBase_df)){
    miRBase_df2 = add_seed(miRBase_df, seedLength=7, startPosition=2, seqColName="seq")
    mirhash=hash()
    for(i in 1:nrow(miRBase_df2)){
      mirsV = c()
      if(has.key(as.character(miRBase_df2$seed7[i]), mirhash)) mirsV=mirhash[[as.character(miRBase_df2$seed7[i])]]
      mirsV = c(mirsV, as.character(miRBase_df2$miRNA[i]) )
      mirhash[[as.character(miRBase_df2$seed7[i])]]=mirsV
    }
    miRNAcol=rep(NA, nrow(seeds7))
    for(i in 1:nrow(seeds7)){
      if(has.key(as.character(seeds7[,seedColName][i]),mirhash)){
        miRNAcol[i] = paste(mirhash[[as.character(seeds7[,seedColName][i])]], collapse=", ")
      }
    }
    seeds7=data.frame(seeds7, miRNA=miRNAcol, stringsAsFactors=FALSE)
  }
  
  res1 = arrange(seeds7, pvalue)
  res1$score = round(res1$score, digits=4)
  res1$ratiosHitsVsCount = round(res1$ratiosHitsVsCount, digits=4)
  res1$sd = round(res1$sd, digits=4)
  return(res1)
}


seed_correction = function(screen, seedColName="seed7", scoreColName="score", 
                           geneColName="GeneID", fixed_correction_coeff=0.4, 
                           sd_correction_coeff=0.6, min_siRNAs_x_seed=3, progress_bar=FALSE){
  if( exists("progress_bar") && progress_bar) pb <- txtProgressBar(min = 0, max = nrow(screen)/1000, style = 3, width=100)
  previous_order_temp<-NULL
  screen = data.frame(screen, previous_order_temp = seq(1:nrow(screen)))
  sd_matrix = create_sd_matrix(screen, seedColName="seed7", scoreColName="score")
  screen = arrange(screen, screen[,seedColName])
  medianScore = median(screen[,scoreColName])
  zscoresV = c()
  scoresV = c()
  currentSeed=NULL
  for(i in 1:nrow(screen)){
    if(( exists("progress_bar") && progress_bar ) && i%%1000==0) setTxtProgressBar(pb, i/1000)
    if(is.null(currentSeed)) currentSeed = screen[,seedColName][i]
    if(currentSeed != screen[,seedColName][i]){
      zscores2V = c()
      for(j in 1:length(scoresV)){
        scores2V = scoresV[-j]
        if(length(scores2V) >= min_siRNAs_x_seed ){
          deltaS = medianScore - median(scores2V)
          zscores2V[j] = scoresV[j] +  ( fixed_correction_coeff + (sd_correction_coeff/20)*(20-get_sd_quant(sd(scores2V), mean(scores2V), sd_matrix)) )*deltaS
        }else{
          zscores2V[j] = NA
        }
      }
      zscoresV = c(zscoresV, zscores2V)
      currentSeed = screen[,seedColName][i]
      scoresV=c()
    } 
    scoresV = c(scoresV, screen[,scoreColName][i])
  }
  zscores2V = c()
  for(j in 1:length(scoresV)){
    scores2V = scoresV[-j]
    if(length(scores2V) >= min_siRNAs_x_seed ){
      deltaS = medianScore - median(scores2V)
      zscores2V[j] = scoresV[j] +  ( fixed_correction_coeff + (sd_correction_coeff/20)*(20-get_sd_quant(sd(scores2V), mean(scores2V), sd_matrix)) )*deltaS
    }else zscores2V[j] = NA
  }
  zscoresV = c(zscoresV, zscores2V)
  screen$score = replace_non_null_elements(screen$score, zscoresV)
  screen = arrange(screen, previous_order_temp)
  screen = delColDf(screen, "previous_order_temp")
  return(screen)
}



seed_correction_pooled = function(screen, seedColName="seed7", scoreColName="score", 
                                  geneColName="GeneID", fixed_correction_coeff=0.4, 
                                  sd_correction_coeff=0.6, min_siRNAs_x_seed=4, poolSize=4, enhancer_analysis=NULL, 
                                  use_all_seeds=TRUE, progress_bar=FALSE){
  
  if( exists("progress_bar") && progress_bar) pb <- txtProgressBar(min = 0, max = nrow(screen)/1000, style = 3, width=100)
  previous_order_temp<-NULL
  screen = data.frame(screen, previous_order_temp = seq(1:nrow(screen)))
  sd_matrix = create_sd_matrix(screen, seedColName="seed7", scoreColName="score")
  screen = arrange(screen, screen[,seedColName])
  medianScore = median(screen[,scoreColName])
  meanScore = mean(screen[,scoreColName])
  zscoresV = c()
  zdeltasV = c()
  scoresV = c()
  currentSeed=NULL
  for(i in 1:nrow(screen)){
    if(( exists("progress_bar") && progress_bar ) && i%%1000==0) setTxtProgressBar(pb, i/1000)
    if(is.null(currentSeed)) currentSeed = screen[,seedColName][i]
    if(currentSeed != screen[,seedColName][i]){
      zdeltas2V = c()
      zscores2V = c()
      for(j in 1:length(scoresV)){
        scores2V = scoresV[-j]
        if(length(scores2V) >= min_siRNAs_x_seed ){
          deltaS = medianScore - median(scores2V)
          zscores2V[j] = scoresV[j] +  ( fixed_correction_coeff + (sd_correction_coeff/20)*(20-get_sd_quant(sd(scores2V), mean(scores2V), sd_matrix)) )*deltaS
          zdeltas2V[j] =  ( fixed_correction_coeff + (sd_correction_coeff/20)*(20-get_sd_quant(sd(scores2V), mean(scores2V), sd_matrix)) )*deltaS
        }else{
          zscores2V[j] = NA
          zdeltas2V[j] = NA
        }
      }
      zscoresV = c(zscoresV, zscores2V)
      zdeltasV = c(zdeltasV, zdeltas2V)
      currentSeed = screen[,seedColName][i]
      scoresV=c()
    } 
    scoresV = c(scoresV, screen[,scoreColName][i])
  }
  zscores2V = c()
  zdeltas2V = c()
  for(j in 1:length(scoresV)){
    scores2V = scoresV[-j]
    if(length(scores2V) >= min_siRNAs_x_seed ){
      deltaS = medianScore - median(scores2V)
      zscores2V[j] = scoresV[j] +  ( fixed_correction_coeff + (sd_correction_coeff/20)*(20-get_sd_quant(sd(scores2V), mean(scores2V), sd_matrix)) )*deltaS
      zdeltas2V[j] = ( fixed_correction_coeff + (sd_correction_coeff/20)*(20-get_sd_quant(sd(scores2V), mean(scores2V), sd_matrix)) )*deltaS
    }else{ 
      zscores2V[j] = NA
      zdeltas2V[j] = NA
    }
  }
  zscoresV = c(zscoresV, zscores2V)
  zdeltasV = c(zdeltasV, zdeltas2V)
  screen = data.frame(screen, deltas=zdeltasV)
  screen = data.frame(screen, scores_old=screen$score)
  screen = bydfa(screen, "GeneID", "deltas", function(x){return(.pooled_weightened_mean(x, medianScore, poolSize=poolSize, enhancer_analysis=enhancer_analysis, use_all_seeds=use_all_seeds))}, "deltaSum")
  screen$score = screen$score + screen$deltaSum
  
  screen = arrange(screen, previous_order_temp)
  screen = delColDf(screen, "previous_order_temp")
  screen = delColDf(screen, "deltaSum")
  screen = delColDf(screen, "scores_old")
  
  return(screen)
}


.pooled_weightened_mean <- function(vect, baseVal, enhancer_analysis, poolSize=4, use_all_seeds=TRUE){
  vect = vect[!is.na(vect)] 
  if(length(vect)==0) return(0)
  maxVal = max(vect, na.rm=TRUE)
  minVal = min(vect, na.rm=TRUE)
  dist = 0
  delta = 0
  if(!use_all_seeds){
    if(is.null(enhancer_analysis)){
      cat("\nERROR: If you want to apply a minmax analysis than you need to set the enhancer_analysis variable. \nLook at the documentation\n") 
      stop()
    }
    if(enhancer_analysis) {  if(minVal < baseVal){delta = minVal}}else{if(maxVal > baseVal){delta = maxVal}}
  }else{
    for(i in 1:length(vect)){
      x = vect[i] - baseVal
      dist = dist + abs(x)
    }
    if(length(vect) < poolSize){
      dist = (dist * poolSize) / length(vect)
    }
    for(i in 1:length(vect)){
      x = vect[i] - baseVal
      safetyFactor=1
      #if(enhancer_analysis && x>0){safetyFactor=0.5 }
      #if(!enhancer_analysis && x<0){safetyFactor=0.5 }
      delta = delta + safetyFactor * (abs(x)*x)/(dist)
    }  
  }
  return(delta)
}



seed_removal = function(screen,seedColName="seed7", scoreColName="score", geneColName="GeneID",
                             min_siRNAs_x_seed=4, remove_unrepresented_seeds=TRUE, lower_bound_threshold = -0.5,
                             higher_bound_threshold = 0.5, min_oligos_x_gene_threshold = 2, 
                        useMedian=FALSE, removeGenes=FALSE, include_current_gene=FALSE, progress_bar=FALSE){
 
  if(include_current_gene){
    if(removeGenes){
      cat("ERROR: We currently don't support both 'removeGenes' and 'include_current_gene' activated.\n ")
      cat("Please change the input parameters.")
      stop()
    }
    #TODO: add support for the median
    seeds = seeds_analysis(screen, seedColName="seed7", scoreColName="score")
    seedsToRemove = subset(seeds, count>=min_siRNAs_x_seed & (score<lower_bound_threshold | score>higher_bound_threshold))
    screen2 = subset(screen, !(screen[,seedColName] %in% seedsToRemove[,seedColName]))
    if(remove_unrepresented_seeds) screen2 = subset(screen2, screen2[,seedColName] %in% subset(seeds, count>=min_siRNAs_x_seed)[,seedColName] )
    if(!is.null(min_oligos_x_gene_threshold)){
      genesToRemove = sqldf(paste('select ', geneColName, ', count(', geneColName, ') from screen2 group by ', geneColName, ' having count(', geneColName,') <', min_oligos_x_gene_threshold, sep=""))[, geneColName]
      screen2 = subset(screen2, !(screen2[, geneColName] %in% genesToRemove) )
    }
    return(screen2)
    
  }else{
    previous_order_temp<-NULL
    genesToRemove=rep(NA, nrow(screen))
    seedsToRemove=rep(NA, nrow(screen))
    linesToRemove=rep(NA, nrow(screen))
    if( exists("progress_bar") && progress_bar) pb <- txtProgressBar(min = 0, max = nrow(screen)/1000, style = 3, width=100)
    screen = data.frame(screen, previous_order_temp = seq(1:nrow(screen)))
    screen = arrange(screen, screen[,seedColName])
    medianScore = median(screen[,scoreColName])
    meanScore = mean(screen[,scoreColName])
    scoresV = c()
    genesV=c()
    currentSeed=NULL
    for(i in 1:nrow(screen)){
      if((exists("progress_bar") && progress_bar ) && i%%1000==0) setTxtProgressBar(pb, i/1000)
      if(is.null(currentSeed)) currentSeed = screen[,seedColName][i]
      if(currentSeed != screen[,seedColName][i]){
        for(j in 1:length(scoresV)){
          scores2V = scoresV[-j]
          if(length(scores2V) >= min_siRNAs_x_seed ){
            if(
              (!useMedian && (mean(scores2V) < lower_bound_threshold || mean(scores2V) > higher_bound_threshold) ) ||
                (useMedian && (median(scores2V) < lower_bound_threshold || median(scores2V) > higher_bound_threshold))
            )
            {
              genesToRemove[i-j]=genesV[j]
              linesToRemove[(i-length(scoresV))+j]=(i-length(scoresV))+j
            }
          }else{
            seedsToRemove[i-j] = currentSeed
          }
        }
        currentSeed = screen[,seedColName][i]
        scoresV=c()
        genesV=c()
      } 
      scoresV = c(scoresV, screen[,scoreColName][i])
      genesV = c(genesV, screen[,geneColName][i])
    }
    for(j in 1:length(scoresV)){
      scores2V = scoresV[-j]
      if(length(scores2V) >= min_siRNAs_x_seed ){
        if(
          (!useMedian && (mean(scores2V) < lower_bound_threshold || mean(scores2V) > higher_bound_threshold) ) ||
            (useMedian && (median(scores2V) < lower_bound_threshold || median(scores2V) > higher_bound_threshold))
        )
        {
          genesToRemove[i-j]=genesV[j]
          linesToRemove[(i-length(scoresV))+j]=(i-length(scoresV))+j
        }
      }else {
        seedsToRemove[i-j] = currentSeed
      }
    }
    screen = screen[-linesToRemove[!is.na(linesToRemove)],]
    screen = arrange(screen, previous_order_temp)
    screen = delColDf(screen, "previous_order_temp")
    if(removeGenes){ 
      if(remove_unrepresented_seeds) {
        genesOfSeedsToRemove = unique(subset(screen, screen[,seedColName] %in% seedsToRemove )[,geneColName])
        genesToRemove=c(genesToRemove, genesOfSeedsToRemove)
      }
      screen=subset(screen, !(screen[,geneColName] %in% genesToRemove ))
    }else{ 
      if(remove_unrepresented_seeds) screen=subset(screen, !(screen[,seedColName] %in% seedsToRemove )) 
      else{}
    }
    if(!is.null(min_oligos_x_gene_threshold)){
      genesToRemove2 = sqldf(paste('select ', geneColName, ', count(', geneColName, ') from screen group by ', geneColName, ' having count(', geneColName,') <', min_oligos_x_gene_threshold, sep=""))[, geneColName]
      screen = subset(screen, !(screen[, geneColName] %in% genesToRemove2) )
    }
    return(screen)
  }
  
}

