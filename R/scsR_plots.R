benchmark_shared_hits = function(glA, glB, col, avoidIntersectL=FALSE, 
                                 output_file=NULL, npoints=400, title="", scaleAXPoint = 1, 
                                 scaleBXPoint = NULL, fixedBXPoint=400, displayRandomMultipleLines=TRUE, 
                                 nrandom=20, intersectGenes=TRUE, visualize_pval=FALSE, max_ylim=NULL, xlab=NULL, ylab="shared hits"){
  
  thlog=200
  
  if(!is.null(fixedBXPoint) && !is.null(scaleBXPoint)){
    cat("ERROR: Either fixedBXPoint or scaleBXPoint must be set to NULL.  ")
    stop()
  } 
  if(length(glA)!=length(glB) && length(glB)!=1 ) {
    cat("ERROR: The length of the two input lists (glA and glB) differs ! Please check the input")
    stop()
  }
  if(length(glA)!=length(avoidIntersectL) && length(avoidIntersectL)!=1 ) {
    cat("ERROR: The length of avoidIntersectL differs with the length of glA ! Please check the input")
    stop()
  }
  if(length(glA)!=length(col)) {
    cat("ERROR: The length of col differs with the length of glA ! Please check the input")
    stop()
  }
  
  
  sAgl = list()
  sBgl = list()
  for(i in 1:length(glA)){
    if(length(glB)==1){ sAgl[[i]] = unique(glA[[i]][glA[[i]] %in% glB[[1]]]) }
    else{ sAgl[[i]] = unique(glA[[i]][glA[[i]] %in% glB[[i]]]) }
    
    if( (length(avoidIntersectL)>1 && !avoidIntersectL[i]) || (length(avoidIntersectL)==1 && !avoidIntersectL[1]) ){
      if(length(glB)==1 ){
        sBgl[[i]] = unique(glB[[1]][glB[[1]] %in% glA[[i]]]) 
      }else{sBgl[[i]] = unique(glB[[i]][glB[[i]] %in% glA[[i]]]) }
    }else{ 
      if(length(glB)==1 ){
        sBgl[[i]] = unique(glB[[1]]) 
      }else{
        sBgl[[i]] = unique(glB[[i]]) 
      }
    }
    
  }
  
  if(intersectGenes){
    intersectionSet = sBgl[[1]]
    for(i in 1:length(glA)){
      if((length(avoidIntersectL)>1 && !avoidIntersectL[i]) || (length(avoidIntersectL)==1 && !avoidIntersectL[1])){
        intersectionSet = intersectAll(intersectionSet, sAgl[[i]]) 
      }
    }
    for(i in 1:length(glB)){ intersectionSet = intersectAll(intersectionSet, sBgl[[i]]) }
    for(i in 1:length(glA)){
      sAgl[[i]] = sAgl[[i]][sAgl[[i]] %in% intersectionSet]
      sBgl[[i]] = sBgl[[i]][sBgl[[i]] %in% intersectionSet]
    }
  }
  
  backgroundForHyper=NULL
  if(exists("intersectionSet")) backgroundForHyper=intersectionSet
  
  rmat=matrix(NA,(nrandom*2), length( sAgl[[1]] ))
  for(i in 1:(nrandom*2)){  rmat[i,] = sample( sAgl[[1]] ) }
  
  bvrowsNum = length(glA)
  if(displayRandomMultipleLines){ 
    bvrowsNum = bvrowsNum + nrandom
  }else{bvrowsNum = bvrowsNum + 1}
  bv = matrix(NA,npoints,bvrowsNum)
  for(i in 1:npoints){
    
    if(is.null(fixedBXPoint)) {
      bxpoint = (i*scaleBXPoint)
    }else{bxpoint = fixedBXPoint}
    
    for(j in 1:length(glA)){
      if(length(sAgl[[j]])>=i*scaleAXPoint){
        bv[i,j] = length(intersect( sAgl[[j]][1:(i*scaleAXPoint)] , sBgl[[j]][1:bxpoint] ))
        if(visualize_pval) bv[i,j] = min(thlog, (-1)*log(enrichment_geneSet(sAgl[[j]][1:(i*scaleAXPoint)], sBgl[[j]][1:bxpoint], backgroundForHyper)))
      }
    }
    
    if(!is.null(nrandom)){
      tempV = c()
      for(j in 1:(nrandom)){
        tempInter =  length(intersect( rmat[j,][1:(i*scaleAXPoint)] , rmat[(j+(nrandom)),][1:bxpoint] ))
        tempV=c(tempV, tempInter)
        if(displayRandomMultipleLines) bv[i,(j+length(sAgl))] = tempInter
      }
      if(!displayRandomMultipleLines) bv[i,(1+length(sAgl))] = floor(mean(tempV))
    }
  }
  
  maxYval = max(bv, na.rm=TRUE)
  if(!is.null(max_ylim)){maxYval=max_ylim}
  gSeq = seq(1:npoints)
  if(fixedBXPoint) {gSeq=seq(1,npoints*scaleAXPoint, by=scaleAXPoint)}
  
  
  if(!is.null(output_file)) pdf(output_file)
  if(is.null(fixedBXPoint) & is.null(xlab)){  xlab=paste("number of best hits compared (A = x * ", scaleAXPoint, ", B = x * ", scaleBXPoint, ")" , sep="")
  }else{ if(is.null(xlab)) {xlab=paste("number of best hits compared (A = x ", ", B = ", fixedBXPoint, ")" , sep="") }}
  plot(x=gSeq, y=bv[,1], xlab=xlab, ylab=ylab, main=title, ylim=c(0,maxYval), col=col[1], type="l")
  
  if(!is.null(nrandom) && !visualize_pval){
    grayScale = 80/nrandom
    if(!displayRandomMultipleLines){
      points(x=gSeq, y=bv[,ncol(bv)], col="gray48", pch = 8)
    }else{
      for(i in (length(glA)+1):ncol(bv)){
        grayCol = floor((i-(length(glA)+1))*grayScale) + 10
        grayColStr = paste("gray", grayCol, sep="")
        lines(x=gSeq, y=bv[,i], col=grayColStr, lty=3)
      }
    }
  }
  for(i in 2:length(glA)){ lines(x=gSeq, y=bv[,i], col=col[i]) }
  if(visualize_pval) abline(h=(-1)*log(0.05), col="gray")  
  if(!is.null(output_file)) dev.off()
  
  
  
}


enrichment_heatmap = function(genesVectors, vectorsNames, output_file=NULL, title="", limit=400, species_ncbi_taxonomy_id=9606,
                              enrichmentType="Process", limitMultiPicture=NULL, fdr_threshold=0.05, pvalue_threshold=NULL, 
                              cexRow=NULL, cexCol=1, STRINGversion="9_05", selectTermsVector=NULL, iea = TRUE, sortingMethod="rowMeans", avoidIntersect=FALSE){
  
  string_db <- STRINGdb$new( version=STRINGversion, species=species_ncbi_taxonomy_id, score_threshold=0, input_directory="" )
  enrichList = list()
  
  if(!avoidIntersect) {
    genes=genesVectors[[1]]
    for(i in 2:length(genesVectors)){ genes=intersect(genes, genesVectors[[i]])  }
    for(i in 1:length(genesVectors)){genesVectors[[i]] = intersect(genesVectors[[i]], genes)}
    cat("INFO: we are benchmarking ",length(genes), "genes\n")
  }else{
    cat("WARNING: Your input has not been intersected; hence the datasets could not be perfectly comparable.\n")
    for(i in 1:length(genesVectors)){genesVectors[[i]] = unique(genesVectors[[i]])}
  }
  
  for(i in 1:length(genesVectors)){
    stringGenes = string_db$mp(genesVectors[[i]][0:limit])
    enrichList[[i]] = string_db$get_enrichment(stringGenes, category = enrichmentType, methodMT = "fdr", iea = iea )
  }
  enrichHash = hash()
  for(i in 1:length(genesVectors)){
    for(j in 1:nrow(enrichList[[i]])){
      if((is.null(pvalue_threshold) || enrichList[[i]]$pvalue[j] <= pvalue_threshold) &&
           (is.null(fdr_threshold) || enrichList[[i]]$pvalue_fdr[j] <= fdr_threshold) 
      ){
        enrichVect = rep(NA, length(genesVectors))
        myterm = as.character(enrichList[[i]]$term_description[j])
        
        enter = FALSE
        if(!is.null(selectTermsVector) ){ 
          for(term in selectTermsVector){if(grepl(term, myterm)){enter=TRUE}}
        }else{enter = TRUE}
        
        if(enter){
          if(has.key(myterm, enrichHash)){  enrichVect = enrichHash[[myterm]] }
          enrichVect[i] = -log(enrichList[[i]]$pvalue_fdr[j])
          enrichHash[[myterm]] = enrichVect
        }
      }
    }
  }
  
  enrichMatr = matrix(NA, length(enrichHash), length(genesVectors))
  rownames(enrichMatr) = rep(NA, length(keys(enrichHash)))
  if(length(keys(enrichHash))>0){
    for(i in 1:length(keys(enrichHash))){
      k=keys(enrichHash)[i]
      enrichMatr[i,] = enrichHash[[k]]
      rownames(enrichMatr)[i] = k
    }
  }
  colnames(enrichMatr) = vectorsNames
  
  
  if(is.null(limitMultiPicture)){limitMultiPicture= nrow(enrichMatr)+2}
  if(!is.null(sortingMethod) && sortingMethod=="rowMeans" && nrow(enrichMatr)>1){
    enrichMatr = enrichMatr[order(rowMeans(enrichMatr, na.rm=TRUE), decreasing=TRUE),]
  }
  
  if(is.null(cexRow)){
    if(nrow(enrichMatr) <= 60 ) {
      cexRow=0.4 + 1/log(nrow(enrichMatr))
    }else if(nrow(enrichMatr) <= 120 ) {
      cexRow=0.2 + 1/log(nrow(enrichMatr))
    }else if(nrow(enrichMatr) <= 200 ) {
      cexRow=0.1 + 1/log2(nrow(enrichMatr))
    }else{
      cexRow=0.1
    }
  
  }
  
  lmat = rbind(c(0,3),c(2,1),c(0,4))
  lwid = c(0.7,5)
  lhei = c(1,5,0.7)
  
  if(nrow(enrichMatr)%%limitMultiPicture == 1 && nrow(enrichMatr)!=1) limitMultiPicture=limitMultiPicture+1
  if(nrow(enrichMatr)>1){
    if(!is.null(output_file)) pdf(output_file, paper="a4", width=12, height=12, pagecentre=FALSE)
    suppressWarnings(par(new=TRUE))
    for(i in 1:ceiling(nrow(enrichMatr)/limitMultiPicture)){
      enrichMatrTemp=enrichMatr[(limitMultiPicture*(i-1)+1):min(limitMultiPicture*(i), nrow(enrichMatr)),]
      #heatmap(enrichMatrTemp, Rowv=NA, Colv = NA, col = brewer.pal(6,"Blues"), scale = "none",
      #         margins = c(7,30),  cexRow=cexRow, cexCol=cexCol, main=paste("                                     ",title, sep=""))
      tempV = c()
      for(j in 1:ncol(enrichMatrTemp)){tempV = c(tempV, enrichMatrTemp[,j])}
      tempV=unique(tempV[!is.na(tempV)])
      
      if(length(tempV)==1){
        enrichMatrTemp[is.na(enrichMatrTemp)] <- 0 
      }  
      suppressWarnings(heatmap.2(enrichMatrTemp, density.info="none", trace="none", keysize=1,lmat = lmat, lhei=lhei, lwid=lwid,Rowv=NA, Colv = NA, col = brewer.pal(6,"Blues"), scale = "none",
                margins = c(7,30),  cexRow=cexRow, cexCol=1, main=paste("                ",title, sep="")))
    }  
    if(!is.null(output_file)) dev.off()
    suppressWarnings(par(new=FALSE))
  }else if(nrow(enrichMatr)==1){
    cat("Only one term has been found to be significantly enriched:\n")
    print(enrichMatr)
  }else{cat("No enriched terms below the p-value threshold.")}
  return(enrichMatr)
}


plot_screen_hits = function(screen, output_file=NULL, geneScoreColName="median", seedColName="seed7",  scoreColName="score", geneColName="GeneID", gene_interval = c(1,100), 
                            min_oligos_x_gene=4, min_oligos_x_statistics=4, random=FALSE, kolmogorovSampleSize=5000, ylab="score", xlab="gene", ylim=c(-4,4), 
                            graph_highest_count_thr=16, progress_bar=FALSE){
#   tempDf = subset(suppressWarnings(get_seed_oligos_df(screen, seedColName=seedColName,  scoreColName=scoreColName, geneColName=geneColName, gene_interval = gene_interval, 
#                                                       min_oligos_x_gene=min_oligos_x_gene, min_oligos_x_statistics=min_oligos_x_statistics, 
#                                                       random=random, quiet=quiet, kolmogorovSampleSize=kolmogorovSampleSize)), !is.na(seed_oligos_mean))
  #initialize to not get useless complains from R CMD CHECK 
  seed_oligos_mean <-NULL
  seed_oligos_count<-NULL
  seed_log_pval_ks<-NULL
  geneRank<-NULL
  tempDf =suppressWarnings(get_seed_oligos_df(screen, seedColName=seedColName,  scoreColName=scoreColName, geneColName=geneColName, gene_interval = gene_interval, 
                                              min_oligos_x_gene=min_oligos_x_gene, min_oligos_x_statistics=min_oligos_x_statistics, 
                                              random=random, kolmogorovSampleSize=kolmogorovSampleSize, progress_bar=progress_bar ))
  tempDf$seed_oligos_mean[is.na(tempDf$seed_oligos_mean)] = 0
  for(i in 1:nrow(tempDf)){tempDf$seed_oligos_count[i] = min(tempDf$seed_oligos_count[i], graph_highest_count_thr)}
  none <- element_blank()
  myPlot = suppressMessages(qplot(geneRank, seed_oligos_mean, data=tempDf, xlim=gene_interval, ylim=ylim, na.rm=TRUE,
                                  ylab=ylab, xlab=xlab, size=seed_oligos_count, colour=seed_log_pval_ks )  + 
                              scale_colour_gradient2(midpoint=0, low="darkblue", mid="grey", high="darkred") +
                              geom_point(mapping = aes(geneRank, tempDf[,geneScoreColName], size=NULL), data=tempDf, shape=8, colour=(1)) +
                              geom_vline(aes(xintercept = tempDf$geneRank), size=0.1, alpha=0.3 ) + 
                              scale_x_discrete(breaks = gene_interval[1]:gene_interval[2], labels=unique(tempDf[,geneColName]) ) + 
                              theme(text = element_text(size=8), axis.text.x=element_text(angle=-90, hjust = 0)) +
                              theme(panel.background = element_rect(fill='white', colour='black')) +
                              theme(panel.grid.major = none, panel.grid.minor = none))
  
  if(!is.null(output_file)) ggsave(filename=output_file, plot=myPlot, width=12.65, height=6.18)
  else myPlot
}


plot_screen_seeds_count = function(screen, seedColName="seed7", scoreColName="score", output_file=NULL){
  sa=seeds_analysis(screen, seedColName=seedColName, scoreColName=scoreColName)
  sam=merge(screen, sa, all.x=TRUE, by.x=seedColName, by.y=seedColName)
  if(!is.null(output_file)) pdf(output_file)
  plot(arrange(sam, -count)$count, xlab="oligo", ylab="number of oligos with the same seed ", ylim=c(0,40))
  if(!is.null(output_file)) dev.off()
}


plot_seeds_oligo_count = function(screen, seedColName="seed7", scoreColName="score", output_file=NULL){
  sa=seeds_analysis(screen, seedColName=seedColName, scoreColName=scoreColName)
  if(!is.null(output_file)) pdf(output_file)
  plot(arrange(sa, -count)$count, xlab="seed", ylab="number of oligos", ylim=c(0,40))
  if(!is.null(output_file)) dev.off()
}


plot_seed_score_sd = function(screen, minOligosXSeed=8, seedColName = "seed7", scoreColName = "score", output_file=NULL){
  seeds = seeds_analysis(screen, seedColName=seedColName, scoreColName=scoreColName)
  seedsWidespread = subset(seeds, count>=minOligosXSeed)
  lineSdFix <- lm( seedsWidespread$sd ~ seedsWidespread$score  )
  if(!is.null(output_file)) pdf(output_file)
  plot(seedsWidespread$score, seedsWidespread$sd, ylim=c(0,2.5), xlab="average score of the seeds' oligos", ylab="seeds' oligos standard deviation")
  abline(lineSdFix, col="red")
  mtext(paste("Pearson correlation = ", round(cor.test(seedsWidespread$score, seedsWidespread$sd)$estimate, 5)))
  mtext(paste("Spearman correlation = ", round(cor.test(seedsWidespread$score, seedsWidespread$sd, method="spearman")$estimate, 5)),line = 1)
  if(!is.null(output_file)) dev.off()
}


plot_effective_seeds_head = function(screen, seedColName="seed7",  scoreColName="score", enhancer_analysis=FALSE,
                                     min_oligos_x_seed=10, number_of_seeds=20 , output_file=NULL, color="#CCCCCC33", colorBG="#0000CC11",
                                     xlim=c(-4,4), title=""){
  
  screenRand = screen
  screenRand[,scoreColName]=sample(screen[,scoreColName])
  seedsAnalysisDf = seeds_analysis(screen, seedColName=seedColName, scoreColName=scoreColName, minCount=min_oligos_x_seed, enhancer_analysis=enhancer_analysis)
  seedsAnalysisDf = arrange(seedsAnalysisDf, score)
  seedsAnalysisDfRand = seeds_analysis(screenRand, seedColName=seedColName, scoreColName=scoreColName, minCount=min_oligos_x_seed, enhancer_analysis=enhancer_analysis)
  seedsAnalysisDfRand = arrange(seedsAnalysisDfRand, score)
  
  dataFrame1_h =head(seedsAnalysisDf, number_of_seeds)
  dataFrame2_h =head(seedsAnalysisDfRand, number_of_seeds)
  
  if(!is.null(output_file))
    if(regexpr("pdf", output_file)>0) pdf(output_file)
  else if(regexpr("jpg", output_file)>0) jpeg(output_file)
  
  if(xlim[1]<0){ 
    barplot(dataFrame1_h[scoreColName][,1], names.arg=dataFrame1_h[seedColName][,1], horiz=TRUE, cex.names=0.7, las=1, density=200, lwd = 1, 
            space=0.8, xlim=xlim, col=c(color), main=title )
    par(new=TRUE)
    barplot(dataFrame2_h[scoreColName][,1],  horiz=TRUE, cex.names=0.7, las=1, density=200, lwd = 1, space=0.8, xlim=xlim, col=c(colorBG) )
    #axis(2, at=((1:number_of_seeds)*1.81)-0.6, las=2, lwd = 1, lheight=3, cex.axis=0.7, pos=0, col="green", labels=dataFrame2_h[seedColName][,1]) 
    legend("topleft", c("seed real", "seed random"), fill=c(color, colorBG))
  }else{  
    barplot(dataFrame2_h[scoreColName][,1],  horiz=TRUE, cex.names=0.7, las=1, density=200, lwd = 1, space=0.8, xlim=xlim, col=c(colorBG), main=title )
    par(new=TRUE)
    barplot(dataFrame1_h[scoreColName][,1], names.arg=dataFrame1_h[seedColName][,1], horiz=TRUE, cex.names=0.7, las=1, density=200, lwd = 1, space=0.8, xlim=xlim, col=c(color) )
    #axis(4, at=((1:number_of_seeds)*1.81)-0.6, labels=dataFrame2_h[seedColName][,1], las=2, lwd = 1, lheight=3, cex.axis=0.7, pos=1.5, col="green") 
    legend("topright", c("seed real", "seed random"), fill=c(color, colorBG))
  }
  
  if(!is.null(output_file)) dev.off()
}
