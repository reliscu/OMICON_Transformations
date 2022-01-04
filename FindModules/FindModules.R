## FindModules_0.99 introduces an alternative merging parameter based on the maximum difference of the clustering coefficient between seed features of initial modules vs. seed features of merged modules.  minMEcor is replaced with "merge.by" ("ME" or "CC") and "merge.param" parameters (along with "export.merge.comp" parameter to export scatterplot comparing these two initial module similarity measures for each network).
## FindModules_0.98 parallelizes the calculation of the Jaccard similarity measure using the parallelDist package (also requires Rcpp, RcppParallel, and RcppArmadillo packages).
## FindModules_0.97 introduces an argument ('writeModSnap', default=T) to turn off export of the module snapshot PDFs ("..._Modules.pdf"), since these can be prohibitively huge for huge (e.g. SC) datasets.
## FindModules_0.96 corrects calculation of the Jaccard similarity measure, which previously reported dissimilarity.
## FindModules_0.95 corrects a bug that mislabels the header of column1 in the datkme tables.
## FindModules_0.94 incorporates a new simType ("Jaccard", which is equivalent to 1 - the asymmetric binary dissimilarity calculated by 'dist' with method = "binary"). This should be used for binary data (e.g. some ChIP-seq applications).
## FindModules_0.93 was updated from FM_0.90. Two changes: i) fixed minor bug: (from: ## if(!all.equal(colnames(datExpr1),colnames(simMat))){ ## to: ## if(all.equal(colnames(datExpr1),colnames(simMat))!=TRUE){
## ii) changed ylim form module snapshots (top 15) from ymin-1.5 to ymin-3.
## FindModules_0.90 fixed a bug with the correction of negative values.
## FindModules_0.89 will check for gene expression values = 0; if present, user will have the option to scale all expression values so min = 1 (for entire dataset).
## FindModules_0.88 will check for missing data (NAs) and give the user the option of setting NA = 1 (recommended for RNAseq data). We choose 1 instead of 0 so as not to wreck log transformations.
## FindModules_0.88 will also check for and exclude features (probes, genes, etc.) with 0 variance, instead of retaining these features in kME tables with kME=0. 
## A list of features with 0 variance will be exported as a .csv file.
## FindModules_0.88 has also modified the RemoveNegs function from flagging genes <0 to flagging genes <=0 for module snapshots.
## FindModules_0.88 also can be used locally or on cluster.  On cluster, the "_Modules" directory is no longer automatically created in "/big-data/Shared/", but rather in the R working directory on the cluster.
## FindModules_0.87 has changed the signed vs. unsigned AdjMat calculations for overlapType=="TO".
## FindModules_0.80 is the successor to FindModules_0.65.  This version is meant to take advantage of faster calculations for the similarity matrix.
## Specifically, FindModules_0.65 was written to use the (slow) stats::cor function for calculating Pearson correlations using the default R BLAS.
## This version uses the fast correlation calculations in the WGCNA R package and also allows for WGCNA to use hyperthreading.
## For best results (on a Mac), use in conjunction with the vecLib R BLAS.
## Also note: bicor is now viable as a similarity measure.
## Compared to 0.80, 0.85 does two things: i) adds two checks to ensure that the dimensions and order of a loaded similarity matrix match the expression data, and 
## ii) add another argument (calcBigModStat) to exclude calculation of module statistics for expanded (TopModPosBC) module definitions (highly desirable for datasets with large [>1000] sample sizes).
## 0.86 makes minor modifications for cluster compatibility (time stamp, .libPaths).
## 0.86 also increases the maximum initial module size from 435 (standardColors()) to 657 (colors()).  
## Note that we have standardColors appear in the colorvec first (in the same order as before), followed by the non-standardColors (i.e. the grey colors).

.libPaths("/home/shared/R/x86_64-pc-linux-gnu-library/4.0")

library(WGCNA)
library(flashClust)
library(cluster)
library(svMisc)
library(Biobase)
library(lattice) ## 3D plots
library(qvalue)
library(ellipse) ## plotcorr function
library(purrr) ## required for merging by CC

allowWGCNAThreads()

FindModules=function(
  projectname,
  expr,
  geneinfo,
  sampleindex,
  samplegroups=NULL,
  subset=NULL,
  simMat=NULL,
  saveSimMat=FALSE,
  simType=c("Pearson","Spearman","Bicor","MI","Jaccard"),
  overlapType=c("None","TO","Pearson","Spearman","Bicor"),
  TOtype=c("none","signed","unsigned"),
  TOdenom=c("none","min","mean"),
  beta=1,
  MIestimator="mi.mm",
  MIdisc="equalfreq",
  signumType=c("abs","rel"),
  iterate=FALSE,
  signumvec=c(.999,.99,.98),
  minsizevec=c(8,10,12),
  signum,
  minSize,
  merge.by=c("ME","CC"),
  merge.param,
  export.merge.comp=T,
  ZNCcut=2,
  calcSW=FALSE,
  loadTree=FALSE,
  writeKME=FALSE,
  calcBigModStat=FALSE,
  writeModSnap=TRUE
){
	
	timestamp()
	
## Check to ensure that column 1 contains unique identifiers:
	
	if(length(unique(expr[,1]))!=length(expr[,1])){
		
		stop("Identifiers in column 1 must be unique!")
		
	}
	
## Create root directory, if necessary:
	
	if(length(grep(paste(projectname,"_Modules",sep=""),getwd()))==1){
		
		breakdir=strsplit(getwd(),split="/")
		BNroot=c(1:length(breakdir[[1]]))[is.element(breakdir[[1]],paste(projectname,"_Modules",sep=""))]
		setwd(paste(breakdir[[1]][1:BNroot],collapse="/"))
		
	} else {
		
		dir.create(paste(projectname,"_Modules",sep=""))
		setwd(paste(getwd(),"/",projectname,"_Modules",sep=""))
		
	}
	
	BNrootDir=getwd()	
	
    ## Reorder rows as needed:
    
	if(class(expr[,1])=="integer"|class(expr[,1])=="numeric"){
		
		expr=expr[order(as.numeric(as.character(expr[,1]))),]
		subset=subset[order(as.numeric(as.character(expr[,1])))]
		
	} else {
		
		expr=expr[order(as.character(expr[,1])),]
		subset=subset[order(as.character(expr[,1]))]
		
	}
	
## Check for missing data and set NA = 0 per user (recommended for RNAseq data):

	missing.values=sum(is.na(expr[,sampleindex]))
	total.values=nrow(expr)*length(sampleindex)
	
	if(missing.values>0){
	  cat("\n")
		print(paste(signif(missing.values/total.values*100,2),"% of data are missing. Setting NA = 0 (recommended for RNAseq data)",sep=""),quote=F)
		cat("\n")
		expr[,sampleindex][is.na(expr[,sampleindex])]=0
	}
	
	## Check for values <= 0, and scale so minimum values = 1 (avoids issues with log tranformations)
	
	zero.values=sum(expr[,sampleindex]<=0)
	
	if(zero.values>0){
	  
	  cat("\n")
	  print("Gene expression values <= 0 are present. Scaling all data so minimum expression = 1",quote=F)
	  cat("\n")
	  expr[,sampleindex]=expr[,sampleindex]+abs(min(expr[,sampleindex],na.rm=T))+1
	  
	} 
	
## Note: FM 0.88 excludes genes with 0 variance from kME table:

	datExprall=t(expr[,sampleindex])
	colnames(datExprall)=as.character(expr[,1])
	var.all=apply(datExprall,2,var)
	datExprall=datExprall[,var.all>0]
	
	if(sum(var.all==0)>0){
		
		cat("\n")
		print(paste("Note: ",sum(var.all==0)," features had 0 variance and were excluded",sep=""))
		cat("\n")
		ExcludedFeatures=expr[var.all==0,geneinfo]
		write.table(ExcludedFeatures,file="Excluded_genes_with_0_variance.csv",sep=",",row.names=F,col.names=T)
		
	}
		
	if(!is.null(subset)){
		
		subset=as.logical(subset)
		datExpr1=t(expr[subset,sampleindex])
		colnames(datExpr1)=as.character(expr[subset,1])
		
	} else {
		
		datExpr1=t(expr[,sampleindex])
		colnames(datExpr1)=as.character(expr[,1])
		
	}
	
	var1=apply(datExpr1,2,var)
	
	if(sum(var1==0)>0){
		
		datExpr1=datExpr1[,var1>0]
		
	}
	
	if(is.null(samplegroups)){
		
		samplegroups=c(rep(1,length(sampleindex)))
		
	} else {
		
		if(is.character(samplegroups)){
			
			newsamplegroups=c()
			
			for(z in c(1:length(samplegroups))){
				
				newsamplegroups=c(newsamplegroups,rep(z,length(grep(samplegroups[z],rownames(datExpr1),fixed=T))))
				
			}
			
			samplegroups=newsamplegroups
			
		}
		
		if(sum(is.na(samplegroups))>0){
			
			samplegroups=as.character(samplegroups)
			samplegroups[is.na(samplegroups)]="None"
			
		}
		
	}
	
## findoverlap function:
	
	findoverlap=function(IDs,expr1){
		
		colnames(expr1) %in% IDs
		
	} ## end of findoverlap function
	
## findMEs function:
	
	findMEs=function(overlap1,expr1){
		
		moduleEigengenes(expr=expr1[,overlap1],colors=rep("turqoise",length(overlap1[overlap1])))$eigengenes
		
	} ## end of find MEs function

## findCCs function:
	
	findCCs=function(adjmat,mod){
	  
	  overlap1=is.element(colnames(adjmat),mod)
	  clusterCoef(adjmat[overlap1,overlap1])
	  
	} ## end of find CCs function
		
## mergenodes function:
	
	mergenodes=function(cutreeMods1,tempmods1){
		
		sort(unique(unlist(tempmods1[c(cutreeMods1)])))
		
	}
	
## Calculate modularity and other network statistics:
	
	Modularity=function(datexpr1,modules1){
		
		ModSeeds=unique(unlist(modules1))
		
		ModLabels=c()
		
		for(s in c(1:length(modules1))){
			
			ModLabels=c(ModLabels,rep(names(modules1)[[s]],length(modules1[[s]])))
			
		}
				
			modexpr1=datexpr1[,is.element(colnames(datexpr1),ModSeeds)]
			modadj1=((cor(modexpr1)+1)/2)^2
			diag(modadj1)=0
			sumadj1=sum(modadj1)
			degree1=rowSums(modadj1,na.rm=T)
			rows=row(modadj1)
			cols=col(modadj1)
			modmatrix=matrix(nrow=dim(modadj1)[1],ncol=dim(modadj1)[2],data=NA)
			rownames(modmatrix)=rownames(modadj1)
			colnames(modmatrix)=colnames(modadj1)
		
			for(t in c(1:length(modules1))){
				
				rest1=is.element(ModLabels,names(modules1)[t])
				whichmodseeds=ModSeeds[rest1]
				modmatrix[c(1:dim(modmatrix)[1])[is.element(rownames(modmatrix),whichmodseeds)],c(1:dim(modmatrix)[2])[is.element(colnames(modmatrix),whichmodseeds)]]=as.character(names(modules1)[t])
				
			}
			
			adjmods=data.frame(ID1=colnames(modadj1)[rows[upper.tri(modadj1)]],
							   ID2=colnames(modadj1)[cols[upper.tri(modadj1)]],
							   Adj=modadj1[upper.tri(modadj1)],
							   ModuleSeed=modmatrix[upper.tri(modmatrix)]  
							   )
					
		Neighbors=rep(0,length(adjmods[,1]))
		Neighbors[!is.na(adjmods$ModuleSeed)]=1
		
		findDegree=function(IDvec,degreevec){
			
			degreeout=degreevec[match(IDvec,names(degreevec))]
			
			}
		
		adjmods=data.frame(adjmods,Neighbors,Degree1=findDegree(IDvec=adjmods$ID1,degreevec=degree1),Degree2=findDegree(IDvec=adjmods$ID2,degreevec=degree1))
		adjmods=data.frame(adjmods,Qbase=(adjmods$Adj-((adjmods$Degree1*adjmods$Degree2)/sumadj1))*adjmods$Neighbors)
		ModDegrees=c()
		
		for(u in c(1:length(modules1))){
			
			ModDegrees=c(ModDegrees,sum(degree1[is.element(names(degree1),modules1[[u]])]))
			
			}
		
## Shelving modularity calculations for now; I think the formulas below may be wrong (module order not the same in tapply result versus ModDegrees below)?
		NGMod=NA
#NGMod=sum((tapply(adjmods$Adj,adjmods$ModuleSeed,sum)/(sumadj1/2))-(ModDegrees/sumadj1)^2)
## NGMod reference: Eq.1 from "Performance of modularity maximization in practical contexts" by Good et al. (2010)
		
#GMod=sum(adjmods$Qbase)/sumadj1
## GMod reference: Eq. 7 from "Analysis of weighted networks" by Newman (2004)
		GMod=NA
		
		MeanQbase=mean(adjmods$Qbase[adjmods$Qbase!=0])
		adjmodsout=list(adjmods,NGMod,GMod,MeanQbase)
		names(adjmodsout)=c("Pairwise","NGMod","GMod","MeanQbase")
		adjmodsout
		
	}
		
## Homogeneity score: average correlation b/t ME and module members:
	
	HScore=function(overlap1,overlap2,ME1=MEdf,expr1){
		
		mean(cor(ME1[,overlap2],expr1[,overlap1]))
		
	}
		
## Separation score: 1 - weighted average of signed correlation b/t each ME and every other ME:
	
	SSfx=function(x){
		
		sum(x)/(length(x)-1)
	
	}
	
## PC1 variance explained:
	
	MEPC1=function(overlap1,expr1){
		
		moduleEigengenes(expr=expr1[,overlap1],colors=rep("turqoise",length(overlap1[overlap1])))$varExplained
		
	}
	
## GetValues: retrieve (log-transformed) expression values for all modules.
	
	GetValues=function(overlap1,expr1){
		
		allValues=unlist(log(expr1[,overlap1],2))
		sapply(allValues,c)
		
	}
	
## stderr1: calculate standard error:
	
	stderr1=function(x){sqrt(var(x,na.rm=T)/sum(!is.na(x)))}
	
## MEV: calculate the (log2) mean variance for genes in a module:
	
	MEV=function(overlap1,expr1){
		
		if(length(overlap1[overlap1])==0){
			
			datout=rbind(NA,NA)
			
		}
		
		if(length(overlap1[overlap1])==1){
			
			datout=data.frame(
						  var(log(expr1[,overlap1],2),na.rm=T),
						  stderr1(log(expr1[,overlap1],2))
						  )
			
		}
		
		if(length(overlap1[overlap1])>1){
		
			datout=data.frame(
				   mean(apply(log(expr1[,overlap1],2),2,var,na.rm=T),na.rm=T),
				   stderr1(apply(log(expr1[,overlap1],2),2,var,na.rm=T))
				   )
			
		}
		
		datout
		
	}
	
## MMspec: find the probe sets that are most strongly associated with each module (thresholded by p-value):
	
	CalcMMspecScores=function(kmetable1,seq1a,seq2a,pvalcut1a,direction=c("pos","neg")){
		
		tempdf=data.frame(allcors=unlist(kmetable1[,seq1a]),allpvals=unlist(kmetable1[,seq2a]))
		tempdf=tempdf[order(-tempdf[,1]),]
		tempdf=data.frame(tempdf,Signif=tempdf[,2]<pvalcut1a)
		corthresh=min(abs(tempdf$allcors)[tempdf$Signif])
		
		topmod=function(df1,allmods){
			modno=df1[1]
			topmod1=as.character(allmods[modno])
			topmod1
		}
		
		if(direction=="pos"){
			
			df=data.frame(maxcors=unlist(apply(kmetable1[,seq1a],1,which.max)))
			restUP=apply(kmetable1[,seq1a],1,function(x,y){any(x>y)},y=corthresh)
			df$TopMod=NA
			df$TopMod[restUP]=apply(df[restUP,],1,topmod,allmods=colnames(kmetable1[,seq1a]))
			
		}
		
		if(direction=="neg"){
			
			df=data.frame(mincors=apply(kmetable1[,seq1a],1,which.min))
			restDOWN=apply(kmetable1[,seq1a],1,function(x,y){any(x<(-y))},y=corthresh)
			df$TopMod=NA
			df$TopMod[restDOWN]=apply(df[restDOWN,],1,topmod,allmods=colnames(kmetable1[,seq1a]))
			
		}
		
		gsub("kME","",df$TopMod)
		
	}
		
## TopModIDs will find expanded module sizes or IDs below kme p-value cutoff (denominators for homegrown module specificity calculation):
	
	TopModIDs=function(datkme1,geneinfo1,pvalcut1,justcount){
		
		seq1=seq((length(geneinfo1)+1),dim(datkme1)[2],2)
		seq2=seq((length(geneinfo1)+2),dim(datkme1)[2],2)
		
		bigmods=vector(mode="list",length=length(seq1))
		
		if(justcount==TRUE){
		
			for(i in c(1:length(seq1))){
			
				bigmods[[i]]=length(datkme1[,1][datkme1[,seq1[i]]>0&datkme1[,seq2[i]]<pvalcut1])
			
				collectGarbage()
				
			}
			
		} else {
				
			for(i in c(1:length(seq1))){
					
				bigmods[[i]]=datkme1[,1][datkme1[,seq1[i]]>0&datkme1[,seq2[i]]<pvalcut1]
					
			}
					
		}
		
			bigmods
		
	}
	
## reachFreq: calculate the coverage (reach) and redundancy (frequency) of probe sets assigned to modules.
## Coverage is defined as a) % of all probe sets assigned to a module (Q<.05); b) % of probe sets with Z.K > user-supplied cutoff (e.g. ZNCcut=2) assigned to a module; 
## c) % of probe sets with Z.C > cutoff assigned to a module; d) % of probe sets with Z.K & Z.C > cutoff assigned to a module.
## Redundancy is defined as above, but for % of all probe sets in each case assigned to more than one module.
	
	reachFreq=function(expr1,subset1,ZNCcut1=ZNCcut,modcount1,tstamp1){
		
		modcount1DF=data.frame(table(modcount1))
		modcount1DF[,1]=as.numeric(as.character(modcount1DF[,1]))
		ALL=sum(modcount1DF$Freq)
		ALLreachPcnt=sum(modcount1DF$Freq[modcount1DF$modcount1>0])/sum(modcount1DF$Freq)*100
		ALLredundPcnt=sum(modcount1DF$Freq[modcount1DF$modcount1>1])/sum(modcount1DF$Freq)*100
		ALLratio=(sum(modcount1DF$Freq[modcount1DF$modcount1>0])/sum(modcount1DF$Freq))/(sum(modcount1DF$Freq[modcount1DF$modcount1>1])/sum(modcount1DF$Freq))
		
		reachFreqOut=data.frame(ALL=ALL,ALLreachPcnt=ALLreachPcnt,ALLredundPcnt=ALLredundPcnt,ALLratio=ALLratio)
		
		if(length(grep("Z.K",colnames(expr1)))==1){
			
			tempKDF=data.frame(subset1,expr1$Z.K>ZNCcut1,modcount1,modcount1>0,modcount1>1)
			ZKuniverse=length(tempKDF[,2][!is.na(tempKDF[,2])&tempKDF[,2]])
			ZKreach=sum(apply(tempKDF[,c(1,2,4)],1,all))
			ZKfreq=sum(apply(tempKDF[,c(1,2,5)],1,all))
			reachFreqOut=data.frame(reachFreqOut,
									TopZK=ZKuniverse,
									ZKreachPcnt=(ZKreach/ZKuniverse*100),
									ZKredundPcnt=(ZKfreq/ZKuniverse*100),
									ZKratio=(ZKreach/ZKuniverse)/(ZKfreq/ZKuniverse))
			
		}
		
		if(length(grep("Z.C",colnames(expr1)))==1){
			
			tempCDF=data.frame(subset1,expr1$Z.C>ZNCcut1,modcount1,modcount1>0,modcount1>1)
			ZCuniverse=length(tempCDF[,2][!is.na(tempCDF[,2])&tempCDF[,2]])
			ZCreach=sum(apply(tempCDF[,c(1,2,4)],1,all))
			ZCfreq=sum(apply(tempCDF[,c(1,2,5)],1,all))
			reachFreqOut=data.frame(reachFreqOut,
									TopZC=ZCuniverse,
									ZCreachPcnt=(ZCreach/ZCuniverse*100),
									ZCredundPcnt=(ZCfreq/ZCuniverse*100),
									ZCratio=(ZCreach/ZCuniverse)/(ZCfreq/ZCuniverse))
			
		}
		
		if(length(grep("Z.K",colnames(expr1)))+length(grep("Z.C",colnames(expr1)))==2){
			
			ZCKuniverse=length(tempCDF[,2][!is.na(tempCDF[,2])&tempCDF[,2]&!is.na(tempKDF[,2])&tempKDF[,2]])
			ZCKreach=sum(apply(data.frame(tempCDF[,c(1,2)],tempKDF[,c(2,4)]),1,all))
			ZCKfreq=sum(apply(data.frame(tempCDF[,c(1,2)],tempKDF[,c(2,5)]),1,all))
			reachFreqOut=data.frame(reachFreqOut,
									TopZCK=ZCKuniverse,
									ZCKreachPcnt=(ZCKreach/ZCKuniverse*100),
									ZCKredundPcnt=(ZCKfreq/ZCKuniverse*100),
									ZCKratio=(ZCKreach/ZCKuniverse)/(ZCKfreq/ZCKuniverse)
									)
			
		}
		
		plotvec=sort(c(grep("reach",colnames(reachFreqOut)),grep("redund",colnames(reachFreqOut))))
		
		if(length(plotvec)==2){
			
			main1=""
			
		} else {
			
		main1=paste("Z > ",ZNCcut1,sep="")
			
		}
		
		pdf(file=paste("ReachRedundancyAnalysis_",tstamp1,".pdf",sep=""))
		par(mar=c(7,5,4,2))
		barplot(unlist(reachFreqOut[1,plotvec]),
				names=colnames(reachFreqOut)[plotvec],
				las=3,
				ylab="Percentage",
				cex.lab=1.4,
				main=main1,
				cex.main=1.8,
				ylim=c(0,100))
		colseq=seq(1,length(reachFreqOut[1,]),4)
		universes=paste(colnames(reachFreqOut)[colseq]," (",reachFreqOut[1,colseq],")",sep="")
		mtext(paste(universes,collapse="; "))
		dev.off()
		
		reachFreqOut
		
	}
	
	
## plot.mat function for making heat maps (by Sandrine Dudoit, modified slightly)
	
	plot.mat=function (x, nrgcols = 50, rlabels = FALSE, clabels = FALSE, 
					   rcols = 1, ccols = 1, title = "", ...)
	{
		n <- nrow(x)
		p <- ncol(x)
		image(1:p, 1:n, t(x[n:1, ]), col = rgcolors.func(nrgcols), 
			  axes = FALSE, xlab = "", ylab = "", ...)
		if (length(ccols) == 1) {
			axis(1, at = 1:p, labels = clabels, las = 2, cex.axis = 0.6, 
				 col.axis = ccols)
		}
		if (length(ccols) == p) {
			cols <- unique(ccols)
			for (i in 1:length(cols)) {
				which <- (1:p)[ccols == cols[i]]
				axis(1, at = which, labels = clabels[which], las = 2, 
					 cex.axis = 0.6, col.axis = cols[i])
			}
		}
		if (length(rcols) == 1) {
			axis(2, at = n:1, labels = rlabels, las = 2, cex.axis = 0.6, 
				 col.axis = rcols)
		}
		if (length(rcols) == n) {
			cols <- unique(rcols)
			for (i in 1:length(cols)) {
				which <- (1:n)[rcols == cols[i]]
				axis(2, at = (n:1)[which], labels = rlabels[which], 
					 las = 2, cex.axis = 0.6, col.axis = cols[i])
			}
		}
		mtext(title, side = 3, line = 2, cex=1.25, font=2)
		box()
		
	} ## end of plot.mat	

## Calculate node similarities:
	
	if(!is.null(simMat)){
		
		print("Loading similarity matrix...")
		timestamp()
		load(simMat)
		timestamp()
		
		if(dim(simMat)[2]!=dim(datExpr1)[2]){
			
            stop("Dimensions of simMat and datExpr1 do not match!")
            
		}
		
		if(all.equal(colnames(datExpr1),colnames(simMat))!=TRUE){
            
			stop("Columns of simMat and datExpr1 are out of order!")
			
		}
		
	} else {
	
	print("Calculating node similarities...")
	cat("\n")
						
	if(simType=="Pearson"){
		
		simMat=cor(datExpr1,method="p",use="p")
		collectGarbage()
		diag(simMat)=0
		collectGarbage()
		
	}
				
## Not sure why, but Spearman with "use='p'" does not seem to work;
		
	if(simType=="Spearman"){
		
		simMat=cor(datExpr1,method="s",use="p")
		collectGarbage()
		diag(simMat)=0
		collectGarbage()
		
	}
	
## Note: re-test Bicor after ATLAS/WGCNA situation is resolved.
		
	if(simType=="Bicor"){
			
		simMat=bicor(datExpr1,use="p")
		collectGarbage()
		diag(simMat)=0
		collectGarbage()
		
	}
	
## Note: "Jaccard" calculates the asymmetric binary similarity measure as implemented by the 'dist' function (method = 'binary').
## However, the dist function is too slow for large matrices, so we use the parallelized version (parDist).
## For future reference we can also access additional proximity measures through the 'proxy' R package.
## This can be used when feature data are binary.
		
	if(simType=="Jaccard"){
			
      library(Rcpp)
      library(RcppParallel)
      library(RcppArmadillo)
      library(parallelDist)
            
      ## Note: expression data must be binarized prior to calculation.
      ## By default, the minimum expression value in the matrix is set to 0 and anything greater is set to 1.
      ## But you could also let the user choose a threshold based on quantiles (may need to sample from expression space, as below); implement later.
            
      datExprBin=datExpr1
            
      ## Check for negative values and scale if present to avoid problems:
            
      if(sum(datExprBin<0)>0){
                
          datExprBin=datExprBin+abs(min(datExprBin))
                
      }
            
      datExprBin[datExprBin==min(datExprBin)]=0
      datExprBin[datExprBin>0]=1
            
      simMat=1-as.matrix(parDist(x=t(datExprBin),method="binary"))
			diag(simMat)=0
			collectGarbage()
			
		} ## end of if(simType=="Jaccard"){
		
	if(simType=="MI"){

## WARNING! Non-linear operation...time increases exponentially with # features.  Do not recommend more than ~3000 features or so [revisit with ATLAS BLAS].
		
		timestamp()
		print("Starting build.mim...")
		cat("\n")
		mim=build.mim(datExpr1,estimator=MIestimator,disc=MIdisc)
		timestamp()
		colnames(mim)=as.character(colnames(datExpr1))
		rownames(mim)=as.character(colnames(datExpr1))
		collectGarbage()
		print("Starting mrnet...")
		cat("\n")
		timestamp()
		simMat=mrnet(mim)
		timestamp()
		collectGarbage()
		
		if(max(simMat,na.rm=T)>1){
			
			simMat=simMat/max(simMat)
		
		}
	
		diag(simMat)=0
		collectGarbage()
				
	} ## end of if(simType=="MI"){
	
	if(overlapType=="Pearson"){
			
		simMat=cor(simMat,method="p",use="p")
		collectGarbage()
		diag(simMat)=0
		collectGarbage()
			
		}
			
## Note: CANNOT USE TO WITH ATLAS BLAS RIGHT NOW.
		
	if(overlapType=="TO"){
		
		if(TOtype=="signed"){
			
			 AdjMat=((simMat+1)/2)^beta
			
		} else {
		
		  AdjMat=(abs(simMat)^(beta/2))*sign(simMat)
			
	  }
		
	  collectGarbage()
	  diag(AdjMat)=1
# Note: current TOMdist requires diag(AdjMat)=1
	  rm(simMat)
		collectGarbage()
		timestamp()
		disSimMat=TOMdist(adjMat=AdjMat,TOMDenom=TOdenom,TOMType=TOtype)
		timestamp()
		rm(AdjMat)
		collectGarbage()
		colnames(disSimMat)=colnames(datExpr1)
		rownames(disSimMat)=colnames(datExpr1)
		collectGarbage()
		
		if(loadTree==FALSE){		
		
		  cluster1=flashClust(as.dist(disSimMat),method="complete")
			collectGarbage()
			save(cluster1,file=paste(simType,"-",overlapType,"_",TOtype,"_",TOdenom,"_p",beta,"_",projectname,"_",length(cluster1$order),"_clustering",sep=""))
      dendroname=paste(simType,"-",overlapType,"_",TOtype,"_",TOdenom,"_p",beta,"_",projectname,"_",length(cluster1$order),"_dendrogram.pdf",sep="")
  	  pdf(file=dendroname,bg="transparent",width=24,height=12)
			plot(cluster1,labels=F,main=paste(simType," ",overlapType," ",TOtype," ",TOdenom,", p",beta," ",length(cluster1$order)," features",sep=""))
			dev.off()
			
		} else {
			
			 load(paste(simType,"-",overlapType,"_",TOtype,"_",TOdenom,"_p",beta,"_",projectname,"_",length(simMat[,1]),"_clustering",sep=""))
			
		}
		
		simMat=1-disSimMat
		diag(simMat)=0
		rm(disSimMat)
		collectGarbage()
		
	  } ## end of if(overlapType=="TO"){
		
	  timestamp()
	  print("...done")
	  cat("\n")
	  collectGarbage()
		
## Note: saving large matrices takes forever.  Loading takes a while too.  Seems faster to just recalculate every time (at least for Pearson;
## definitely not true for MI).
		
	  if(saveSimMat==TRUE){
		
		  timestamp()
		  print("Saving similarity matrix...")
		  cat("\n")
		  save(simMat,file=paste(simType,"_",projectname,"_",dim(simMat)[1],sep=""))
		  timestamp()
		
		}
	
	} ## end of if(!is.null(simMat))
		
	collectGarbage()
	
## Load clustering (optional):
	
	if(loadTree==FALSE){
		
		cluster1=flashClust(as.dist(1-simMat),method="complete")
		collectGarbage()
		save(cluster1,file=paste(simType,"_p",beta,"_",projectname,"_",length(cluster1$order),"_clustering",sep=""))
		dendroname=paste(simType,"_p",beta,"_",projectname,"_",length(cluster1$order),"_dendrogram.pdf",sep="")
		pdf(file=dendroname,bg="transparent",width=24,height=12)
		plot(cluster1,labels=F,main=paste(simType,", p",beta," ",length(cluster1$order)," features",sep=""))
		dev.off()
				
	} else {
					
		load(paste(simType,"_p",beta,"_",projectname,"_",length(simMat[,1]),"_clustering",sep=""))
		
	}
			
	standardcolors2=standardColors()
	colorvec=standardcolors2[c(7,6,2,5,3,1,8:(length(unique(samplegroups))+1))]
		
## Iterate (optional):
		
	if(iterate==TRUE){
		
		print("Testing parameters...")
		cat("\n")
		
		if(signumType=="rel"){
		
			if(((dim(datExpr1)[2]*(dim(datExpr1)[2]-1))/2)>1e+06){
			
				relsignumvec=quantile(sample(simMat,1000000,replace=FALSE),probs=signumvec)
				collectGarbage()
			
			} else {
		
				relsignumvec=quantile(simMat[upper.tri(simMat)],probs=signumvec)
				collectGarbage()
			
			}
			
		}
		
		minSizevec=minsizevec
				
		if(merge.by=="CC"){
		  
		  if(simType=="Pearson"|simType=="Spearman"|simType=="Bicor"){
		    
		    adjMat=(simMat+1)/2
		    diag(adjMat)=0
		    
		  }
		  
		  if(simType=="Jaccard"){
		    
		    adjMat=simMat
		    diag(adjMat)=0
		    
		  }
		  
		  if(overlapType=="Pearson"|simType=="Spearman"|simType=="Bicor"){
		    
		    adjMat=(simMat+1)/2
		    diag(adjMat)=0
		    
		  }
		  
		  if(overlapType=="TO"){
		    
		    adjMat=simMat
		    diag(adjMat=0)
		    
		  }
		  
		} ## end of if(merge.by=="CC"){
		
		rm(simMat)
		collectGarbage()
		
		NQS=data.frame(Network=0)
		NQS=NQS[-c(1),]
		count=0

		for(i in c(1:length(signumvec))){
			
			if(signumType=="rel"){
			
				signum=relsignumvec[i]
				
			} else {
				
				signum=signumvec[i]
				
			}
			
			cutree1=cutree(cluster1,h=1-signum)
			
			keptmodsDF=data.frame(table(cutree1))
			
			for(j in c(1:length(minSizevec))){
				
				tstamp1=format(Sys.time(), "%X")
				tstamp1=gsub(" AM","",tstamp1)
				tstamp1=gsub(" PM","",tstamp1)
				tstamp1=gsub(":","-",tstamp1)
				
				if(overlapType=="TO"){
					
					fileheader=paste(simType,"-",overlapType,"_",TOtype,"_",TOdenom,"_signum",signif(signum,3),"_minSize",minSizevec[j],"_merge_",merge.by,"_",merge.param,"_",length(datExpr1[1,]),sep="")
					
				} else {
					
					#fileheader=paste(simType,"-",overlapType,"_signum",signif(signum,3),"_minSize",minSizevec[j],"_merge_",merge.by,"_",merge.param,"_",length(datExpr1[1,]),sep="")
				  fileheader=paste(simType,"-",overlapType,"_signum",signif(signumvec[i],3),"_minSize",minSizevec[j],"_merge_",merge.by,"_",merge.param,sep="")
				  
				}
				
				dir.create(fileheader)
				setwd(paste(getwd(),"/",fileheader,sep=""))
				
				count=count+1

## Find initial modules:
				
				minSize=minSizevec[j]
				keptmodscutree=cutree1[cutree1 %in% keptmodsDF$cutree1[keptmodsDF$Freq>=(minSize)]]
				modules=tapply(as.character(names(keptmodscutree)),factor(keptmodscutree),list)

## Currently the max # of modules is limited by the length of the "colors()" vector (n=657); may need to revisit module naming convention down the road.
				
				if(length(modules)>0&length(modules)<658){
					
					initialModules=length(modules)
					print(paste("signum = ",signif(signum,3),", minSize = ",minSize,": ",length(modules)," initial modules found",sep=""))
					cat("\n")
					
## Calculate initial module eigengenes and pairwise correlations:
					
					allMEs=mapply(FUN=findMEs,overlap1=mapply(FUN=findoverlap,IDs=modules,MoreArgs=list(expr1=datExpr1),SIMPLIFY=FALSE),MoreArgs=list(expr1=datExpr1),SIMPLIFY=TRUE)
					MEdf=matrix(unlist(allMEs),ncol=length(modules))
					MEdf=as.data.frame(MEdf)
					corMEdf=cor(MEdf)
					diag(corMEdf)=0
					collectGarbage()

## Merge similar modules:
					
					if(length(modules)>1){ 
					
					  print("Merging similar modules...")
					  cat("\n")
					
					  if(merge.by=="CC"){
					  
					    ## Strategy for merging by clustering coefficient differentials:
					
					    names(modules)=paste("mod",c(1:length(modules)),sep="")
					  
					    ## Calculate clustering coefficients for initial module seed features:
					    modseed.CC=mapply(FUN=findCCs,mod=modules,MoreArgs=list(adjmat=adjMat),SIMPLIFY=FALSE)
					  
					    ## Merge all possible pairs of module seed features:
					  	
					    merged.seeds=map(cross2(modules,modules),unlist)
					    merged.seeds=lapply(merged.seeds,sort)
					    temp2=unique(merged.seeds)
					    temp3=lapply(temp2,unique)
					    merged.seeds=temp2[unlist(lapply(temp2,length))==unlist(lapply(temp3,length))]
					
					    ## Calculate clustering coefficients for merged module seed features:
					
					    mergedmodseed.CC=mapply(FUN=findCCs,mod=merged.seeds,MoreArgs=list(adjmat=adjMat),SIMPLIFY=FALSE)
					
					    ## Name the merged modules:
					
					    count.mat=matrix(nrow=length(modules),ncol=length(modules),data=-888)
					    count.mat[lower.tri(count.mat)]=c(1:length(mergedmodseed.CC))
					    count.mat=t(count.mat)
					    rownames(count.mat)=paste("mod",c(1:length(modules)),sep="")
					    colnames(count.mat)=paste("mod",c(1:length(modules)),sep="")
					
					    for(k in 1:length(mergedmodseed.CC)){
					  
					      names(mergedmodseed.CC)[k]=paste(rownames(count.mat)[which(count.mat==k,arr.ind=TRUE)[1]],colnames(count.mat)[which(count.mat==k,arr.ind=TRUE)[2]],sep="|")
					  
					    }
					
					   ## Calculate the max difference between original CCs and CCs from merged modules:
					
					    CC.diff.list=vector(mode="list",length=length(mergedmodseed.CC))
					    names(CC.diff.list)=names(mergedmodseed.CC)
					
					    for(k in 1:length(CC.diff.list)){
					  
					      which.mods=unlist(strsplit(names(mergedmodseed.CC[k]),split="|",fixed=T))
					      CC.diff.list[[k]]=data.frame(CC.init=unlist(modseed.CC[is.element(names(modseed.CC),which.mods)]),CC.merged=NA)
					      rownames(CC.diff.list[[k]])=gsub(paste(which.mods[1],".",sep=""),"",rownames(CC.diff.list[[k]]),fixed=T)
					      rownames(CC.diff.list[[k]])=gsub(paste(which.mods[2],".",sep=""),"",rownames(CC.diff.list[[k]]),fixed=T)
					      CC.diff.list[[k]]=CC.diff.list[[k]][order(rownames(CC.diff.list[[k]])),]

					      if(sum(rownames(CC.diff.list[[k]]) == names(mergedmodseed.CC[[k]])) < length(rownames(CC.diff.list[[k]]))){
					    
					        stop(paste("Feature labels do not align for merged module ",k," when comparing clustering coefficients",sep=""))
					    
  					    } else {
					    
					        CC.diff.list[[k]]$CC.merged=mergedmodseed.CC[[k]]
					    
					      }
					  
					      CC.diff.list[[k]]=data.frame(CC.diff.list[[k]],CC.diff=CC.diff.list[[k]][,1]-CC.diff.list[[k]][,2])
					  
					    } ## end of for(k in 1:length(CC.diff.list)){
					
					   ## Now we will create a module similarity matrix by calculating 1 - the  mean CC.diff for every pair of merged modules:
					
					    CC.diff.max=unlist(lapply(CC.diff.list,function(x){max(x$CC.diff)}))
					    names(CC.diff.max)=names(CC.diff.list)
					
					    CCdf=matrix(nrow=length(modules),ncol=length(modules),data=-888)
					    CCdf[lower.tri(CCdf)]=1-CC.diff.max
					    CCdf=t(CCdf)
					    CCdf[lower.tri(CCdf)]=1-CC.diff.max
					    diag(CCdf)=0
					    rownames(CCdf)=paste("mod",c(1:length(modules)),sep="")
					    colnames(CCdf)=paste("mod",c(1:length(modules)),sep="")
					  
					    if(export.merge.comp==T){
					    
					      pdf(file=paste(simType,"-",overlapType,"_signum",signif(signum,3),"_minSize",minSizevec[j],"_compare_merge_parameters.pdf",sep=""))
					      plot(corMEdf[upper.tri(corMEdf)],CCdf[upper.tri(CCdf)],xlab="cor(ME)",ylab="1-max(CC.diff)",main=paste(simType,"-",overlapType,"_signum",signif(signum,3),"_minSize",minSizevec[j],sep=""),xlim=c(-1,1))
					      mtext("Initial module similarities",side=3)
					      abline(v=c(0.85,0.9,0.95))
					      abline(h=c(0.95,0.96,0.97))
					      dev.off()
					    
					    }
					  
					    while(max(CCdf)>=(merge.param)&length(modules)>1){
					    
					      clusterMod=flashClust(as.dist(1-CCdf),method="complete")
					      cutreeMods=cutree(clusterMod,h=1-merge.param)
					      modules=mapply(FUN=mergenodes,cutreeMods1=tapply(c(1:length(cutreeMods)),factor(cutreeMods),list),MoreArgs=list(tempmods1=modules),SIMPLIFY=FALSE)
					    
					      if(length(modules)>1){
					      
					        names(modules)=paste("mod",c(1:length(modules)),sep="")
					    
					        ## Calculate clustering coefficients for initial module seed features:
					        modseed.CC=mapply(FUN=findCCs,mod=modules,MoreArgs=list(adjmat=adjMat),SIMPLIFY=FALSE)
					    
					        ## Merge all possible pairs of module seed features:
					    
					        merged.seeds=map(cross2(modules,modules),unlist)
					        merged.seeds=lapply(merged.seeds,sort)
					        temp2=unique(merged.seeds)
					        temp3=lapply(temp2,unique)
					        merged.seeds=temp2[unlist(lapply(temp2,length))==unlist(lapply(temp3,length))]
					      
					        ## Calculate clustering coefficients for merged module seed features:
					    
					        mergedmodseed.CC=mapply(FUN=findCCs,mod=merged.seeds,MoreArgs=list(adjmat=adjMat),SIMPLIFY=FALSE)
					    
					        ## Name the merged modules:
					    
					        count.mat=matrix(nrow=length(modules),ncol=length(modules),data=-888)
					        count.mat[lower.tri(count.mat)]=c(1:length(mergedmodseed.CC))
					        count.mat=t(count.mat)
					        rownames(count.mat)=paste("mod",c(1:length(modules)),sep="")
					        colnames(count.mat)=paste("mod",c(1:length(modules)),sep="")
					    
					        for(k in 1:length(mergedmodseed.CC)){
					      
					          names(mergedmodseed.CC)[k]=paste(rownames(count.mat)[which(count.mat==k,arr.ind=TRUE)[1]],colnames(count.mat)[which(count.mat==k,arr.ind=TRUE)[2]],sep="|")
					      
					        }
					        
					        ## Calculate the max difference between original CCs and CCs from merged modules:
					        
					        CC.diff.list=vector(mode="list",length=length(mergedmodseed.CC))
					        names(CC.diff.list)=names(mergedmodseed.CC)
					        
					        for(k in 1:length(CC.diff.list)){
					          
					          which.mods=unlist(strsplit(names(mergedmodseed.CC[k]),split="|",fixed=T))
					          CC.diff.list[[k]]=data.frame(CC.init=unlist(modseed.CC[is.element(names(modseed.CC),which.mods)]),CC.merged=NA)
					          rownames(CC.diff.list[[k]])=gsub(paste(which.mods[1],".",sep=""),"",rownames(CC.diff.list[[k]]),fixed=T)
					          rownames(CC.diff.list[[k]])=gsub(paste(which.mods[2],".",sep=""),"",rownames(CC.diff.list[[k]]),fixed=T)
					          CC.diff.list[[k]]=CC.diff.list[[k]][order(rownames(CC.diff.list[[k]])),]
					          
					          if(sum(rownames(CC.diff.list[[k]]) == names(mergedmodseed.CC[[k]])) < length(rownames(CC.diff.list[[k]]))){
					            
					            stop(paste("Feature labels do not align for merged module ",k," when comparing clustering coefficients",sep=""))
					            
					          } else {
					            
					            CC.diff.list[[k]]$CC.merged=mergedmodseed.CC[[k]]
					            
					          }
					          
					          CC.diff.list[[k]]=data.frame(CC.diff.list[[k]],CC.diff=CC.diff.list[[k]][,1]-CC.diff.list[[k]][,2])
					          
					        } ## end of for(k in 1:length(CC.diff.list)){
					        
					        ## Now we will create a module similarity matrix by calculating 1 - the max CC.diff for every pair of merged modules:
					        
					        CC.diff.max=unlist(lapply(CC.diff.list,function(x){max(x$CC.diff)}))
					        names(CC.diff.max)=names(CC.diff.list)
					        
					        CCdf=matrix(nrow=length(modules),ncol=length(modules),data=-888)
					        CCdf[lower.tri(CCdf)]=1-CC.diff.max
					        CCdf=t(CCdf)
					        CCdf[lower.tri(CCdf)]=1-CC.diff.max
					        diag(CCdf)=0
					        rownames(CCdf)=paste("mod",c(1:length(modules)),sep="")
					        colnames(CCdf)=paste("mod",c(1:length(modules)),sep="")
					    
					      } ## end of if(length(modules)>1
					            
					    } ## end of while(max(CCdf)>=(merge.param)&length(modules)>1){
					        
					    ## Calculate module eigengenes:
					  
					    allMEs=mapply(FUN=findMEs,overlap1=mapply(FUN=findoverlap,IDs=modules,MoreArgs=list(expr1=datExpr1),SIMPLIFY=FALSE),MoreArgs=list(expr1=datExpr1),SIMPLIFY=TRUE)
					    MEdf=matrix(unlist(allMEs),ncol=length(modules))
					    MEdf=as.data.frame(MEdf)
					  
					  } ## end of if(merge.by=="CC"){
					
					  if(merge.by=="ME"){
					  
					    while(max(corMEdf)>=(merge.param)){
					    
					      clusterMod=flashClust(as.dist(1-corMEdf),method="complete")
					      cutreeMods=cutree(clusterMod,h=1-merge.param)
					      modules=mapply(FUN=mergenodes,cutreeMods1=tapply(c(1:length(cutreeMods)),factor(cutreeMods),list),MoreArgs=list(tempmods1=modules),SIMPLIFY=FALSE)
					      allMEs=mapply(FUN=findMEs,overlap1=mapply(FUN=findoverlap,IDs=modules,MoreArgs=list(expr1=datExpr1),SIMPLIFY=FALSE),MoreArgs=list(expr1=datExpr1),SIMPLIFY=TRUE)
					      MEdf=matrix(unlist(allMEs),ncol=length(modules))
					      MEdf=as.data.frame(MEdf)
					      corMEdf=cor(MEdf)
					      diag(corMEdf)=0
					      collectGarbage()
					    
					    }
					  
					  } ## end of if(merge.by=="ME"){
					  
					} ## end of if(length(modules)>1){
					
					finalModules=length(modules)
					print(paste(length(modules),"modules found"))
					cat("\n")

## Label modules:
					
					modulesize=c()
					
					for(m in c(1:dim(MEdf)[2])){
						
						modulesize=c(modulesize,length(modules[[m]]))
						
					}
					
					standardcolors=c(standardColors(),setdiff(colors(),standardColors()))
					standardcolors=standardcolors[1:length(modules)]
					names(MEdf)=standardcolors[rank(-modulesize,ties.method="first")]
					names(modules)=names(MEdf)
					corMEdf=cor(MEdf)
					diag(corMEdf)=0

## Re-order modules:
					
					if(length(modules)>1){
					
						whichMethod="average"
						MEcluster=flashClust(as.dist(1-corMEdf),method=whichMethod)
						MEdf=MEdf[,MEcluster$order]
						corMEdf=cor(MEdf)
						diag(corMEdf)=0
						modules=modules[MEcluster$order]
						modulesize=modulesize[MEcluster$order]
						
					}
					
## Export network summary:
					
					print("Exporting network summary...")
					cat("\n")
					
					MEdfout=data.frame(colnames(expr)[sampleindex],MEdf)
					colnames(MEdfout)[1]="Sample"
					
					write.table(MEdfout,file=paste("Module_eigengenes_",tstamp1,".csv",sep=""),sep=",",row.names=F,col.names=T)
										
					if(length(modules)>2){
						
						pdf(file=paste("ME_network_",tstamp1,".pdf",sep=""),width=21,height=12)
						par(mfrow=c(2,1))
						par(mar=c(.5,5,6,2))
						plot(MEcluster,sub="",xlab="",ylab="1 - cor",main=projectname,cex.main=1.5,cex.lab=1.25,cex=0.8)
						mtext(paste(whichMethod," linkage,"," signum = ",signif(signum,3),", minSize = ",minSize,", merge.by = ",merge.by," > ",merge.param,", ",length(datExpr1[1,])," features",sep=""))
						par(mar = c(1,5,6,2))
						plot.mat(MEdf,nrgcols=500,rlabels=colnames(expr)[sampleindex],clabels=colnames(MEdf),rcols=colorvec[samplegroups])
						dev.off()
						
						corMEdfOut=corMEdf
						diag(corMEdfOut)=1
						pdf(file=paste("ME_correlations_",tstamp1,".pdf",sep=""),width=16,height=16)
						colors <- c("#A50F15","#DE2D26","#FB6A4A","#FCAE91","#FEE5D9","white",
									"#EFF3FF","#BDD7E7","#6BAED6","#3182BD","#08519C")
						plotcorr(corMEdfOut,col=colors[5*corMEdfOut+6],main=projectname)
						mtext(paste(whichMethod," linkage,"," signum = ",signif(signum,3),", minSize = ",minSize,", merge.by = ",merge.by," > ",merge.param,", ",length(datExpr1[1,])," features",sep=""))
						dev.off()
						
					} else {
						
						print("Less than three modules; no network summary exported.")
						cat("\n")
						
					}
					
## Calculate kME table:
					
					print("Building kME table...")
					cat("\n")

					collectGarbage()
					kMEtable=data.frame(expr[var.all>0,geneinfo])
					
					modseeds=data.frame(unlist(modules),rep(names(modules),modulesize))
					
					if(class(expr[,1])=="integer"|class(expr[,1])=="numeric"){
						
						modseeds=modseeds[order(as.numeric(as.character(modseeds[,1]))),]
						
					} else {
						
						modseeds=modseeds[order(as.character(modseeds[,1])),]
						
					}
					
					ModSeeds=rep(NA,length(dim(datExprall)[2]))
					ModSeeds[is.element(colnames(datExprall),modseeds[,1])]=as.character(modseeds[is.element(modseeds[,1],colnames(datExprall)),2])
										
					count1=0
					
					for(k in c(1:dim(MEdf)[2])){
						
						kME1=cor(MEdf[,k],datExprall)
						kME1pval=corPvalueStudent(kME1,length(sampleindex))
## Convention: p-values for perfect positive or negative correlations will be set to 1e-300 so as not to break downstream code.  Very, very rare event.
						kME1pval[is.na(kME1pval)]=1e-300
						kMEtable=data.frame(kMEtable,t(kME1),t(kME1pval))
						colnames(kMEtable)[length(geneinfo)+(2*count1)+1]=paste("kME",colnames(MEdf)[k],sep="")
						colnames(kMEtable)[length(geneinfo)+(2*count1)+2]=paste("kME",colnames(MEdf)[k],".pval",sep="")
						count1=count1+1
						
					}
					
					collectGarbage()
					
	## Find p-value for determining module specificity.  Will use Bonferroni p-value as one (more stringent) measure of specificity.
	## Will also use FDR <.05 as another (less stringent) measure of module specificity.
	## Both measures will be exported to kME table; but BC measure will be used for module specificity calculations that are present
	## in module statistics and network statistics outputs.
					
					pvalseq=grep(".pval",colnames(kMEtable))
					corseq=pvalseq-1
					
					allPvalues=c()
					
					if(length(pvalseq)==1){
						
						allPvalues=kMEtable[,pvalseq]
						
					} else {
					
					for(n in c(1:length(pvalseq))){
						
						allPvalues=c(allPvalues,kMEtable[,pvalseq][,n])
						
					}
						
					}
					
					allQvalues=qvalue(allPvalues)
					specPvalCutQ=max(allQvalues$pvalues[allQvalues$qvalues<.05])
					
					specPvalCutB=.05/(length(kMEtable[,1])*length(pvalseq))
						
					collectGarbage()
					
	## Find probe sets that are most strongly associated with each module (thresholded by p-value):
					
					print("Finding top modules for all genes...")
					cat("\n")
					
					pvalseq=grep(".pval",colnames(kMEtable))
					corseq=pvalseq-1
					
					if(length(corseq)==1){
						
						MMspecScoresBpos=rep(NA,length(kMEtable[,1]))
						MMspecScoresQpos=rep(NA,length(kMEtable[,1]))
						MMspecScoresBneg=rep(NA,length(kMEtable[,1]))
						MMspecScoresQneg=rep(NA,length(kMEtable[,1]))
						
					} else {
					
					MMspecScoresBpos=CalcMMspecScores(kmetable1=kMEtable,seq1a=corseq,seq2a=pvalseq,pvalcut1a=specPvalCutB,direction="pos")
					MMspecScoresQpos=CalcMMspecScores(kmetable1=kMEtable,seq1a=corseq,seq2a=pvalseq,pvalcut1a=specPvalCutQ,direction="pos")
					MMspecScoresBneg=CalcMMspecScores(kmetable1=kMEtable,seq1a=corseq,seq2a=pvalseq,pvalcut1a=specPvalCutB,direction="neg")
					MMspecScoresQneg=CalcMMspecScores(kmetable1=kMEtable,seq1a=corseq,seq2a=pvalseq,pvalcut1a=specPvalCutQ,direction="neg")

					}
					
## Calculate mean (log2) expression levels and percentiles:
					
					FeatureMean=apply(log(datExprall,2),2,mean,na.rm=T)
					FeaturePercentile=(FeatureMean/max(FeatureMean,na.rm=T))*100
										
					kMEtable=data.frame(kMEtable[,1:length(geneinfo)],ModSeed=ModSeeds,MeanExpr=FeatureMean,MeanExprPcnt=FeaturePercentile,TopModPosBC=MMspecScoresBpos,TopModPosFDR=MMspecScoresQpos,TopModNegBC=MMspecScoresBneg,TopModNegFDR=MMspecScoresQneg,kMEtable[,(length(geneinfo)+1):dim(kMEtable)[2]])
                    colnames(kMEtable)[1:length(geneinfo)]=colnames(expr)[1:length(geneinfo)]
                    colnames(kMEtable)[length(geneinfo)+2]="MeanExpr"
					colnames(kMEtable)[length(geneinfo)+3]="MeanExprPercentile"
					colnames(kMEtable)[length(geneinfo)+4]=paste("TopModPosBC_",signif(specPvalCutB,3),sep="")
					colnames(kMEtable)[length(geneinfo)+5]=paste("TopModPosFDR_",signif(specPvalCutQ,3),sep="")
					colnames(kMEtable)[length(geneinfo)+6]=paste("TopModNegBC_",signif(specPvalCutB,3),sep="")
					colnames(kMEtable)[length(geneinfo)+7]=paste("TopModNegFDR_",signif(specPvalCutQ,3),sep="")

					if(writeKME==TRUE){
					
						write.table(kMEtable,file=paste("kME_table_",tstamp1,".csv",sep=""),sep=",",row.names=F,col.names=T)
						
					}
					
					collectGarbage()
					
## Export module heat maps, ME barplots, and lineplots:
					
                    if(writeModSnap==TRUE){
                    
                        print("Exporting heat maps & barplots & lineplots...")
                        cat("\n")
					
                        pdf(file=paste("Modules_",tstamp1,".pdf",sep=""),width=10,height=7.5)
                        nf=layout(matrix(c(1,3,2,3), 2, 2, byrow=TRUE), respect=FALSE)
					
                        for(p in c(1:dim(MEdf)[2])){
						
                            restmod=is.element(colnames(datExprall),modules[[p]])
                            datmod=datExprall[,restmod]
                            kmerank=rank(-cor(MEdf[,p],datmod))
                            datmod=datmod[,order(kmerank)]
                            par(mar=c(0,2,5,2))
                            plot.mat(t(scale(datmod)),title=names(MEdf)[p],cex.main=1.5)
                            par(mar=c(10,2,1,2))
                            barplot1=barplot(MEdf[,p],col=gsub("ME","",names(MEdf)[p]),xaxs="i",yaxs="i")
                            axis(side=1,at=barplot1[,1],labels=F,tick=T)
                            box()
                            mtext(rownames(datmod),1,at=barplot1[,1],adj=1.3,cex=0.7,las=3,col=colorvec[samplegroups])
                            par(mar=c(10,4,5,2))
                            nogenes=15
                            restmod2=is.element(colnames(kMEtable),paste("kME",names(modules[p]),sep=""))
                            overlapTOP=is.element(colnames(datExprall),kMEtable[,1][order(-kMEtable[,restmod2])][1:nogenes])
                            datmod2=datExprall[,overlapTOP]
                            ordering=order(rank(-cor(MEdf[,p],datmod2)))
						
                            datMat1=scale(datmod2[,1:nogenes][,ordering])
                            Gene1=as.character(expr$Gene[is.element(expr[,1],colnames(datMat1))])[ordering]
                            multiples=grep("/",Gene1)
						
                            if(length(multiples)>0){
						
                                for(m in c(1:length(multiples))){
								
                                    Gene1[multiples[m]]=paste(strsplit(Gene1[multiples[m]],split=" /")[[1]][1],"*",sep="")
							
                                }
						
                            }
                            
                            if(max(datMat1)>10){
                              subtr=max(datMat1)/2
                            } else {
                              subtr=5
                            }
						    
                            matplot(data.frame(c(1:dim(datmod2)[1])),data.frame(datMat1),main="Top genes by kME", ylab="Expression level (z-score)", xlab="Sample", type="l", lwd=2, cex.main=1.5, cex.lab=1.25, axes=F,col=standardColors()[1:nogenes],lty=1,ylim=c(min(datMat1)-subtr,max(datMat1)+.1))
                            legend("bottom",Gene1,cex=1,lty=1,col=standardColors()[1:nogenes],lwd=2,ncol=3)
                            axis(side=2,tick=T)
                            box()
					
                        } ## end of for(p in c(1:dim(MEdf)[2]))
					
                        dev.off()
					
                        collectGarbage()
                        
                    } ## end of if(writeModSnap==TRUE){
									
## Export module quality statistics:
					
					print("Calculating module statistics...")
					cat("\n")	
					
	## Calculate expanded module sizes below kme p-value cutoff (denominators for module specificity calculation):
					
					bigmodules=unlist(TopModIDs(datkme1=kMEtable,geneinfo1=c(1:(length(geneinfo)+7)),pvalcut1=specPvalCutB,justcount=TRUE))
					
	## Calculate and assemble module quality statistics (note: currently using ModSeeds for calculations; could also use bigmoduleIDs):	
					
					bigmoduleIDs=TopModIDs(datkme1=kMEtable,geneinfo1=c(1:(length(geneinfo)+7)),pvalcut1=specPvalCutB,justcount=FALSE)
					names(bigmoduleIDs)=names(modules)
					
					MMspecOut=data.frame(table(MMspecScoresBpos))
					noSpecMods=setdiff(names(modules),MMspecOut[,1])
					
					if(length(noSpecMods)>0){
						
						noSpecModsDF=data.frame(MMspecScoresBpos=noSpecMods,Freq=0)
						MMspecOut=rbind(MMspecOut,noSpecModsDF)
						MMspecOut=MMspecOut[order(as.character(MMspecOut[,1])),]
						
					}
					
					MMspecOut2=MMspecOut[order(match(MMspecOut[,1],names(modules))),]
					
					HScores=mapply(FUN=HScore,overlap1=mapply(FUN=findoverlap,IDs=modules,MoreArgs=list(expr1=datExpr1),SIMPLIFY=FALSE),
								   overlap2=mapply(FUN=findoverlap,IDs=names(modules),MoreArgs=list(expr1=MEdf),SIMPLIFY=FALSE),
								   MoreArgs=list(expr1=datExpr1),SIMPLIFY=TRUE)
					
					signCorMEdf=(corMEdf+1)/2
					diag(signCorMEdf)=0
					
					## Note: in a signed network with beta=2, mean Separation should be ~0.5 with bimodal distribution (< 0.5 = positive ME correlations; > 0.5 = negative ME correlations).
					SScores=1-apply(signCorMEdf,1,SSfx)
					
					PC1VE=unlist(mapply(FUN=MEPC1,overlap1=mapply(FUN=findoverlap,IDs=modules,MoreArgs=list(expr1=datExpr1),SIMPLIFY=FALSE),MoreArgs=list(expr1=datExpr1),SIMPLIFY=TRUE))
					
					adjME=((cor(MEdf)+1)/2)^2
					diag(adjME)=0
					
					if(length(modules)>2){
						
						FNCME=fundamentalNetworkConcepts(adjME)
						
						if(length(unique(FNCME$Connectivity))>2&length(unique(FNCME$ClusterCoef))>2){
							
							C.k=cor.test(scale(FNCME$Connectivity),scale(FNCME$ClusterCoef),method="s")$estimate
						
						} else {
							
							C.k=NA
							
						}
						
					} else {
						
						FNCME=data.frame(Connectivity=rep(NA,length(modules)),ClusterCoef=rep(NA,length(modules)))
						C.k=NA
						
					}
					
					overlap1=mapply(FUN=findoverlap,IDs=modules,MoreArgs=list(expr1=datExprall),SIMPLIFY=TRUE)
					
					if(length(grep("Z.K",colnames(expr)))==1){
					
						meanZKmodgenes=c()
						
						for(r in c(1:length(modules))){
							
							meanZKmodgenes=c(meanZKmodgenes,mean(expr$Z.K[overlap1[,r]],na.rm=T))
							
						}
						
					} else {
						
						meanZKmodgenes=rep(NA,length(modules))
						
					}
					
					if(length(grep("Z.C",colnames(expr)))==1){
						
						meanZCmodgenes=c()
						
						for(r in c(1:length(modules))){
							
							meanZCmodgenes=c(meanZCmodgenes,mean(expr$Z.C[overlap1[,r]],na.rm=T))
							
						}
						
					} else {
						
						meanZCmodgenes=rep(NA,length(modules))
						
					}
					
					modValues=mapply(FUN=GetValues,overlap1=mapply(FUN=findoverlap,IDs=modules,MoreArgs=list(expr1=datExprall),SIMPLIFY=FALSE),
									 MoreArgs=list(expr1=datExprall),SIMPLIFY=FALSE)
					
					if(length(modules)>1){
					
						MeanExpr=sapply(modValues,mean,na.rm=T)
						Z.MeanExpr=scale(MeanExpr)
						SEExpr=sapply(modValues,stderr1)
						MEVout=mapply(FUN=MEV,overlap1=mapply(FUN=findoverlap,IDs=modules,MoreArgs=list(expr1=datExprall),SIMPLIFY=FALSE),
								MoreArgs=list(expr1=datExprall),SIMPLIFY=TRUE)
						
						} else {
						
						modValues=unlist(modValues)
						MeanExpr=mean(modValues,na.rm=T)
						Z.MeanExpr=scale(MeanExpr)
						SEExpr=stderr1(modValues)
						overlap1=is.element(colnames(datExprall),modules[[1]])
						MEVout=rbind(mean(apply(log(datExprall[,overlap1],2),2,var,na.rm=T),na.rm=T),stderr1(apply(log(datExprall[,overlap1],2),2,var)))
						
					}
					
					MeanExprVar=unlist(MEVout[1,])
					SEMeanExprVar=unlist(MEVout[2,])					
					Z.MeanExprVar=scale(MeanExprVar)
					
					if(calcBigModStat==TRUE){
					
						bigmodValues=mapply(FUN=GetValues,overlap1=mapply(FUN=findoverlap,IDs=bigmoduleIDs,MoreArgs=list(expr1=datExprall),SIMPLIFY=FALSE),
										MoreArgs=list(expr1=datExprall),SIMPLIFY=FALSE)
					
						if(length(bigmoduleIDs)>1){
					
							MeanExprBCMod=sapply(bigmodValues,mean,na.rm=T)
							Z.MeanExprBCMod=scale(MeanExprBCMod)
							MEVoutBCMod=mapply(FUN=MEV,overlap1=mapply(FUN=findoverlap,IDs=bigmoduleIDs,MoreArgs=list(expr1=datExprall),SIMPLIFY=FALSE),
										   MoreArgs=list(expr1=datExprall),SIMPLIFY=TRUE)
						
							if(sum(is.na(MeanExprBCMod))>0){
					
								SEExprBCMod=rep(NA,length(MeanExprBCMod))
								SEExprBCMod[!is.na(MeanExprBCMod)]=sapply(bigmodValues[!is.na(MeanExprBCMod)],stderr1)
						
								} else {
								
									SEExprBCMod=sapply(bigmodValues,stderr1)
								
								}
						
						} else {
						
							bigmodValues=unlist(bigmodValues)
							MeanExprBCMod=mean(bigmodValues,na.rm=T)
							Z.MeanExprBCMod=scale(MeanExprBCMod)
							SEExprBCMod=stderr1(bigmodValues)
							overlap1=is.element(colnames(datExprall),bigmoduleIDs[[1]])
							MEVoutBCMod=rbind(mean(apply(log(datExprall[,overlap1],2),2,var,na.rm=T),na.rm=T),stderr1(apply(log(datExprall[,overlap1],2),2,var)))
								
							}
						
						MeanExprVarBCMod=unlist(MEVoutBCMod[1,])
						SEMeanExprVarBCMod=unlist(MEVoutBCMod[2,])				
						Z.MeanExprVarBCMod=scale(MeanExprVarBCMod)
						
					} else {
						
						MeanExprBCMod=NA
						Z.MeanExprBCMod=NA
						SEExprBCMod=NA
						MeanExprVarBCMod=NA
						Z.MeanExprVarBCMod=NA
						SEMeanExprVarBCMod=NA
						
					}
					
					MEdfRank=apply(abs(MEdf),2,sort,decreasing=TRUE)
					MEmax=MEdfRank[1,]
					PosOrNeg=apply(MEdf,2,function(x){abs(max(x))>abs(min(x))})
					MEstep1=rep(-88888,length(MEmax))
					MEstep1[PosOrNeg]=MEmax[PosOrNeg]-MEdfRank[2,PosOrNeg]
					MEdfRank=apply(MEdf,2,sort,decreasing=TRUE)
					MEstep1[!PosOrNeg]=MEdfRank[dim(MEdfRank)[1],!PosOrNeg]-MEdfRank[dim(MEdfRank)[1]-1,!PosOrNeg]
					MEstep1=abs(MEstep1)
					
					if(calcSW==TRUE){
						
						QM=data.frame(Module=names(modules),InitialSize=modulesize,ExpandedSize=bigmodules,Unique=MMspecOut2$Freq,Homogeneity=HScores,ZHomogeneity=scale(HScores),Separation=SScores,ZSeparation=scale(SScores),Silhouette=SWScores,ZSilhouette=scale(SWScores),PC1VE=PC1VE)
						QM=data.frame(QM[,1:4],Specificity=QM$Unique/QM$ExpandedSize*100,QM[,5:11],Z.Kmod=scale(FNCME$Connectivity),Z.Cmod=scale(FNCME$ClusterCoef),MeanZKgenemodgenes=meanZKmodgenes,MeanZCmodgenes=meanZCmodgenes,MeanExpr=MeanExpr,ZMeanExpr=Z.MeanExpr,SEExpr=SEExpr,MeanExprVar=MeanExprVar,ZMeanExprVar=Z.MeanExprVar,SEMeanExprVar=SEMeanExprVar,MeanExprBCMod=MeanExprBCMod,ZMeanExprBCMod=Z.MeanExprBCMod,SEExprBCMod=SEExprBCMod,MeanExprVarBCMod=MeanExprVarBCMod,ZMeanExprVarBCMod=Z.MeanExprVarBCMod,SEMeanExprVarBCMod=SEMeanExprVarBCMod,MEmax=MEmax,MEstep1=MEstep1)
						colnames(QM)[grep("Unique",colnames(QM))]=paste("Unique_",signif(specPvalCutB,3),sep="")
						colnames(QM)[grep("ExpandedSize",colnames(QM))]=paste("ExpandedSize_",signif(specPvalCutB,3),sep="")
						
					} else {
						
						QM=data.frame(Module=names(modules),InitialSize=modulesize,ExpandedSize=bigmodules,Unique=MMspecOut2$Freq,Homogeneity=HScores,ZHomogeneity=scale(HScores),Separation=SScores,ZSeparation=scale(SScores),PC1VE=PC1VE)
						QM=data.frame(QM[,1:4],Specificity=QM$Unique/QM$ExpandedSize*100,QM[,5:9],Z.Kmod=scale(FNCME$Connectivity),Z.Cmod=scale(FNCME$ClusterCoef),MeanZKmodgenes=meanZKmodgenes,MeanZCmodgenes=meanZCmodgenes,MeanExpr=MeanExpr,ZMeanExpr=Z.MeanExpr,SEExpr=SEExpr,MeanExprVar=MeanExprVar,ZMeanExprVar=Z.MeanExprVar,SEMeanExprVar=SEMeanExprVar,MeanExprBCMod=MeanExprBCMod,ZMeanExprBCMod=Z.MeanExprBCMod,SEExprBCMod=SEExprBCMod,MeanExprVarBCMod=MeanExprVarBCMod,ZMeanExprVarBCMod=Z.MeanExprVarBCMod,SEMeanExprVarBCMod=SEMeanExprVarBCMod,MEmax=MEmax,MEstep1=MEstep1)
						colnames(QM)[grep("Unique",colnames(QM))]=paste("Unique_",signif(specPvalCutB,3),sep="")
						colnames(QM)[grep("ExpandedSize",colnames(QM))]=paste("ExpandedSize_",signif(specPvalCutB,3),sep="")
						
					}
					
					write.table(QM,file=paste("Module_statistics_",tstamp1,".csv",sep=""),sep=",",row.names=F,col.names=T)
					
					pdf(file=paste("Module_statistics_",tstamp1,".pdf",sep=""),width=12,height=12)
											
						par(mar=c(10,5,4,2))
						barplot(QM[,2],main="Initial module size",col=as.character(QM$Module),ylab="Members",cex.main=2,cex.lab=1.5,names=as.character(QM$Module),las=3)
					
						par(mar=c(10,5,4,2))
						barplot(QM[,3],main=paste("Expanded module (P < ",gsub("ExpandedSize_","",colnames(QM)[3]),") size",sep=""),col=as.character(QM$Module),ylab="Members",cex.main=2,cex.lab=1.5,names=as.character(QM$Module),las=3)
					
						par(mar=c(5,5,4,2))
						plot(QM[,2],QM[,3],xlab="Initial size",ylab=paste(gsub("_"," (P < ",colnames(QM)[3]),")",sep=""),main="Module size",cex.main=2,cex.lab=1.5,pch=19,col=as.character(QM$Module))
						if(dim(QM)[1]>2&length(unique(QM[,2]))>1){
							abline(lm(QM[,3]~QM[,2]),col="red",lwd=2)
						}
						
						par(mar=c(5,5,4,2))
						plot(QM[,3],QM[,4],xlab=paste(gsub("_"," (P < ",colnames(QM)[3]),")",sep=""),ylab=paste(gsub("_"," (P < ",colnames(QM)[4]),")",sep=""),main="Module size and specificity",cex.main=2,cex.lab=1.5,pch=19,col=as.character(QM$Module))
						if(dim(QM)[1]>2&sd(QM[,3])>0&sd(QM[,4])>0){
						abline(lm(QM[,4]~QM[,3]),col="red",lwd=2)
						}
						abline(0,1)
						
						if(!all(is.nan(QM$Specificity))){
							par(mar=c(10,5,4,2))
							barplot(QM$Specificity,main="Module specificity",col=as.character(QM$Module),ylab="Percentage",cex.main=2,cex.lab=1.5,names=as.character(QM$Module),las=3)
						} else {
							par(mar=c(10,5,4,2))
							barplot(rep(0,length(QM$Specificity)),main="Module specificity",col=as.character(QM$Module),ylab="Percentage",cex.main=2,cex.lab=1.5,names=as.character(QM$Module),las=3)
						}
						
						par(mar=c(5,5,4,2))
						plot(QM$Homogeneity,QM$PC1VE,xlab="Homogeneity",ylab="PC1VE",main="Module purity",cex.main=2,cex.lab=1.5,pch=19,col=as.character(QM$Module))
						if(dim(QM)[1]>2&!is.na(lm(QM$PC1VE~QM$Homogeneity)[[1]][2])){
							abline(lm(QM$PC1VE~QM$Homogeneity),col="red",lwd=2)
						}
					
						par(mar=c(5,5,4,2))
						plot(QM$InitialSize,QM$PC1VE,xlab="Initial size",ylab="PC1VE",main="Module purity vs size",cex.main=2,cex.lab=1.5,pch=19,col=as.character(QM$Module))
					
						if(sum(!is.na(QM$Z.Kmod))!=0){
					
							par(mar=c(5,5,4,2))
							plot(QM$Z.Kmod,QM$Z.Cmod,xlab="Module Z.K",ylab="Module Z.C",main="Module cor(K,C)",cex.main=2,cex.lab=1.5,pch=19,col=as.character(QM$Module))
							mtext(paste("rho = ",signif(cor.test(QM$Z.Kmod,QM$Z.Cmod,method="s")$estimate,2),sep=""))
							if(dim(QM)[1]>2){
								abline(lm(QM$Z.Cmod~QM$Z.Kmod),col="red",lwd=2)
							}
					
						}
						
						par(mar=c(10,5,4,2))
						sdME=sd(QM$MeanExpr)
						meanME=mean(QM$MeanExpr)
						barplot(QM$MeanExpr,main="Module expression levels",col=as.character(QM[,1]),ylab="Mean (log2) expression",cex.main=2,cex.lab=1.5,names=as.character(QM$Module),las=3,ylim=c(min(0,min(QM$MeanExpr,na.rm=T)-abs(min(QM$MeanExpr,na.rm=T)*0.1)),max(QM$MeanExpr,na.rm=T)+max(QM$MeanExpr,na.rm=T)*0.1))
						addErrorBars(means=QM$MeanExpr,errors=QM$SEExpr,two.side=TRUE)
						abline(h=meanME,col="black",lwd=2)
						abline(h=meanME+(2*sdME),col="red",lwd=2)
						abline(h=meanME-(2*sdME),col="red",lwd=2)
						
						par(mar=c(10,5,4,2))
						sdMEV=sd(QM$MeanExprVar)
						meanMEV=mean(QM$MeanExprVar)
						barplot(QM$MeanExprVar,main="Module expression variance",col=as.character(QM[,1]),ylab="Mean (log2) variance",cex.main=2,cex.lab=1.5,names=as.character(QM$Module),las=3,ylim=c(min(0,min(QM$MeanExprVar,na.rm=T)-abs(min(QM$MeanExprVar,na.rm=T)*0.1)),max(QM$MeanExprVar,na.rm=T)+max(QM$MeanExprVar,na.rm=T)*0.1))
						addErrorBars(means=QM$MeanExprVar,errors=QM$SEMeanExprVar,two.side=TRUE)
						abline(h=meanMEV,col="black",lwd=2)
						abline(h=meanMEV+(2*sdMEV),col="red",lwd=2)
						abline(h=meanMEV-(2*sdMEV),col="red",lwd=2)
						
						if(!all(is.nan(QM$Specificity))&calcBigModStat==TRUE){
							par(mar=c(10,5,4,2))
							sdME=sd(QM$MeanExprBCMod,na.rm=T)
							meanME=mean(QM$MeanExprBCMod,na.rm=T)
							barplot(QM$MeanExprBCMod,main=paste("Expanded module (P < ",gsub("ExpandedSize_","",colnames(QM)[3]),") expression levels",sep=""),col=as.character(QM[,1]),ylab="Mean (log2) expression",cex.main=2,cex.lab=1.5,names=as.character(QM$Module),las=3,ylim=c(min(0,min(QM$MeanExprBCMod,na.rm=T)-abs(min(QM$MeanExprBCMod,na.rm=T)*0.1)),max(QM$MeanExprBCMod,na.rm=T)+max(QM$MeanExprBCMod,na.rm=T)*0.1))
							addErrorBars(means=QM$MeanExprBCMod,errors=QM$SEExprBCMod,two.side=TRUE)
							abline(h=meanME,col="black",lwd=2)
							abline(h=meanME+(2*sdME),col="red",lwd=2)
							abline(h=meanME-(2*sdME),col="red",lwd=2)
						} else {
							par(mar=c(10,5,4,2))
							barplot(rep(0,length(QM$MeanExprBCMod)),main=paste("Expanded module (P < ",gsub("ExpandedSize_","",colnames(QM)[3]),") expression levels",sep=""),col=as.character(QM[,1]),ylab="Mean (log2) expression",cex.main=2,cex.lab=1.5,names=as.character(QM$Module),las=3)
						}
					
						if(!all(is.nan(QM$Specificity))&calcBigModStat==TRUE){
							par(mar=c(10,5,4,2))
							sdMEV=sd(QM$MeanExprVarBCMod,na.rm=T)
							meanMEV=mean(QM$MeanExprVarBCMod,na.rm=T)
							barplot(QM$MeanExprVarBCMod,main=paste("Expanded module (P < ",gsub("ExpandedSize_","",colnames(QM)[3]),") variance",sep=""),col=as.character(QM[,1]),ylab="Mean (log2) variance",cex.main=2,cex.lab=1.5,names=as.character(QM$Module),las=3,ylim=c(min(0,min(QM$MeanExprVarBCMod,na.rm=T)-abs(min(QM$MeanExprVarBCMod,na.rm=T)*0.1)),max(QM$MeanExprVarBCMod,na.rm=T)+max(QM$MeanExprVarBCMod,na.rm=T)*0.1))
							addErrorBars(means=QM$MeanExprVarBCMod,errors=QM$SEMeanExprVarBCMod,two.side=TRUE)
							abline(h=meanMEV,col="black",lwd=2)
							abline(h=meanMEV+(2*sdMEV),col="red",lwd=2)
							abline(h=meanMEV-(2*sdMEV),col="red",lwd=2)
						} else {
							par(mar=c(10,5,4,2))
							barplot(rep(0,length(QM$MeanExprVarBCMod)),main=paste("Expanded module (P < ",gsub("ExpandedSize_","",colnames(QM)[3]),") variance",sep=""),col=as.character(QM[,1]),ylab="Mean (log2) variance",cex.main=2,cex.lab=1.5,names=as.character(QM$Module),las=3)
						}
							
						par(mar=c(10,5,4,2))
						barplot(QM$MEmax,main="|Module eigengene maximum|",col=as.character(QM$Module),ylab="Highest value",cex.main=2,cex.lab=1.5,names=as.character(QM$Module),las=3)
						
						par(mar=c(10,5,4,2))
						barplot(QM$MEstep1,main="Module eigengene first step",col=as.character(QM$Module),ylab="Difference between top two values",cex.main=2,cex.lab=1.5,names=as.character(QM$Module),las=3)
							
						dev.off()
		
## Export network statistics:				
					
					print("Calculating network statistics...")
					cat("\n")
										
	## Calculate % probe sets that are present in multiple modules and reachFreq:
					
					pvalseq=grep(".pval",colnames(kMEtable))
					corseq=pvalseq-1
					
					LogiMods=matrix(nrow=dim(kMEtable)[1],ncol=length(pvalseq),data=FALSE)
					
					if(length(pvalseq)==1){
						
						LogiMods[,p][kMEtable[,corseq]>0&kMEtable[,pvalseq]<specPvalCutB]=TRUE
						
					} else {
					
					for(p in c(1:length(pvalseq))){
						
						LogiMods[,p][kMEtable[,corseq][,p]>0&kMEtable[,pvalseq][,p]<specPvalCutB]=TRUE
						
						}
						
					}
					
					No.Modules=rowSums(LogiMods,na.rm=T)
					No.ModulesDF=data.frame(table(No.Modules))
					names(No.Modules)=colnames(datExprall)
					
					reachFreqOut=reachFreq(expr1=expr,subset1=subset,ZNCcut1=ZNCcut,modcount1=No.Modules,tstamp1=tstamp1)
					
	##	Calculate modularity and mean Z.K / Z.C:
					
					ModOut=Modularity(datexpr1=datExpr1,modules1=modules)
					
					if(length(grep("Z.K",colnames(expr)))==1){
					
						meanZKmodgenes=mean(expr$Z.K[subset][is.element(expr[subset,1],unique(unlist(modules)))],na.rm=T)
						meanZKbigmodgenes=mean(expr$Z.K[subset][is.element(expr[subset,1],unique(unlist(bigmoduleIDs)))],na.rm=T)
						
						
					} else {
						
						meanZKmodgenes=NA
						meanZKbigmodgenes=NA
						
					}
					
					if(length(grep("Z.C",colnames(expr)))==1){
					
						meanZCmodgenes=mean(expr$Z.C[subset][is.element(expr[subset,1],unique(unlist(modules)))],na.rm=T)
						meanZCbigmodgenes=mean(expr$Z.C[subset][is.element(expr[subset,1],unique(unlist(bigmoduleIDs)))],na.rm=T)
						
					} else {
						
						meanZCmodgenes=NA
						meanZCbigmodgenes=NA
						
					}
					
					collectGarbage()
					
					#NQS=rbind(NQS,data.frame(Network=count,Similarity=simType,RelSignum=signumvec[i],Signum=signum,MinSize=minSize,MergeBy=merge.by,MergeParam=merge.param,InitialModules=initialModules,FinalModules=finalModules,ChangeModules=0,Members=length(unique(unlist(modules))),MeanSpecificity=mean(QM$Specificity),NGModularity=ModOut$NGMod,GModularity=ModOut$GMod,MeanQbase=ModOut$MeanQbase,MeanHomogeneity=mean(HScores),MeanPC1VE=mean(PC1VE),MeanSeparation=mean(SScores),MeanZKmodgenes=meanZKmodgenes,MeanZKbigmodgenes=meanZKbigmodgenes,MeanZCmodgenes=meanZCmodgenes,MeanZCbigmodgenes=meanZCbigmodgenes,corModKC=C.k,reachFreqOut))	
					
					NQS$ChangeModules=(NQS$FinalModules-NQS$InitialModules)/NQS$InitialModules*100
					NQS$ChangeModules[NQS$InitialModules==0]=NA
					NQS=NQS[order(NQS$MinSize,-NQS$Signum),]
	
					#write.table(data.frame(t(NQS)),file=paste(simType,"-",overlapType,"_p",beta,"_",projectname,"_",length(cluster1$order),"_network_statistics_",tstamp1,".csv",sep=""),sep=",",row.names=T,col.names=F)
					write.table(data.frame(t(NQS)),file=paste(simType,"-",overlapType,"_signum",signum[i],"_minSize",minSizevec[j],"_merge_",merge.by,"_",merge.param,"_network_statistics.csv",sep=""),sep=",",row.names=T,col.names=F)
					
				} else { ## end of if(length(modules)>0&length(modules)<436)
				
					KCcount=length(grep("Z.K",colnames(expr)))+length(grep("Z.C",colnames(expr)))
					if(KCcount==0){
						RFcount=1
					}
					if(KCcount==1){
						RFcount=2
					}
					if(KCcount==2){
						RFcount=4
					}
					
					reachFreqOut=t(rep(NA,RFcount*4))
					RFOnames=c("ALL","ALLreachPcnt","ALLredundPcnt","ALLratio","TopZK","ZKreachPcnt","ZKredundPcnt","ZKratio",
							   "TopZC","ZCreachPcnt","ZCredundPcnt","ZCratio","TopZCK","ZCKreachPcnt","ZCKredundPcnt","ZCKratio")
					colnames(reachFreqOut)=RFOnames[1:length(reachFreqOut)]
					
					#NQS=rbind(NQS,data.frame(Network=count,Similarity=simType,RelSignum=signumvec[i],Signum=signum,MinSize=minSize,MergeBy=merge.by,MergeParam=merge.param,InitialModules=NA,FinalModules=NA,ChangeModules=NA,Members=NA,MeanSpecificity=NA,NGModularity=NA,GModularity=NA,MeanQbase=NA,MeanHomogeneity=NA,MeanPC1VE=NA,MeanSeparation=NA,MeanZKmodgenes=NA,MeanZKbigmodgenes=NA,MeanZCmodgenes=NA,MeanZCbigmodgenes=NA,corModKC=NA,reachFreqOut))
					
					NQS$ChangeModules=(NQS$FinalModules-NQS$InitialModules)/NQS$InitialModules*100
					NQS$ChangeModules[NQS$InitialModules==0]=NA
					NQS=NQS[order(NQS$MinSize,-NQS$Signum),]
					
					#write.table(data.frame(t(NQS)),file=paste(simType,"-",overlapType,"_p",beta,"_",projectname,"_",length(cluster1$order),"_network_statistics_",tstamp1,".csv",sep=""),sep=",",row.names=T,col.names=F)
					write.table(data.frame(t(NQS)),file=paste(simType,"-",overlapType,"_signum",signum[i],"_minSize",minSizevec[j],"_merge_",merge.by,"_",merge.param,"_network_statistics.csv",sep=""),sep=",",row.names=T,col.names=F)
					
				}
				
				setwd(BNrootDir)
				
			}## end of for(j in c(1:length(minSizevec)))
				
		}  ## end of for(i in c(1:length(signumvec)))
		
		# NQS$ChangeModules=(NQS$FinalModules-NQS$InitialModules)/NQS$InitialModules*100
		# NQS$ChangeModules[NQS$InitialModules==0]=NA
		# NQS=NQS[order(NQS$MinSize,-NQS$Signum),]
		# setwd(BNrootDir)		
		# write.table(data.frame(t(NQS)),file=paste(simType,"-",overlapType,"_p",beta,"_",projectname,"_",length(cluster1$order),"_network_statistics_",tstamp1,".csv",sep=""),sep=",",row.names=T,col.names=F)
		# write.table(data.frame(t(NQS)),file=paste(simType,"-",overlapType,"_p",beta,"_",projectname,"_network_statistics.csv",sep=""),sep=",",row.names=T,col.names=F)

	} else { #if(iterate==TRUE)
				  
						print("Finding initial modules...")
						cat("\n")
	
						tstamp1=format(Sys.time(), "%X")
						tstamp1=gsub(" AM","",tstamp1)
						tstamp1=gsub(" PM","",tstamp1)
						tstamp1=gsub(":","-",tstamp1)
		
						if(signumType=="rel"){
		
							relsignum=signum
		
							if(((dim(datExpr1)[2]*(dim(datExpr1)[2]-1))/2)>1e+06){
			
								signum=quantile(sample(simMat,1000000,replace=FALSE),probs=signum)
								collectGarbage()
			
							  } else {
			
								signum=quantile(simMat[upper.tri(simMat)],probs=signum)
								collectGarbage()
			
								}
		
							}
						  
							  rm(simMat)
							  collectGarbage()	
							  
							  cutree1=cutree(cluster1,h=1-signum)
							  
							  keptmodsDF=data.frame(table(cutree1))
							  							  
						if(overlapType=="TO"){
							  
							fileheader=paste(simType,"-",overlapType,"_",TOtype,"_",TOdenom,"_signum",signif(signum,3),"_minSize",minSize,"_mergeBy_",merge.by,merge.param,"_",length(datExpr1[1,]),sep="")
							  
							} else {
							  
							fileheader=paste(simType,"-",overlapType,"_signum",signif(signum,3),"_minSize",minSize,"_mergeBy_",merge.by,merge.param,"_",length(datExpr1[1,]),sep="")
							  
							}
							  
						 dir.create(fileheader)
						 setwd(paste(getwd(),"/",fileheader,sep=""))
							  
## Find initial modules:
							  
						 keptmodscutree=cutree1[cutree1 %in% keptmodsDF$cutree1[keptmodsDF$Freq>=(minSize)]]
						 modules=tapply(as.character(names(keptmodscutree)),factor(keptmodscutree),list)
							  
## Currently the max # of modules is limited by the length of the "colors()" vector (n=657); may need to revisit module naming convention down the road.
							  
						if(length(modules)>0&length(modules)<658){
							  
							  initialModules=length(modules)
							  print(paste("signum = ",signif(signum,3),", minSize = ",minSize,": ",length(modules)," initial modules found",sep=""))
							  cat("\n")
							  
## Calculate module eigengenes and merge similar modules:
							  
							  print("Merging similar modules...")
							  cat("\n")
							  allMEs=mapply(FUN=findMEs,overlap1=mapply(FUN=findoverlap,IDs=modules,MoreArgs=list(expr1=datExpr1),SIMPLIFY=FALSE),MoreArgs=list(expr1=datExpr1),SIMPLIFY=TRUE)
							  MEdf=matrix(unlist(allMEs),ncol=length(modules))
							  MEdf=as.data.frame(MEdf)
							  corMEdf=cor(MEdf)
							  diag(corMEdf)=0
							  collectGarbage()
							  
							  if(length(modules)>1){
							  
							    if(merge.by=="CC"){
							    
							      ## Strategy for merging by clustering coefficient differentials:
							    
							      names(modules)=paste("mod",c(1:length(modules)),sep="")
							    
							      ## Calculate clustering coefficients for initial module seed features:
							      modseed.CC=mapply(FUN=findCCs,mod=modules,MoreArgs=list(adjmat=adjMat),SIMPLIFY=FALSE)
							    
							      ## Merge all possible pairs of module seed features:
							    
							      merged.seeds=map(cross2(modules,modules),unlist)
							      merged.seeds=lapply(merged.seeds,sort)
							      temp2=unique(merged.seeds)
							      temp3=lapply(temp2,unique)
							      merged.seeds=temp2[unlist(lapply(temp2,length))==unlist(lapply(temp3,length))]
							    
							      ## Calculate clustering coefficients for merged module seed features:
							    
							      mergedmodseed.CC=mapply(FUN=findCCs,mod=merged.seeds,MoreArgs=list(adjmat=adjMat),SIMPLIFY=FALSE)
							    
							      ## Name the merged modules:
							    
							      count.mat=matrix(nrow=length(modules),ncol=length(modules),data=-888)
							      count.mat[lower.tri(count.mat)]=c(1:length(mergedmodseed.CC))
							      count.mat=t(count.mat)
							      rownames(count.mat)=paste("mod",c(1:length(modules)),sep="")
							      colnames(count.mat)=paste("mod",c(1:length(modules)),sep="")
							    
							      for(k in 1:length(mergedmodseed.CC)){
							      
							        names(mergedmodseed.CC)[k]=paste(rownames(count.mat)[which(count.mat==k,arr.ind=TRUE)[1]],colnames(count.mat)[which(count.mat==k,arr.ind=TRUE)[2]],sep="|")
							      
							      }
							    
							      ## Calculate the max difference between original CCs and CCs from merged modules:
							    
							      CC.diff.list=vector(mode="list",length=length(mergedmodseed.CC))
							      names(CC.diff.list)=names(mergedmodseed.CC)
							    
							      for(k in 1:length(CC.diff.list)){
							      
							        which.mods=unlist(strsplit(names(mergedmodseed.CC[k]),split="|",fixed=T))
							        CC.diff.list[[k]]=data.frame(CC.init=unlist(modseed.CC[is.element(names(modseed.CC),which.mods)]),CC.merged=NA)
							        rownames(CC.diff.list[[k]])=gsub(paste(which.mods[1],".",sep=""),"",rownames(CC.diff.list[[k]]),fixed=T)
							        rownames(CC.diff.list[[k]])=gsub(paste(which.mods[2],".",sep=""),"",rownames(CC.diff.list[[k]]),fixed=T)
							        CC.diff.list[[k]]=CC.diff.list[[k]][order(rownames(CC.diff.list[[k]])),]
							      
							        if(!all.equal(rownames(CC.diff.list[[k]]),names(mergedmodseed.CC[[k]]))){
							        
							          stop(paste("Feature labels do not align for merged module ",k," when comparing clustering coefficients",sep=""))
							        
							        } else {
							          
							          CC.diff.list[[k]]$CC.merged=mergedmodseed.CC[[k]]
							        
  							      }
							      
	  						      CC.diff.list[[k]]=data.frame(CC.diff.list[[k]],CC.diff=CC.diff.list[[k]][,1]-CC.diff.list[[k]][,2])
							      
							      } ## end of for(k in 1:length(CC.diff.list)){
							    
							      ## Now we will create a module similarity matrix by calculating 1 - the  max CC.diff for every pair of merged modules:
							    
							      CC.diff.max=unlist(lapply(CC.diff.list,function(x){max(x$CC.diff)}))
							      names(CC.diff.max)=names(CC.diff.list)
							    
							      CCdf=matrix(nrow=length(modules),ncol=length(modules),data=-888)
							      CCdf[lower.tri(CCdf)]=1-CC.diff.max
							      CCdf=t(CCdf)
							      CCdf[lower.tri(CCdf)]=1-CC.diff.max
							      diag(CCdf)=0
							      rownames(CCdf)=paste("mod",c(1:length(modules)),sep="")
							      colnames(CCdf)=paste("mod",c(1:length(modules)),sep="")
							    
							      if(export.merge.comp==T){
							      
							        pdf(file=paste(simType,"-",overlapType,"_signum",signif(signum,3),"_minSize",minSizevec[j],"_compare_merge_parameters.pdf",sep=""))
							        plot(corMEdf[upper.tri(corMEdf)],CCdf[upper.tri(CCdf)],xlab="cor(ME)",ylab="1-mean(CC.diff)",main=paste(simType,"-",overlapType,"_signum",signif(signum,3),"_minSize",minSizevec[j],sep=""),xlim=c(-1,1))
							        mtext("Initial module similarities",side=3)
							        abline(v=c(0.85,0.9,0.95))
							        abline(h=c(0.95,0.96,0.97))
							        dev.off()
							      
							      }
							    
							      while(max(CCdf)>=(merge.param)&length(modules)>1){
							      
							        clusterMod=flashClust(as.dist(1-CCdf),method="complete")
							        cutreeMods=cutree(clusterMod,h=1-merge.param)
							        modules=mapply(FUN=mergenodes,cutreeMods1=tapply(c(1:length(cutreeMods)),factor(cutreeMods),list),MoreArgs=list(tempmods1=modules),SIMPLIFY=FALSE)
							      
							        if(length(modules)>1){
							        
							          names(modules)=paste("mod",c(1:length(modules)),sep="")
							      
							          ## Calculate clustering coefficients for initial module seed features:
							          modseed.CC=mapply(FUN=findCCs,mod=modules,MoreArgs=list(adjmat=adjMat),SIMPLIFY=FALSE)
							      
							          ## Merge all possible pairs of module seed features:
							      
							          merged.seeds=map(cross2(modules,modules),unlist)
							          merged.seeds=lapply(merged.seeds,sort)
							          temp2=unique(merged.seeds)
							          temp3=lapply(temp2,unique)
							          merged.seeds=temp2[unlist(lapply(temp2,length))==unlist(lapply(temp3,length))]
							      
							          ## Calculate clustering coefficients for merged module seed features:
							      
							          mergedmodseed.CC=mapply(FUN=findCCs,mod=merged.seeds,MoreArgs=list(adjmat=adjMat),SIMPLIFY=FALSE)
							      
							          ## Name the merged modules:
							      
  							        count.mat=matrix(nrow=length(modules),ncol=length(modules),data=-888)
							          count.mat[lower.tri(count.mat)]=c(1:length(mergedmodseed.CC))
							          count.mat=t(count.mat)
							          rownames(count.mat)=paste("mod",c(1:length(modules)),sep="")
							          colnames(count.mat)=paste("mod",c(1:length(modules)),sep="")
							      
							          for(k in 1:length(mergedmodseed.CC)){
							        
							            names(mergedmodseed.CC)[k]=paste(rownames(count.mat)[which(count.mat==k,arr.ind=TRUE)[1]],colnames(count.mat)[which(count.mat==k,arr.ind=TRUE)[2]],sep="|")
							          
							          }
							      
							          ## Calculate the max difference between original CCs and CCs from merged modules:
							      
							          CC.diff.list=vector(mode="list",length=length(mergedmodseed.CC))
							          names(CC.diff.list)=names(mergedmodseed.CC)
							      
							          for(k in 1:length(CC.diff.list)){
							        
							            which.mods=unlist(strsplit(names(mergedmodseed.CC[k]),split="|",fixed=T))
							            CC.diff.list[[k]]=data.frame(CC.init=unlist(modseed.CC[is.element(names(modseed.CC),which.mods)]),CC.merged=NA)
							            rownames(CC.diff.list[[k]])=gsub(paste(which.mods[1],".",sep=""),"",rownames(CC.diff.list[[k]]),fixed=T)
							            rownames(CC.diff.list[[k]])=gsub(paste(which.mods[2],".",sep=""),"",rownames(CC.diff.list[[k]]),fixed=T)
							            CC.diff.list[[k]]=CC.diff.list[[k]][order(rownames(CC.diff.list[[k]])),]
							        
							            if(!all.equal(rownames(CC.diff.list[[k]]),names(mergedmodseed.CC[[k]]))){
							          
							              stop(paste("Feature labels do not align for merged module ",k," when comparing clustering coefficients",sep=""))
							          
							            } else {
							          
							              CC.diff.list[[k]]$CC.merged=mergedmodseed.CC[[k]]
							          
							            }
							        
							            CC.diff.list[[k]]=data.frame(CC.diff.list[[k]],CC.diff=CC.diff.list[[k]][,1]-CC.diff.list[[k]][,2])
							        
							          } ## end of for(k in 1:length(CC.diff.list)){
							      
							          ## Now we will create a module similarity matrix by calculating 1 - the  max CC.diff for every pair of merged modules:
							      
							          CC.diff.max=unlist(lapply(CC.diff.list,function(x){max(x$CC.diff)}))
							          names(CC.diff.max)=names(CC.diff.list)
							      
							          CCdf=matrix(nrow=length(modules),ncol=length(modules),data=-888)
							          CCdf[lower.tri(CCdf)]=1-CC.diff.max
							          CCdf=t(CCdf)
							          CCdf[lower.tri(CCdf)]=1-CC.diff.max
							          diag(CCdf)=0
							          rownames(CCdf)=paste("mod",c(1:length(modules)),sep="")
							          colnames(CCdf)=paste("mod",c(1:length(modules)),sep="")
							      
							        } ## end of if(length(modules)>1){
							          
							      } ## end of while(max(CCdf)>=(merge.param)){
							    
							    ## Calculate module eigengenes:
							    
							    allMEs=mapply(FUN=findMEs,overlap1=mapply(FUN=findoverlap,IDs=modules,MoreArgs=list(expr1=datExpr1),SIMPLIFY=FALSE),MoreArgs=list(expr1=datExpr1),SIMPLIFY=TRUE)
							    MEdf=matrix(unlist(allMEs),ncol=length(modules))
							    MEdf=as.data.frame(MEdf)
							    
							  } ## end of if(merge.by=="CC"){
							  
							  if(merge.by=="ME"){
							    
							    while(max(corMEdf)>=(merge.param)){
							      
							      clusterMod=flashClust(as.dist(1-corMEdf),method="complete")
							      cutreeMods=cutree(clusterMod,h=1-merge.param)
							      modules=mapply(FUN=mergenodes,cutreeMods1=tapply(c(1:length(cutreeMods)),factor(cutreeMods),list),MoreArgs=list(tempmods1=modules),SIMPLIFY=FALSE)
							      allMEs=mapply(FUN=findMEs,overlap1=mapply(FUN=findoverlap,IDs=modules,MoreArgs=list(expr1=datExpr1),SIMPLIFY=FALSE),MoreArgs=list(expr1=datExpr1),SIMPLIFY=TRUE)
							      MEdf=matrix(unlist(allMEs),ncol=length(modules))
							      MEdf=as.data.frame(MEdf)
							      corMEdf=cor(MEdf)
							      diag(corMEdf)=0
							      collectGarbage()
							      
							      } ## end of while(max(corMEdf)>=(merge.param)&length(modules)>1){
							    
							    } ## end of if(merge.by=="ME"){
							    
							  } ## end of if(length(modules)>1){  
							  
							  finalModules=length(modules)
							  print(paste(length(modules),"modules found"))
							  cat("\n")
							  
## Label modules:
							  
							  modulesize=c()
							  
							  for(m in c(1:dim(MEdf)[2])){
							  
								  modulesize=c(modulesize,length(modules[[m]]))
							  
								}
							  
							  standardcolors=c(standardColors(),setdiff(colors(),standardColors()))
							  standardcolors=standardcolors[1:length(modules)]
							  names(MEdf)=standardcolors[rank(-modulesize,ties.method="first")]
							  names(modules)=names(MEdf)
							  corMEdf=cor(MEdf)
							  diag(corMEdf)=0
							  
## Re-order modules:
							  
							  if(length(modules)>1){
							
								  whichMethod="average"
								  MEcluster=flashClust(as.dist(1-corMEdf),method=whichMethod)
								  MEdf=MEdf[,MEcluster$order]
								  corMEdf=cor(MEdf)
								  diag(corMEdf)=0
								  modules=modules[MEcluster$order]
								  modulesize=modulesize[MEcluster$order]
								  
							  }
							  
## Export network summary:
							  
							  print("Exporting network summary...")
							  cat("\n")
							  
							  MEdfout=data.frame(colnames(expr)[sampleindex],MEdf)
							  colnames(MEdfout)[1]="Sample"
							  
							  write.table(MEdfout,file=paste("Module_eigengenes_",tstamp1,".csv",sep=""),sep=",",row.names=F,col.names=T)
							  
							  if(length(modules)>2){
							  
								pdf(file=paste("ME_network_",tstamp1,".pdf",sep=""),width=21,height=12)
								par(mfrow=c(2,1))
								par(mar=c(.5,5,6,2))
								plot(MEcluster,sub="",xlab="",ylab="1 - cor",main=projectname,cex.main=1.5,cex.lab=1.25,cex=0.8)
								mtext(paste(whichMethod," linkage,"," signum = ",signif(signum,3),", minSize = ",minSize,", merge.by = ",merge.by," > ",merge.param,", ",length(datExpr1[1,])," features",sep=""))
								par(mar = c(1,5,6,2))
								plot.mat(MEdf,nrgcols=500,rlabels=colnames(expr)[sampleindex],clabels=colnames(MEdf),rcols=colorvec[samplegroups])
								dev.off()
								  
								corMEdfOut=corMEdf
								diag(corMEdfOut)=1
								pdf(file=paste("ME_correlations_",tstamp1,".pdf",sep=""),width=16,height=16)
								colors <- c("#A50F15","#DE2D26","#FB6A4A","#FCAE91","#FEE5D9","white",
								            "#EFF3FF","#BDD7E7","#6BAED6","#3182BD","#08519C")
								plotcorr(corMEdfOut,col=colors[5*corMEdfOut+6],main=projectname)
								mtext(paste(whichMethod," linkage,"," signum = ",signif(signum,3),", minSize = ",minSize,", merge.by = ",merge.by," > ",merge.param,", ",length(datExpr1[1,])," features",sep=""))
							    dev.off()
							  
								} else {
							  
								print("Less than three modules; no network summary exported.")
								cat("\n")
							  
							  }

## Calculate kME table:
							  
							  print("Building kME table...")
							  cat("\n")
							  
							  collectGarbage()
							  kMEtable=data.frame(expr[var.all>0,geneinfo])
							  
							  modseeds=data.frame(unlist(modules),rep(names(modules),modulesize))
							
							  if(class(expr[,1])=="integer"|class(expr[,1])=="numeric"){
								
								  modseeds=modseeds[order(as.numeric(as.character(modseeds[,1]))),]
								
								  } else {
								
								  modseeds=modseeds[order(as.character(modseeds[,1])),]
								
								  }
							
							  ModSeeds=rep(NA,length(dim(datExprall)[2]))
							  ModSeeds[is.element(colnames(datExprall),modseeds[,1])]=as.character(modseeds[is.element(modseeds[,1],colnames(datExprall)),2])
							
							  count1=0
							  
							  for(k in c(1:dim(MEdf)[2])){
							  
							  kME1=cor(MEdf[,k],datExprall)
							  kME1pval=corPvalueStudent(kME1,length(sampleindex))
## Convention: p-values for perfect positive or negative correlations will be set to 1e-300 so as not to break downstream code.  Very, very rare event.
							  kME1pval[is.na(kME1pval)]=1e-300
							  kMEtable=data.frame(kMEtable,t(kME1),t(kME1pval))
							  colnames(kMEtable)[length(geneinfo)+(2*count1)+1]=paste("kME",colnames(MEdf)[k],sep="")
							  colnames(kMEtable)[length(geneinfo)+(2*count1)+2]=paste("kME",colnames(MEdf)[k],".pval",sep="")
							  count1=count1+1
							  
							  }
							  
							  collectGarbage()
							  
## Find p-value for determining module specificity.  Will use Bonferroni p-value as one (more stringent) measure of specificity.
## Will also use FDR <.05 as another (less stringent) measure of module specificity.
## Both measures will be exported to kME table; but BC measure will be used for module specificity calculations that are present
## in module statistics and network statistics outputs.
							
							  pvalseq=grep(".pval",colnames(kMEtable))
							  corseq=pvalseq-1
							  
							  allPvalues=c()
							
							  if(length(pvalseq)==1){
								
								allPvalues=kMEtable[,pvalseq]
								
							    } else {
								
								for(n in c(1:length(pvalseq))){
							  
								  allPvalues=c(allPvalues,kMEtable[,pvalseq][,n])
							  
								}
								
							  }
							  
							  allQvalues=qvalue(allPvalues)
							  specPvalCutQ=max(allQvalues$pvalues[allQvalues$qvalues<.05])
							  
							  specPvalCutB=.05/(length(kMEtable[,1])*length(pvalseq))
							  
							  collectGarbage()
							  
## Find probe sets that are most strongly associated with each module (thresholded by p-value):
							  
							  print("Finding top modules for all genes...")
							  cat("\n")
							  
							  pvalseq=grep(".pval",colnames(kMEtable))
							  corseq=pvalseq-1
							  
							  if(length(corseq)==1){
							
							  MMspecScoresBpos=rep(NA,length(kMEtable[,1]))
							  MMspecScoresQpos=rep(NA,length(kMEtable[,1]))
							  MMspecScoresBneg=rep(NA,length(kMEtable[,1]))
							  MMspecScoresQneg=rep(NA,length(kMEtable[,1]))
								  
							  } else {
								  
							  MMspecScoresBpos=CalcMMspecScores(kmetable1=kMEtable,seq1a=corseq,seq2a=pvalseq,pvalcut1a=specPvalCutB,direction="pos")
							  MMspecScoresQpos=CalcMMspecScores(kmetable1=kMEtable,seq1a=corseq,seq2a=pvalseq,pvalcut1a=specPvalCutQ,direction="pos")
							  MMspecScoresBneg=CalcMMspecScores(kmetable1=kMEtable,seq1a=corseq,seq2a=pvalseq,pvalcut1a=specPvalCutB,direction="neg")
							  MMspecScoresQneg=CalcMMspecScores(kmetable1=kMEtable,seq1a=corseq,seq2a=pvalseq,pvalcut1a=specPvalCutQ,direction="neg")
								  
							  }
							
## Calculate mean (log2) expression levels and percentiles:
							
							  FeatureMean=apply(log(datExprall,2),2,mean,na.rm=T)
							  FeaturePercentile=(FeatureMean/max(FeatureMean,na.rm=T))*100
							
							  kMEtable=data.frame(kMEtable[,1:length(geneinfo)],ModSeed=ModSeeds,MeanExpr=FeatureMean,MeanExprPcnt=FeaturePercentile,TopModPosBC=MMspecScoresBpos,TopModPosFDR=MMspecScoresQpos,TopModNegBC=MMspecScoresBneg,TopModNegFDR=MMspecScoresQneg,kMEtable[,(length(geneinfo)+1):dim(kMEtable)[2]])
							  colnames(kMEtable)[1:length(geneinfo)]=colnames(expr)[1:length(geneinfo)]
                              colnames(kMEtable)[length(geneinfo)+2]="MeanExpr"
							  colnames(kMEtable)[length(geneinfo)+3]="MeanExprPercentile"
							  colnames(kMEtable)[length(geneinfo)+4]=paste("TopModPosBC_",signif(specPvalCutB,3),sep="")
							  colnames(kMEtable)[length(geneinfo)+5]=paste("TopModPosFDR_",signif(specPvalCutQ,3),sep="")
							  colnames(kMEtable)[length(geneinfo)+6]=paste("TopModNegBC_",signif(specPvalCutB,3),sep="")
							  colnames(kMEtable)[length(geneinfo)+7]=paste("TopModNegFDR_",signif(specPvalCutQ,3),sep="")
							
							  if(writeKME==TRUE){
							
								write.table(kMEtable,file=paste("kME_table_",tstamp1,".csv",sep=""),sep=",",row.names=F,col.names=T)
								
							  }
							
							  collectGarbage()
							  
## Export module heat maps, ME barplots, and lineplots:
							  
                              if(writeModSnap==TRUE){
                              
                                print("Exporting heat maps & barplots & lineplots...")
                                cat("\n")
							  
                                pdf(file=paste("Modules_",tstamp1,".pdf",sep=""),width=10,height=7.5)
                                nf=layout(matrix(c(1,3,2,3), 2, 2, byrow=TRUE), respect=FALSE)
							  
                                for(p in c(1:dim(MEdf)[2])){
							  
                                    restmod=is.element(colnames(datExprall),modules[[p]])
                                    datmod=datExprall[,restmod]
                                    kmerank=rank(-cor(MEdf[,p],datmod))
                                    datmod=datmod[,order(kmerank)]
                                    par(mar=c(0,2,5,2))
                                    plot.mat(t(scale(datmod)),title=names(MEdf)[p],cex.main=1.5)
                                    par(mar=c(10,2,1,2))
                                    barplot1=barplot(MEdf[,p],col=gsub("ME","",names(MEdf)[p]),xaxs="i",yaxs="i")
                                    axis(side=1,at=barplot1[,1],labels=F,tick=T)
                                    box()
                                    mtext(rownames(datmod),1,at=barplot1[,1],adj=1.3,cex=0.7,las=3,col=colorvec[samplegroups])
                                    par(mar=c(10,4,5,2))
                                    nogenes=15
                                    restmod2=is.element(colnames(kMEtable),paste("kME",names(modules[p]),sep=""))
                                    overlapTOP=is.element(colnames(datExprall),kMEtable[,1][order(-kMEtable[,restmod2])][1:nogenes])
                                    datmod2=datExprall[,overlapTOP]
                                    ordering=order(rank(-cor(MEdf[,p],datmod2)))
							   
                                    datMat1=scale(datmod2[,1:nogenes][,ordering])
                                    Gene1=as.character(expr$Gene[is.element(expr[,1],colnames(datMat1))])[ordering]
                                    multiples=grep("/",Gene1)
                                  
                                    if(length(multiples)>0){
                                  
                                        for(m in c(1:length(multiples))){
                                  
                                            Gene1[multiples[m]]=paste(strsplit(Gene1[multiples[m]],split=" /")[[1]][1],"*",sep="")
                                  
                                            }
                                  
                                        }
                                  
                                    if(max(datMat1)>10){
                                      subtr=max(datMat1)/2
                                    } else {
                                      subtr=5
                                    }
                                    
                                    matplot(data.frame(c(1:dim(datmod2)[1])),data.frame(datMat1),main="Top genes by kME", ylab="Expression level (z-score)", xlab="Sample", type="l", lwd=2, cex.main=1.5, cex.lab=1.25, axes=F,col=standardColors()[1:nogenes],lty=1,ylim=c(min(datMat1)-subtr,max(datMat1)+.1))
                                    legend("bottom",Gene1,cex=1,lty=1,col=standardColors()[1:nogenes],lwd=2,ncol=3)
                                    axis(side=2,tick=T)
                                    box()
                                  
                                    } ## end of for(p in c(1:dim(MEdf)[2]))
                                
                                  dev.off()
                                  
                                  collectGarbage()
                                  
                              } ## end of if(writeModSnap==TRUE){
							  
## Export module quality statistics:
							  
							  print("Calculating module statistics...")
							  cat("\n")	
							  
## Calculate expanded module sizes below kme p-value cutoff (denominators for module specificity calculation):
							  
							bigmodules=unlist(TopModIDs(datkme1=kMEtable,geneinfo1=c(1:(length(geneinfo)+7)),pvalcut1=specPvalCutB,justcount=TRUE))
							  
## Calculate and assemble module quality statistics:	
							  
							bigmoduleIDs=TopModIDs(datkme1=kMEtable,geneinfo1=c(1:(length(geneinfo)+7)),pvalcut1=specPvalCutB,justcount=FALSE)
							names(bigmoduleIDs)=names(modules)
							
							MMspecOut=data.frame(table(MMspecScoresBpos))
							noSpecMods=setdiff(names(modules),MMspecOut[,1])
							
							if(length(noSpecMods)>0){
								
								noSpecModsDF=data.frame(MMspecScoresBpos=noSpecMods,Freq=0)
								MMspecOut=rbind(MMspecOut,noSpecModsDF)
								MMspecOut=MMspecOut[order(as.character(MMspecOut[,1])),]
								
							}
							
							MMspecOut2=MMspecOut[order(match(MMspecOut[,1],names(modules))),]
							
							HScores=mapply(FUN=HScore,overlap1=mapply(FUN=findoverlap,IDs=modules,MoreArgs=list(expr1=datExpr1),SIMPLIFY=FALSE),
										   overlap2=mapply(FUN=findoverlap,IDs=names(modules),MoreArgs=list(expr1=MEdf),SIMPLIFY=FALSE),
										   MoreArgs=list(expr1=datExpr1),SIMPLIFY=TRUE)
							
							signCorMEdf=(corMEdf+1)/2
							diag(signCorMEdf)=0
							
## Note: in a signed network with beta=2, mean Separation should be ~0.5 with bimodal distribution (< 0.5 = positive ME correlations; > 0.5 = negative ME correlations).
							SScores=1-apply(signCorMEdf,1,SSfx)
							
							PC1VE=unlist(mapply(FUN=MEPC1,overlap1=mapply(FUN=findoverlap,IDs=modules,MoreArgs=list(expr1=datExpr1),SIMPLIFY=FALSE),MoreArgs=list(expr1=datExpr1),SIMPLIFY=TRUE))
							
							adjME=((cor(MEdf)+1)/2)^2
							diag(adjME)=0
							
							if(length(modules)>2){
								
								FNCME=fundamentalNetworkConcepts(adjME)
								
								if(length(unique(FNCME$Connectivity))>2&length(unique(FNCME$ClusterCoef))>2){
								
									C.k=cor.test(scale(FNCME$Connectivity),scale(FNCME$ClusterCoef),method="s")$estimate
									
								} else {
									
									C.k=NA
									
								}
								
							} else {
								
								FNCME=data.frame(Connectivity=rep(NA,length(modules)),ClusterCoef=rep(NA,length(modules)))
								C.k=NA
								
							}
							
							overlap1=mapply(FUN=findoverlap,IDs=modules,MoreArgs=list(expr1=datExprall),SIMPLIFY=TRUE)
							
							if(length(grep("Z.K",colnames(expr)))==1){
								
								meanZKmodgenes=c()
								
								for(r in c(1:length(modules))){
									
									meanZKmodgenes=c(meanZKmodgenes,mean(expr$Z.K[overlap1[,r]],na.rm=T))
									
								}
								
							} else {
								
								meanZKmodgenes=rep(NA,length(modules))
								
							}
							
							if(length(grep("Z.C",colnames(expr)))==1){
								
								meanZCmodgenes=c()
								
								for(r in c(1:length(modules))){
									
									meanZCmodgenes=c(meanZCmodgenes,mean(expr$Z.C[overlap1[,r]],na.rm=T))
									
								}
								
							} else {
								
								meanZCmodgenes=rep(NA,length(modules))
								
							}
							
							modValues=mapply(FUN=GetValues,overlap1=mapply(FUN=findoverlap,IDs=modules,MoreArgs=list(expr1=datExprall),SIMPLIFY=FALSE),
											 MoreArgs=list(expr1=datExprall),SIMPLIFY=FALSE)
							
							if(length(modules)>1){
								
								MeanExpr=sapply(modValues,mean,na.rm=T)
								Z.MeanExpr=scale(MeanExpr)
								SEExpr=sapply(modValues,stderr1)
								MEVout=mapply(FUN=MEV,overlap1=mapply(FUN=findoverlap,IDs=modules,MoreArgs=list(expr1=datExprall),SIMPLIFY=FALSE),
											  MoreArgs=list(expr1=datExprall),SIMPLIFY=TRUE)
								
							} else {
								
								modValues=unlist(modValues)
								MeanExpr=mean(modValues,na.rm=T)
								Z.MeanExpr=scale(MeanExpr)
								SEExpr=stderr1(modValues)
								overlap1=is.element(colnames(datExprall),modules)
								MEVout=rbind(mean(apply(log(datExprall[,overlap1],2),2,var)),stderr1(apply(log(datExprall[,overlap1],2),2,var)))
								
							}
							
							MeanExprVar=unlist(MEVout[1,])
							SEMeanExprVar=unlist(MEVout[2,])
							Z.MeanExprVar=scale(MeanExprVar)
							
							if(calcBigModStat==TRUE){
							
								bigmodValues=mapply(FUN=GetValues,overlap1=mapply(FUN=findoverlap,IDs=bigmoduleIDs,MoreArgs=list(expr1=datExprall),SIMPLIFY=FALSE),
												MoreArgs=list(expr1=datExprall),SIMPLIFY=FALSE)
							
								if(length(bigmoduleIDs)>1){
								
									MeanExprBCMod=sapply(bigmodValues,mean,na.rm=T)
									Z.MeanExprBCMod=scale(MeanExprBCMod)
									MEVoutBCMod=mapply(FUN=MEV,overlap1=mapply(FUN=findoverlap,IDs=bigmoduleIDs,MoreArgs=list(expr1=datExprall),SIMPLIFY=FALSE),
												   MoreArgs=list(expr1=datExprall),SIMPLIFY=TRUE)
								
									if(sum(is.na(MeanExprBCMod))>0){
									
										SEExprBCMod=rep(NA,length(MeanExprBCMod))
										SEExprBCMod[!is.na(MeanExprBCMod)]=sapply(bigmodValues[!is.na(MeanExprBCMod)],stderr1)
									
									} else {
									
										SEExprBCMod=sapply(bigmodValues,stderr1)
									
									}
								
								} else {
								
									bigmodValues=unlist(bigmodValues)
									MeanExprBCMod=mean(bigmodValues,na.rm=T)
									Z.MeanExprBCMod=scale(MeanExprBCMod)
									SEExprBCMod=stderr1(bigmodValues)
									overlap1=is.element(colnames(datExprall),bigmoduleIDs)
									MEVoutBCMod=rbind(mean(apply(log(datExprall[,overlap1],2),2,var)),stderr1(apply(log(datExprall[,overlap1],2),2,var)))
								
								}
							
								MeanExprVarBCMod=unlist(MEVoutBCMod[1,])
								SEMeanExprVarBCMod=unlist(MEVoutBCMod[2,])				
								Z.MeanExprVarBCMod=scale(MeanExprVarBCMod)
								
							} else {
								
								MeanExprBCMod=NA
								Z.MeanExprBCMod=NA
								SEExprBCMod=NA
								MeanExprVarBCMod=NA
								Z.MeanExprVarBCMod=NA
								SEMeanExprVarBCMod=NA
								
							}
								
							MEdfRank=apply(abs(MEdf),2,sort,decreasing=TRUE)
							MEmax=MEdfRank[1,]
							PosOrNeg=apply(MEdf,2,function(x){abs(max(x))>abs(min(x))})
							MEstep1=rep(-88888,length(MEmax))
							MEstep1[PosOrNeg]=MEmax[PosOrNeg]-MEdfRank[2,PosOrNeg]
							MEdfRank=apply(MEdf,2,sort,decreasing=TRUE)
							MEstep1[!PosOrNeg]=MEdfRank[dim(MEdfRank)[1],!PosOrNeg]-MEdfRank[dim(MEdfRank)[1]-1,!PosOrNeg]
							MEstep1=abs(MEstep1)
							
							if(calcSW==TRUE){
								
								QM=data.frame(Module=names(modules),InitialSize=modulesize,ExpandedSize=bigmodules,Unique=MMspecOut2$Freq,Homogeneity=HScores,ZHomogeneity=scale(HScores),Separation=SScores,ZSeparation=scale(SScores),Silhouette=SWScores,ZSilhouette=scale(SWScores),PC1VE=PC1VE)
								QM=data.frame(QM[,1:4],Specificity=QM$Unique/QM$ExpandedSize*100,QM[,5:11],Z.Kmod=scale(FNCME$Connectivity),Z.Cmod=scale(FNCME$ClusterCoef),MeanZKgenemodgenes=meanZKmodgenes,MeanZCmodgenes=meanZCmodgenes,MeanExpr=MeanExpr,ZMeanExpr=Z.MeanExpr,SEExpr=SEExpr,MeanExprVar=MeanExprVar,ZMeanExprVar=Z.MeanExprVar,SEMeanExprVar=SEMeanExprVar,MeanExprBCMod=MeanExprBCMod,ZMeanExprBCMod=Z.MeanExprBCMod,SEExprBCMod=SEExprBCMod,MeanExprVarBCMod=MeanExprVarBCMod,ZMeanExprVarBCMod=Z.MeanExprVarBCMod,SEMeanExprVarBCMod=SEMeanExprVarBCMod,MEmax=MEmax,MEstep1=MEstep1)
								colnames(QM)[grep("Unique",colnames(QM))]=paste("Unique_",signif(specPvalCutB,3),sep="")
								colnames(QM)[grep("ExpandedSize",colnames(QM))]=paste("ExpandedSize_",signif(specPvalCutB,3),sep="")
								
							} else {
								
								QM=data.frame(Module=names(modules),InitialSize=modulesize,ExpandedSize=bigmodules,Unique=MMspecOut2$Freq,Homogeneity=HScores,ZHomogeneity=scale(HScores),Separation=SScores,ZSeparation=scale(SScores),PC1VE=PC1VE)
								QM=data.frame(QM[,1:4],Specificity=QM$Unique/QM$ExpandedSize*100,QM[,5:9],Z.Kmod=scale(FNCME$Connectivity),Z.Cmod=scale(FNCME$ClusterCoef),MeanZKmodgenes=meanZKmodgenes,MeanZCmodgenes=meanZCmodgenes,MeanExpr=MeanExpr,ZMeanExpr=Z.MeanExpr,SEExpr=SEExpr,MeanExprVar=MeanExprVar,ZMeanExprVar=Z.MeanExprVar,SEMeanExprVar=SEMeanExprVar,MeanExprBCMod=MeanExprBCMod,ZMeanExprBCMod=Z.MeanExprBCMod,SEExprBCMod=SEExprBCMod,MeanExprVarBCMod=MeanExprVarBCMod,ZMeanExprVarBCMod=Z.MeanExprVarBCMod,SEMeanExprVarBCMod=SEMeanExprVarBCMod,MEmax=MEmax,MEstep1=MEstep1)
								colnames(QM)[grep("Unique",colnames(QM))]=paste("Unique_",signif(specPvalCutB,3),sep="")
								colnames(QM)[grep("ExpandedSize",colnames(QM))]=paste("ExpandedSize_",signif(specPvalCutB,3),sep="")
								
							}
							
							write.table(QM,file=paste("Module_statistics_",tstamp1,".csv",sep=""),sep=",",row.names=F,col.names=T)
							
							pdf(file=paste("Module_statistics_",tstamp1,".pdf",sep=""),width=12,height=12)
							
							par(mar=c(10,5,4,2))
							barplot(QM[,2],main="Initial module size",col=as.character(QM$Module),ylab="Members",cex.main=2,cex.lab=1.5,names=as.character(QM$Module),las=3)
							
							par(mar=c(10,5,4,2))
							barplot(QM[,3],main=paste("Expanded module (P < ",gsub("ExpandedSize_","",colnames(QM)[3]),") size",sep=""),col=as.character(QM$Module),ylab="Members",cex.main=2,cex.lab=1.5,names=as.character(QM$Module),las=3)
							
							par(mar=c(5,5,4,2))
							plot(QM[,2],QM[,3],xlab="Initial size",ylab=paste(gsub("_"," (P < ",colnames(QM)[3]),")",sep=""),main="Module size",cex.main=2,cex.lab=1.5,pch=19,col=as.character(QM$Module))
							if(dim(QM)[1]>2&length(unique(QM[,2]))>1){
								abline(lm(QM[,3]~QM[,2]),col="red",lwd=2)
							}
							
							par(mar=c(5,5,4,2))
							plot(QM[,3],QM[,4],xlab=paste(gsub("_"," (P < ",colnames(QM)[3]),")",sep=""),ylab=paste(gsub("_"," (P < ",colnames(QM)[4]),")",sep=""),main="Module size and specificity",cex.main=2,cex.lab=1.5,pch=19,col=as.character(QM$Module))
							if(dim(QM)[1]>2){
								abline(lm(QM[,4]~QM[,3]),col="red",lwd=2)
							}
							abline(0,1)
							
							par(mar=c(10,5,4,2))
							barplot(QM$Specificity,main="Module specificity",col=as.character(QM$Module),ylab="Percentage",cex.main=2,cex.lab=1.5,names=as.character(QM$Module),las=3)
							
							par(mar=c(5,5,4,2))
							plot(QM$Homogeneity,QM$PC1VE,xlab="Homogeneity",ylab="PC1VE",main="Module purity",cex.main=2,cex.lab=1.5,pch=19,col=as.character(QM$Module))
							if(dim(QM)[1]>2&!is.na(lm(QM$PC1VE~QM$Homogeneity)[[1]][2])){
								abline(lm(QM$PC1VE~QM$Homogeneity),col="red",lwd=2)
							}
							
							par(mar=c(5,5,4,2))
							plot(QM$InitialSize,QM$PC1VE,xlab="Initial size",ylab="PC1VE",main="Module purity vs size",cex.main=2,cex.lab=1.5,pch=19,col=as.character(QM$Module))
							
							if(sum(!is.na(QM$Z.Kmod))!=0){
								
								par(mar=c(5,5,4,2))
								plot(QM$Z.Kmod,QM$Z.Cmod,xlab="Module Z.K",ylab="Module Z.C",main="Module cor(K,C)",cex.main=2,cex.lab=1.5,pch=19,col=as.character(QM$Module))
								mtext(paste("rho = ",signif(cor.test(QM$Z.Kmod,QM$Z.Cmod,method="s")$estimate,2),sep=""))
								if(dim(QM)[1]>2){
									abline(lm(QM$Z.Cmod~QM$Z.Kmod),col="red",lwd=2)
								}
								
							}
							
							par(mar=c(10,5,4,2))
							sdME=sd(QM$MeanExpr)
							meanME=mean(QM$MeanExpr)
							barplot(QM$MeanExpr,main="Module expression levels",col=as.character(QM[,1]),ylab="Mean (log2) expression",cex.main=2,cex.lab=1.5,names=as.character(QM$Module),las=3,ylim=c(min(0,min(QM$MeanExpr,na.rm=T)-abs(min(QM$MeanExpr,na.rm=T)*0.1)),max(QM$MeanExpr,na.rm=T)+max(QM$MeanExpr,na.rm=T)*0.1))
							addErrorBars(means=QM$MeanExpr,errors=QM$SEExpr,two.side=TRUE)
							abline(h=meanME,col="black",lwd=2)
							abline(h=meanME+(2*sdME),col="red",lwd=2)
							abline(h=meanME-(2*sdME),col="red",lwd=2)
							
							par(mar=c(10,5,4,2))
							sdMEV=sd(QM$MeanExprVar)
							meanMEV=mean(QM$MeanExprVar)
							barplot(QM$MeanExprVar,main="Module expression variance",col=as.character(QM[,1]),ylab="Mean (log2) variance",cex.main=2,cex.lab=1.5,names=as.character(QM$Module),las=3,ylim=c(min(0,min(QM$MeanExprVar,na.rm=T)-abs(min(QM$MeanExprVar,na.rm=T)*0.1)),max(QM$MeanExprVar,na.rm=T)+max(QM$MeanExprVar,na.rm=T)*0.1))
							addErrorBars(means=QM$MeanExprVar,errors=QM$SEMeanExprVar,two.side=TRUE)
							abline(h=meanMEV,col="black",lwd=2)
							abline(h=meanMEV+(2*sdMEV),col="red",lwd=2)
							abline(h=meanMEV-(2*sdMEV),col="red",lwd=2)
							
							if(!all(is.nan(QM$Specificity))&calcBigModStat==TRUE){
								par(mar=c(10,5,4,2))
								sdME=sd(QM$MeanExprBCMod,na.rm=T)
								meanME=mean(QM$MeanExprBCMod,na.rm=T)
								barplot(QM$MeanExprBCMod,main=paste("Expanded module (P < ",gsub("ExpandedSize_","",colnames(QM)[3]),") expression levels",sep=""),col=as.character(QM[,1]),ylab="Mean (log2) expression",cex.main=2,cex.lab=1.5,names=as.character(QM$Module),las=3,ylim=c(min(0,min(QM$MeanExprBCMod,na.rm=T)-abs(min(QM$MeanExprBCMod,na.rm=T)*0.1)),max(QM$MeanExprBCMod,na.rm=T)+max(QM$MeanExprBCMod,na.rm=T)*0.1))
								addErrorBars(means=QM$MeanExprBCMod,errors=QM$SEExprBCMod,two.side=TRUE)
								abline(h=meanME,col="black",lwd=2)
								abline(h=meanME+(2*sdME),col="red",lwd=2)
								abline(h=meanME-(2*sdME),col="red",lwd=2)
							} else {
								par(mar=c(10,5,4,2))
								barplot(rep(0,length(QM$MeanExprBCMod)),main=paste("Expanded module (P < ",gsub("ExpandedSize_","",colnames(QM)[3]),") expression levels",sep=""),col=as.character(QM[,1]),ylab="Mean (log2) expression",cex.main=2,cex.lab=1.5,names=as.character(QM$Module),las=3)
							}
							
							if(!all(is.nan(QM$Specificity))&calcBigModStat==TRUE){
								par(mar=c(10,5,4,2))
								sdMEV=sd(QM$MeanExprVarBCMod,na.rm=T)
								meanMEV=mean(QM$MeanExprVarBCMod,na.rm=T)
								barplot(QM$MeanExprVarBCMod,main=paste("Expanded module (P < ",gsub("ExpandedSize_","",colnames(QM)[3]),") variance",sep=""),col=as.character(QM[,1]),ylab="Mean (log2) variance",cex.main=2,cex.lab=1.5,names=as.character(QM$Module),las=3,ylim=c(min(0,min(QM$MeanExprVarBCMod,na.rm=T)-abs(min(QM$MeanExprVarBCMod,na.rm=T)*0.1)),max(QM$MeanExprVarBCMod,na.rm=T)+max(QM$MeanExprVarBCMod,na.rm=T)*0.1))
								addErrorBars(means=QM$MeanExprVarBCMod,errors=QM$SEMeanExprVarBCMod,two.side=TRUE)
								abline(h=meanMEV,col="black",lwd=2)
								abline(h=meanMEV+(2*sdMEV),col="red",lwd=2)
								abline(h=meanMEV-(2*sdMEV),col="red",lwd=2)
							} else {
								par(mar=c(10,5,4,2))
								barplot(rep(0,length(QM$MeanExprVarBCMod)),main=paste("Expanded module (P < ",gsub("ExpandedSize_","",colnames(QM)[3]),") variance",sep=""),col=as.character(QM[,1]),ylab="Mean (log2) variance",cex.main=2,cex.lab=1.5,names=as.character(QM$Module),las=3)
							}
														
							par(mar=c(10,5,4,2))
							barplot(QM$MEmax,main="|Module eigengene maximum|",col=as.character(QM$Module),ylab="Highest value",cex.main=2,cex.lab=1.5,names=as.character(QM$Module),las=3)
							
							par(mar=c(10,5,4,2))
							barplot(QM$MEstep1,main="Module eigengene first step",col=as.character(QM$Module),ylab="Difference between top two values",cex.main=2,cex.lab=1.5,names=as.character(QM$Module),las=3)
							
							dev.off()
														  
							write.table(QM,file=paste("Module_statistics_",tstamp1,".csv",sep=""),sep=",",row.names=F,col.names=T)
							  
## Export network statistics:				
							  
							  print("Calculating network statistics...")
							  cat("\n")
							  
## Calculate % probe sets that are present in multiple modules and reachFreq:
							  
							  pvalseq=grep(".pval",colnames(kMEtable))
							  corseq=pvalseq-1
							  
							  LogiMods=matrix(nrow=dim(kMEtable)[1],ncol=length(pvalseq),data=FALSE)
							  
							  for(p in c(1:length(pvalseq))){
							  
								LogiMods[,p][kMEtable[,corseq][,p]>0&kMEtable[,pvalseq][,p]<specPvalCutB]=TRUE
							  
								}
							  
							  No.Modules=rowSums(LogiMods,na.rm=T)
							  No.ModulesDF=data.frame(table(No.Modules))
							  names(No.Modules)=colnames(datExprall)
							  
							  reachFreqOut=reachFreq(expr1=expr,subset1=subset,ZNCcut1=ZNCcut,modcount1=No.Modules,tstamp1=tstamp1)
							  
##	Calculate modularity and mean Z.K / Z.C:
							  
							  ModOut=Modularity(datexpr1=datExpr1,modules1=modules)
							  
							  if(length(grep("Z.K",colnames(expr)))==1){
							  
								meanZKmodgenes=mean(expr$Z.K[subset][is.element(expr[subset,1],unique(unlist(modules)))],na.rm=T)
								meanZKbigmodgenes=mean(expr$Z.K[subset][is.element(expr[subset,1],unique(unlist(bigmoduleIDs)))],na.rm=T)
							  
							  
								} else {
							  
								meanZKmodgenes=NA
								meanZKbigmodgenes=NA
							  
								}
							  
							  if(length(grep("Z.C",colnames(expr)))==1){
							  
								meanZCmodgenes=mean(expr$Z.C[subset][is.element(expr[subset,1],unique(unlist(modules)))],na.rm=T)
								meanZCbigmodgenes=mean(expr$Z.C[subset][is.element(expr[subset,1],unique(unlist(bigmoduleIDs)))],na.rm=T)
							  
								} else {
							  
								meanZCmodgenes=NA
								meanZCbigmodgenes=NA
							  
								}
							  
							  collectGarbage()
							  
							  NQS=data.frame(Network=1,Similarity=simType,RelSignum=relsignum,Signum=signum,MinSize=minSize,MergeBy=merge.by,MergeParam=merge.param,,InitialModules=initialModules,FinalModules=finalModules,ChangeModules=0,Members=length(unique(unlist(modules))),MeanSpecificity=mean(QM$Specificity),NGModularity=ModOut$NGMod,GModularity=ModOut$GMod,MeanQbase=ModOut$MeanQbase,MeanHomogeneity=mean(HScores),MeanPC1VE=mean(PC1VE),MeanSeparation=mean(SScores),MeanZKmodgenes=meanZKmodgenes,MeanZKbigmodgenes=meanZKbigmodgenes,MeanZCmodgenes=meanZCmodgenes,MeanZCbigmodgenes=meanZCbigmodgenes,corModKC=C.k,reachFreqOut)	
							  NQS$ChangeModules=(NQS$FinalModules-NQS$InitialModules)/NQS$InitialModules*100
							  							  
							  } else { ## end of if(length(modules)>0)
							  
								NQS=data.frame(Network=1,Similarity=simType,RelSignum=relsignum,Signum=signum,MinSize=minSize,MergeBy=merge.by,MergeParam=merge.param,InitialModules=0)
										
								}
										
							  #setwd(BNrootDir)
																				
							  write.table(data.frame(t(NQS)),file=paste("Network_statistics_",tstamp1,".csv",sep=""),sep=",",row.names=T,col.names=F)
										
							  collectGarbage()
		
	} ## end of if(iterate==TRUE
				   
	timestamp()
	
} ## end of function

FindModules(projectname,
            expr,
            geneinfo,
            sampleindex,
            samplegroups=NULL,
            subset=NULL,
            simMat=NULL,
            saveSimMat=FALSE,
            simType,
            overlapType,
            TOtype,
            TOdenom,
            beta=1,
            MIestimator="mi.mm",
            MIdisc="equalfreq",
            signumType=c("abs","rel"),
            iterate=FALSE,
            signumvec=c(.999,.99,.98),
            minsizevec=c(8,10,12),
            signum,
            minSize,
            merge.by=c("ME","CC"),
            merge.param,
            export.merge.comp=T,
            ZNCcut=2,
            calcSW=FALSE,
            loadTree=FALSE,
            writeKME=FALSE,
            calcBigModStat=FALSE,
            writeModSnap=TRUE
            )


