

convert2Matrix = function()
{
    arguments = commandArgs(trailingOnly = TRUE)
    inputFile = as.character(arguments[1])
    metadata = as.character(arguments[2])
    maf_cutoff = as.numeric(arguments[3])
    outputDir = as.character(arguments[4])
    prefix = as.character(arguments[5])
    
    require(tictoc)
    require(data.table)
    require(doParallel)
    
    tic('time to check the input data')
    if(!file.exists(inputFile))
        stop(inputFile, ' no such file or directory')
    if(!file.exists(metadata))
        stop(metadata, ' no such file or directory')
    if(!dir.exists(outputDir))
    {
        system(sprintf('mkdir -p %s', outputDir))
        message(outputDir, ' has been created.')
    }
    if(maf_cutoff>=0 & maf_cutoff<=0.4)
        message('The relatedness values obtained after filtering loci with a MAF below ',maf_cutoff,' will be considered.')
    else
        stop('The MAF cut-off should be one of these values [0, 0.1, 0.2, 0.3, 0.4]')
    toc()
    
    relatednessData = readRDS(inputFile)
    meta = fread(metadata)
    
    if(maf_cutoff==0)
    {
        names(relatednessData)[3]='ofInterest'
        relatednessData$ofInterest = as.numeric(as.character(relatednessData$ofInterest))
        # dataOfInterest = relatednessData[which(relatednessData$ofInterest>=r_cutoff),]
        # rm(relatednessData)
    }
    else if(maf_cutoff==0.1)
    {
        names(relatednessData)[4]='ofInterest'
        relatednessData$ofInterest = as.numeric(as.character(relatednessData$ofInterest))
        # dataOfInterest = relatednessData[which(relatednessData$ofInterest>=r_cutoff),]
        # rm(relatednessData)
    }
    else if(maf_cutoff==0.2)
    {
        names(relatednessData)[5]='ofInterest'
        relatednessData$ofInterest = as.numeric(as.character(relatednessData$ofInterest))
        # dataOfInterest = relatednessData[which(relatednessData$ofInterest>=r_cutoff),]
        # rm(relatednessData)
    }
    else if(maf_cutoff==0.3)
    {
        names(relatednessData)[6]='ofInterest'
        relatednessData$ofInterest = as.numeric(as.character(relatednessData$ofInterest))
        # dataOfInterest = relatednessData[which(relatednessData$ofInterest>=r_cutoff),]
        # rm(relatednessData)
    }
    else if(maf_cutoff==0.4)
    {
        names(relatednessData)[7]='ofInterest'
        relatednessData$ofInterest = as.numeric(as.character(relatednessData$ofInterest))
        # dataOfInterest = relatednessData[which(relatednessData$ofInterest>=r_cutoff),]
        # rm(relatednessData)
    }
    
    
    #surLigne = unique(relatednessData$iid1)
    modifier = function(x){gsub('-','.',x)}
    meta$SampleID=lapply(meta$SampleID,modifier)
    # relatednessMatrix = matrix(NA,nrow = length(unique(relatednessData$iid1)), ncol = length(unique(relatednessData$iid2)))
    # rownames(relatednessMatrix) = unique(relatednessData$iid1)   #meta$SampleID
    # colnames(relatednessMatrix) = unique(relatednessData$iid2)   #meta$SampleID
    relatednessMatrix = matrix(NA,nrow = length(unique(meta$SampleID)), ncol = length(unique(meta$SampleID)))
    rownames(relatednessMatrix) = unique(meta$SampleID)   #meta$SampleID
    colnames(relatednessMatrix) = unique(meta$SampleID) 
    surLigne = unique(relatednessData$iid1)     #meta$SampleID
    # cl = makeCluster(10)
    # registerDoParallel(cl)
    for(i in 1:length(surLigne)){
        # print(paste0('i=',i))
        ligne = surLigne[i]
        l = match(ligne,rownames(relatednessMatrix))
        target = relatednessData[which(relatednessData$iid1==ligne),]
        t = unique(target$iid2)
        for(j in 1:length(t)){
            tt = target[which(target$iid2==t[j]),]
            k = match(t[j],colnames(relatednessMatrix)) 
            if(nrow(tt) == 0)
                relatednessMatrix[l,k] = 0
            else if(nrow(tt)==1)
                relatednessMatrix[l,k] = round(as.numeric(tt$ofInterest), digits = 2)
            else if(nrow(tt)>1 & length(unique(tt$iid1))==1 & length(unique(tt$iid2))==1)
                relatednessMatrix[l,k] = round(as.numeric(tt$ofInterest[1]), digits = 2)
            else if(nrow(tt)>1 & length(unique(tt$iid1))>1 | length(unique(tt$iid2))>1)
                relatednessMatrix[l,k] = mean(round(as.numeric(tt$ofInterest), digits = 2), na.rm = TRUE)
        }
    }
    # stopCluster(cl)
    saveRDS(relatednessMatrix, file = paste0(outputDir,'/',prefix,'.RDS'))
}


convert2Matrix()