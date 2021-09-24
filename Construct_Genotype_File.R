#' @title 
#' Construct genotype file
#'
#' @description 
#' this function filter out SNPs and isolates that do not satisfy the given MAF and missingness cut-offs 
#' @param pathToVCFfile full path to the filtered VCF file
#' @param pathToOutDir full path to the output files directory
#' @param metadataFile full path to the file with sample metadata file
#' @param sweepRegions full path to the file with genomic locations of the regions to be discarded 
#' @param type the mixed genotypes recoding method 
#' @return output files will be stored in output directory
#' @examples
#' Rscript Construct_Genotype_File.R Results/MyGenotype.txt Data/MyDelgeme4.vcf.gz Data/SampleMetadata.txt Data/Pop_Order.txt

constructGenotypeFile = function()
{
    arguments = commandArgs(trailingOnly = TRUE)
    pathToVCFfile = as.character(arguments[1])
    pathToOutDir = as.character(arguments[2])  
    metadataFile = as.character(arguments[3])
    sweepRegions = as.character(arguments[4])
    type = as.character(arguments[5])  
    
    require(data.table)
    require(tictoc)
    require(statip)
    
    tic("Time to recode data")
    metadata = fread(metadataFile)
    if(file.exists(sweepRegions))
        selectiveRegions = fread(sweepRegions)
    genotype = paste0(pathToOutDir,'/Genotype.txt')
    expression = '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'
    system(sprintf("bcftools query -f'%s' %s > %s", expression, pathToVCFfile, genotype))
    allelicDepth = paste0(pathToOutDir,'/AllelicDepth.txt')
    expression = '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n'
    system(sprintf("bcftools query -f'%s' %s > %s", expression, pathToVCFfile, allelicDepth))
    sampleList = paste0(pathToOutDir,'/SampleList.txt')
    system(sprintf("bcftools query -l %s > %s", pathToVCFfile, sampleList))
    gt = fread(genotype, nThread = 15, header = FALSE)
    ad = fread(allelicDepth, nThread = 15, header = FALSE)
    sampleIDs = fread(sampleList, header = FALSE)
    if(type=='raw' | type=='Raw' | type=='RAW')
    {
        outDir = paste0(pathToOutDir,'/Raw_Genotypes')
        system(sprintf("mkdir -p %s", outDir))
        f4c = subset(gt, select=c(1:4));gt = subset(gt, select=-c(1:4));gt = as.matrix(gt)
        recodedData = matrix(NA, nrow = dim(gt)[1], ncol = dim(gt)[2])
        for(i in 1:dim(gt)[1])
        {
            m=which(gt[i,]=='0/0');recodedData[i,m]=0
            m=which(gt[i,]=='0/1');recodedData[i,m]=2
            m=which(gt[i,]=='1/1');recodedData[i,m]=1
            #m=which(gt[i,]=='./.');recodedData[i,m]=NA
        }
        recodedData = as.data.frame(recodedData)
        names(recodedData) = sampleIDs$V1
        pops = unique(metadata$Location)
        #print(paste0('verdict=',sweepRegions))
        if(file.exists(sweepRegions))
        {
            rtd = NULL
            for(j in 1:nrow(selectiveRegions))
            {
                rtd = c(rtd, which(f4c$V1==selectiveRegions$Chrom[j] & (f4c$V2>=selectiveRegions$Start[j] & f4c$V2<=selectiveRegions$End[j])))
                # rtd = c(rtd, which(f4c$Chrom==selectiveRegions$Chrom[j] & (f4c$Pos>=selectiveRegions$Start[j] & f4c$Pos<=selectiveRegions$End[j])))
                
            }
            rtd = unique(rtd)
            f4c = f4c[-rtd,]
            recodedData = recodedData[-rtd,]
            
        }
        for(i in 1:length(pops))
        {
            target = metadata[which(metadata$Location==pops[i]),]$SampleID
            m = match(target, names(recodedData))
            X=subset(recodedData, select=m)
            chroms=f4c$V1
            pos=f4c$V2
            samps=sampleIDs$V1
            L = list(X, chroms, pos, samps)
            names(L)=c('X','chroms','pos','samps')
            saveRDS(L, paste0(outDir,'/',pops[i],'_Raw_Genotype.RDS'))
        }
        system(sprintf("rm -f %s %s %s", genotype, allelicDepth, sampleList))
    }
    else if(type=='minor_Karim' | type=='Minor_Karim' | type=='MINOR_Karim')
    {
        outDir = paste0(pathToOutDir,'/Minor_Karim_Genotypes')
        system(sprintf("mkdir -p %s", outDir))
        f4c = subset(gt, select=c(1:4));gt = subset(gt, select=-c(1:4));gt = as.matrix(gt)
        ad = subset(ad, select=-c(1:4));ad = as.matrix(ad)
        
        recodedData = matrix(NA, nrow = dim(gt)[1], ncol = dim(gt)[2])
        for(i in 1:dim(gt)[1])
        {
            m=which(gt[i,]=='0/0');recodedData[i,m]=0; m=which(gt[i,]=='1/1');recodedData[i,m]=1
            m=which(gt[i,]=='./.');recodedData[i,m]=NA
            m=which(gt[i,]=='0/1')
            if(length(m)==0)
                next
            for(j in 1:length(m))
            {
                if(is.character(ad[i,m[j]]))
                {
                    r = as.numeric(unlist(strsplit(ad[i,m[j]],','))[1]);a = as.numeric(unlist(strsplit(ad[i,m[j]],','))[2]) 
                    if(r==0 & a==0)
                    {
                        rr=length(which(recodedData[i,]==0)); aa=length(which(recodedData[i,]==1))
                        if(rr<aa)
                            recodedData[i,m[j]] = 0
                        else if(rr>aa)
                            recodedData[i,m[j]] = 1
                        else
                            recodedData[i,m[j]] = rbern(1, rr/(rr+aa))
                    }
                    else if(r!=0 & a!=0)
                    {
                        if(r<a)
                            recodedData[i,m[j]] = 0
                        else if(r>a)
                            recodedData[i,m[j]] = 1
                        else
                        {
                            rr=length(which(recodedData[i,]==0)); aa=length(which(recodedData[i,]==1))
                            if(rr<aa)
                                recodedData[i,m[j]] = 0
                            else if(rr>aa)
                                recodedData[i,m[j]] = 1
                            else
                                recodedData[i,m[j]] = rbern(1, rr/(rr+aa))
                        }
                    }
                    else if(r==0 | a==0)
                    {
                        if(r<a)
                            recodedData[i,m[j]] = 0
                        else if(r>a)
                            recodedData[i,m[j]] = 1
                        else
                        {
                            rr=length(which(recodedData[i,]==0)); aa=length(which(recodedData[i,]==1))
                            if(rr<aa)
                                recodedData[i,m[j]] = 0
                            else if(rr>aa)
                                recodedData[i,m[j]] = 1
                            else
                                recodedData[i,m[j]] = rbern(1, rr/(rr+aa))
                        }
                    }
                }
                
            }
        }
        
        recodedData = as.data.frame(recodedData)
        names(recodedData) = sampleIDs$V1
        pops = unique(metadata$Location)
        if(sweepRegions != 'none')
        {
            rtd = NULL
            for(j in 1:nrow(selectiveRegions))
                rtd = c(rtd, which(f4c$V1==selectiveRegions$Chrom[j] & (f4c$V2>=selectiveRegions$Start[j] & f4c$V2<=selectiveRegions$End[j])))
            rtd = unique(rtd)
            f4c = f4c[-rtd,]
            recodedData = recodedData[-rtd,]
        }
        for(i in 1:length(pops))
        {
            target = metadata[which(metadata$Location==pops[i]),]$SampleID
            m = match(target, names(recodedData))
            X=subset(recodedData, select=m)
            chroms=f4c$V1; pos=f4c$V2; samps=sampleIDs$V1
            L = list(X, chroms, pos, samps)
            names(L)=c('X','chroms','pos','samps')
            saveRDS(L, paste0(outDir,'/',pops[i],'_Minor_Karim_Genotype.RDS'))
        }
        system(sprintf("rm -f %s %s %s", genotype, allelicDepth, sampleList))
    }
    else if(type=='minor_David' | type=='Minor_David' | type=='MINOR_David')
    {
        outDir = paste0(pathToOutDir,'/Minor_David_Genotypes')
        system(sprintf("mkdir -p %s", outDir))
        f4c = subset(gt, select=c(1:4));gt = subset(gt, select=-c(1:4));gt = as.matrix(gt)
        ad = subset(ad, select=-c(1:4));ad = as.matrix(ad)
        
        recodedData = matrix(NA, nrow = dim(gt)[1], ncol = dim(gt)[2])
        for(i in 1:dim(gt)[1])
        {
            m=which(gt[i,]=='0/0');recodedData[i,m]=0
            m=which(gt[i,]=='1/1');recodedData[i,m]=1
            m=which(gt[i,]=='./.');recodedData[i,m]=NA
            m=which(gt[i,]=='0/1')
            for(j in 1:length(m))
            {
                rr=length(which(gt[i,]=='0/0')); aa=length(which(gt[i,]=='1/1'))
                if(rr<aa)
                    recodedData[i,m[j]] = 0
                else if(rr>aa)
                    recodedData[i,m[j]] = 1
                else
                    recodedData[i,m[j]] = rbern(1, rr/(rr+aa))
            }
        }
        
        recodedData = as.data.frame(recodedData)
        names(recodedData) = sampleIDs$V1
        pops = unique(metadata$Location)
        if(sweepRegions != 'none')
        {
            rtd = NULL
            for(j in 1:nrow(selectiveRegions))
                rtd = c(rtd, which(f4c$V1==selectiveRegions$Chrom[j] & (f4c$V2>=selectiveRegions$Start[j] & f4c$V2<=selectiveRegions$End[j])))
            rtd = unique(rtd)
            f4c = f4c[-rtd,]
            recodedData = recodedData[-rtd,]
        }
        for(i in 1:length(pops))
        {
            target = metadata[which(metadata$Location==pops[i]),]$SampleID
            m = match(target, names(recodedData))
            X=subset(recodedData, select=m)
            chroms=f4c$V1; pos=f4c$V2; samps=sampleIDs$V1
            L = list(X, chroms, pos, samps)
            names(L)=c('X','chroms','pos','samps')
            saveRDS(L, paste0(outDir,'/',pops[i],'_Minor_David_Genotype.RDS'))
        }
        system(sprintf("rm -f %s %s %s", genotype, allelicDepth, sampleList))
    }
    else if(type=='missing' | type=='Missing' | type=='MISSING')
    {
        outDir = paste0(pathToOutDir,'/Missing_Genotypes')
        system(sprintf("mkdir -p %s", outDir))
        f4c = subset(gt, select=c(1:4));gt = subset(gt, select=-c(1:4));gt = as.matrix(gt)
        recodedData = matrix(NA, nrow = dim(gt)[1], ncol = dim(gt)[2])
        for(i in 1:dim(gt)[1])
        {
            m=which(gt[i,]=='0/0');recodedData[i,m]=0
            m=which(gt[i,]=='0/1');recodedData[i,m]=NA
            m=which(gt[i,]=='1/1');recodedData[i,m]=1
            m=which(gt[i,]=='./.');recodedData[i,m]=NA
        }
        
        recodedData = as.data.frame(recodedData)
        names(recodedData) = sampleIDs$V1
        pops = unique(metadata$Location)
        if(sweepRegions != 'none')
        {
            rtd = NULL
            for(j in 1:nrow(selectiveRegions))
                rtd = c(rtd, which(f4c$V1==selectiveRegions$Chrom[j] & (f4c$V2>=selectiveRegions$Start[j] & f4c$V2<=selectiveRegions$End[j])))
            rtd = unique(rtd)
            f4c = f4c[-rtd,]
            recodedData = recodedData[-rtd,]
        }
        for(i in 1:length(pops))
        {
            target = metadata[which(metadata$Location==pops[i]),]$SampleID
            m = match(target, names(recodedData))
            X=subset(recodedData, select=m)
            chroms=f4c$V1; pos=f4c$V2; samps=sampleIDs$V1
            L = list(X, chroms, pos, samps)
            names(L)=c('X','chroms','pos','samps')
            saveRDS(L, paste0(outDir,'/',pops[i],'_Missing_Genotype.RDS'))
        }
        system(sprintf("rm -f %s %s %s", genotype, allelicDepth, sampleList))
    }
    
    toc()
}







constructGenotypeFile()

