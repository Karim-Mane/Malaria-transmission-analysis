#' @title 
#' filter SNPs and isolates
#'
#' @description 
#' this function filter out SNPs and isolates that do not satisfy the given MAF and missingness cut-offs 
#' @param inputFile the input VCF file
# @param sampleIds file with the isolate IDs
#' @param mafCutOff numeric value. SNPs with a MAF below this value will be removed
#' @param missingnessCutOff numeric value. SNPs and isolates with a missingness greater than this value will be discarded
#' @param outputDir path to the folder where to store the output files
#' @param numThreads the number of threads to use during the processing
#' @return output files will be stored in output directory
#' @examples
#' Rscript Filter_data.R Data/MyDelgeme4.vcf.gz 0.01 20 Results/ 3
filtration = function()
{
    arguments = commandArgs(trailingOnly=TRUE)
    inputFile = as.character(arguments[1])
    # sampleIds = as.character(arguments[2])
    mafCutOff = as.numeric(arguments[3])
    missingnessCutOff = as.numeric(arguments[4])
    outputDir = as.character(arguments[5])
    numThreads = as.numeric(arguments[6])
    
    require(data.table)
    require(foreach)
    require(doParallel)
    require(tictoc)
    
    tic("Filtration time")
    if(!file.exists(inputFile))
        stop(inputFile," no such file or directory!")
    # if(!file.exists(sampleIds))
    #     stop(sampleIds," no such file or directory!")
    if(!dir.exists(outputDir))
        stop(outputDir, " no such file or directory!")
    if(mafCutOff<0 | mafCutOff>1)
        stop("incorrect MAF cut-off! value should be between 0 and 1")
    if(missingnessCutOff<0 | missingnessCutOff>100)
        stop("incorrect missingness cut-off! value should be between 0 and 100")

    
    # sampleIds = paste0(outputDir,'sampleIds.txt')
    # system(sprintf("bcftools query -l %s > %s", inputFile, sampleIds))
    sampleIds = paste0(pathToOutDir,'/SampleList.txt')
    system(sprintf("bcftools query -l %s > %s", inputFile, sampleIds))
    genotypeFile = paste0(pathToOutDir,'/Genotype.txt')
    expression = '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'
    system(sprintf("bcftools query -f'%s' %s > %s", expression, inputFile, genotypeFile))
    
    # genotypeFile = inputFile
    genotypes = fread(genotypeFile, header = FALSE)
    listOfIsolates = fread(sampleIds, header = FALSE)
    
    message("filtering SNPs with missingness > ", missingnessCutOff, "%")
    f4c = subset(genotypes, select=c(1:4))
    names(f4c) = c('CHROM','POS',"REF",'ALT')
    genotypes = subset(genotypes, select=-c(1:4))
    # if(ncol(genotypes)==(nrow(listOfIsolates)+1))
    #     genotypes = subset(genotypes, select=-ncol(genotypes))
    names(genotypes) = listOfIsolates$V1
    genotypes = as.matrix(genotypes)
    cl = makeCluster(numThreads)
    registerDoParallel(cl)
    percentMissingGenotype = foreach(i=1:nrow(genotypes), .combine = "c") %dopar% {
        m = length(which(genotypes[i,]=='./.' | genotypes[i,]=='.|.'))/ncol(genotypes)
        m
    }
    p = hist(percentMissingGenotype, 100, main="distribution of SNPs missingness", plot=FALSE)
    pdf(paste0(outputDir,"/SNPs_missingness.pdf"))
    plot(p)
    dev.off()
    stopCluster(cl)
    m = which(percentMissingGenotype>(missingnessCutOff/100))  #(missingnessCutOff*ncol(genotypes))/100
    if(length(m)>0)
    {
        message(length(m)," SNPs have missingness > ", missingnessCutOff/100)
        filteredSnps = genotypes[m,]
        filteredF4c = f4c[m,]
        filtered = cbind(filteredF4c, as.data.frame(filteredSnps))
        fwrite(filtered, paste0(outputDir,'/Filtered_SNPs_OnMissingness.txt'), col.names = TRUE, row.names = FALSE, quote = FALSE, nThread = numThreads)
        genotypes = genotypes[-m,]
        f4c = f4c[-m,]
    }
    else
        message(length(m)," SNPs have missingness > ", missingnessCutOff,"%")
    
    
    message("filtering isolates with missingness > ", missingnessCutOff, "%")
    cl = makeCluster(numThreads)
    registerDoParallel(cl)
    percentMissingGenotype = foreach(i=1:ncol(genotypes), .combine = "c") %dopar% {
        m = length(which(genotypes[,i]=='./.' | genotypes[,i]=='.|.'))/nrow(genotypes)
        m
    }
    stopCluster(cl)
    p = hist(percentMissingGenotype, 100, main="distribution of isolates missingness", plot=FALSE)
    pdf(paste0(outputDir,"/Isolates_missingness.pdf"))
    plot(p)
    dev.off()
    m = which(percentMissingGenotype>(missingnessCutOff/100))
    if(length(m)>0)
    {
        message(length(m)," isolates have missingness > ", missingnessCutOff, "%")
        filteredIsolates = listOfIsolates$V1[m]
        fwrite(filteredIsolates, paste0(outputDir,'/Filtered_Isolates_OnMissingness.txt'), col.names = FALSE, row.names = FALSE, quote = FALSE, nThread = numThreads, sep = "\t")
        genotypes = as.data.frame(genotypes, stringAsFactor=FALSE)
        genotypes = subset(genotypes, select=-c(m))
    }
    else
        message(length(m)," isolates have missingness > ", missingnessCutOff,"%")
    
    message("filtering SNPs with MAF < ", mafCutOff)
    outMaf = compute_MAF(genotypes, numThreads)
    message("range of maf: ",range(outMaf))
    p = hist(outMaf, 100, main="distribution of SNPs MAF", plot=FALSE)
    pdf(paste0(outputDir,"/SNPs_MAF.pdf"))
    plot(p)
    dev.off()
    m = which(outMaf<mafCutOff)
    if(length(m)>0)
    {
        message(length(m)," SNPs have MAF < ", mafCutOff)
        filteredSnps = genotypes[m,]
        filteredF4c = f4c[m,]
        filtered = cbind(filteredF4c, as.data.frame(filteredSnps))
        fwrite(filtered, paste0(outputDir,'/Filtered_SNPs_OnMAF.txt'), col.names = TRUE, row.names = FALSE, quote = FALSE, nThread = numThreads, sep = "\t")
        genotypes = genotypes[-m,]
        f4c = f4c[-m,]
        
    }
    
    genotypes = cbind(f4c, as.data.frame(genotypes, stringAsFactor=FALSE))
    fwrite(genotypes, paste0(outputDir,'/MyGenotype.txt'), col.names = TRUE, row.names = FALSE, quote = FALSE, nThread = numThreads, sep = "\t")
    
    toc()
}


#' @title 
#' compute SNPs MAF
#'
#' @description 
#' This function calculates the MAF for each SNPs in the input matrix
#' @param genotypes Matrix of genotype data with SNPs in row and isolates in column. Data should be in raw VCF format i.e. 0/0, 1/1, 0/1, ./.
#' @param numThreads number of threads to use during calculation
#' @return numeric vector
#' @examples
#' maf = compute_MAF(myMatrix, 10)
compute_MAF = function(genotypes, numThreads)
{
    require(foreach)
    require(doParallel)
    
    genotypes = as.matrix(genotypes)
    cl = makeCluster(numThreads)
    registerDoParallel(cl)
    maf = foreach(i = 1:dim(genotypes)[1], .combine = "c") %dopar%
    {
        ref = length(which(genotypes[i,]=='0/0' | genotypes[i,]=='0|0'))
        alt = length(which(genotypes[i,]=='1/1' | genotypes[i,]=='1|1'))
        if(ref<alt)
            maf = ref/(ref+alt)
        else
            maf = alt/(ref+alt)
        maf    
    }
    stopCluster(cl)
    return(maf)
}















filtration()



