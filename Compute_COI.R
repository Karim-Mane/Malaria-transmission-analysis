#' @title 
#' Compute Fws
#'
#' @description 
#' this function filter out SNPs and isolates that do not satisfy the given MAF and missingness cut-offs 
#' @param inputFile the genotypes file generated from the Filter_Data.R script
#' @param vcfFile the input VCF file
#' @param metaFile the file with sample metadata 
#' @param popOrder the file with population order, 1 per line 
#' @return output files will be stored in output directory
#' @examples
#' Rscript Compute_COI.R Results/MyGenotype.txt Data/MyDelgeme4.vcf.gz Data/SampleMetadata.txt Data/Pop_Order.txt
calculateCOI = function()
{
    arguments = commandArgs(trailingOnly=TRUE)
    inputFile = as.character(arguments[1])
    vcfFile = as.character(arguments[2])
    metaFile = as.character(arguments[3])
    popOrder = as.character(arguments[4])
    
    # devtools::install_github("bahlolab/moimix")
    # BiocManager::install("SeqArray")
    require(data.table)
    require(tictoc)
    require(SeqArray)
    require(moimix)
    require(ggplot2)
    

    tic("Time to calculate COI")
    message("\nfiltering out the SNPs and isolates that did not pass the QC")
    if(!file.exists(inputFile))
        stop(inputFile," no such file or directory!")
    if(!file.exists(vcfFile))
        stop(vcfFile," no such file or directory!")
    if(!file.exists(metaFile))
        stop(metaFile," no such file or directory!")
    
    outputDir = dirname(inputFile)
    header = paste0(outputDir,'/','Header.txt')
    body = paste0(outputDir,'/','Body.txt')
    correctRows = paste0(outputDir,'/','Good_snps.txt')
    filteredVcf = paste0(outputDir,'/','Filtered_genotypes.vcf')
    failedIsolates = paste0(outputDir,'/Filtered_Isolates_OnMissingness.txt')
    passQCgenotypes = paste0(outputDir,'/','Passed_QC_genotypes.vcf')
        
    system(sprintf("bcftools view -h %s > %s", vcfFile, header))
    system(sprintf("bcftools view -H %s > %s", vcfFile, body))
    system(sprintf("awk -F'\t' 'NR==FNR{c[$1$2]++;next};c[$1$2] > 0' %s %s > %s",inputFile,body,correctRows))
    system(sprintf("cat %s %s > %s", header, correctRows, filteredVcf))
    if(file.exists(failedIsolates))
    {
        inFile = fread(inputFile)
        gIds = names(inFile)[5:ncol(inFile)]
        goodIds = paste0(outputDir,'/','Good_ids.txt')
        fwrite(gIds, goodIds, col.names = FALSE, row.names = FALSE, quote = FALSE)
        system(sprintf("bcftools query -S %s %s > %s", goodIds, filteredVcf, passQCgenotypes))
        system(sprintf("bgzip %s", passQCgenotypes))
        system(sprintf("rm -f %s %s %s %s",header, body, correctRows, goodIds))
    }
    else
    {
        system(sprintf("mv %s %s", filteredVcf, passQCgenotypes))
        system(sprintf("bgzip %s", passQCgenotypes))
        system(sprintf("rm -f %s %s %s",header, body, correctRows))
    }
    
    
    message("\ncalculating the COI")
    passQCgenotypes = paste0(outputDir,'/','Passed_QC_genotypes.vcf.gz')
    gdsFile = paste0(outputDir,'/','Mydata.gds')
    seqVCF2GDS(passQCgenotypes, gdsFile)  ##Converting the VCF to GDS format (the VCF file should be tabixed, gzipped)
    my_vcf = seqOpen(gdsFile)     #read the GDS file into R and check that its output matches the orginal VCF file
    seqSummary(my_vcf)
    sample.id = seqGetData(my_vcf, "sample.id")     # save sample IDs
    coords = getCoordinates(my_vcf)      # get genomic coordinates of all variants
    seqSetFilter(my_vcf, variant.id = coords$variant.id[coords$chromosome != "Pf3D7_API_v3"])     # filter variant.ids not on apicoplast
    seqSetFilter(my_vcf, variant.id = coords$variant.id[coords$chromosome != "Pf_M76611"])      # filter variant.ids not on mitochondrion
    fws_overall = getFws(my_vcf)      #estimate the MOI
    fws_overall = cbind(as.character(sample.id), as.numeric(fws_overall))
    fws_overall = as.data.frame(fws_overall)
        
    message("\nploting COI as a boxplot")
    popOrder = fread(popOrder, header = FALSE)
    metadata = fread(metaFile)
    m = match(fws_overall$V1, metadata$SampleID)
    fws_overall$V3=metadata$Location[m]
    fws_overall$V3 = factor(fws_overall$V3, levels = popOrder$V1)
    fws_overall$V2 = as.numeric(as.character(fws_overall$V2))
    names(fws_overall) = c('sampleID', 'Fws', 'Location')
    p = ggplot(fws_overall, aes(x=Location, y=Fws, fill=Location)) + 
        geom_boxplot() +
        scale_fill_manual(values=c('yellow', 'gray','pink','red','orange','green','brown','blue','seagreen','cyan','black')) +  #,'magenta','gold'
        theme(axis.text.x = element_text(angle = 60, hjust = 1))+
        theme( axis.line = element_line(colour = "black", 
                                        size = 0.5, linetype = "solid"),
               legend.position = "none") +
        theme(
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()
        ) +
        theme(axis.text.x = element_text(face="bold", color="black"), #, size=14
              axis.text.y = element_text(face="bold", color="black")) #, size=14
    pdf(paste0(outputDir,'/','Fws_COI.pdf'))
    print(p)
    dev.off()
    coiFile = paste0(outputDir,'/','Fws_COI.RDS')
    saveRDS(fws_overall, coiFile)
    
    message("determining the proportion of monoclonal isolates by location")
    fws_overall$status="polyclonal"
    fws_overall[which(fws_overall$Fws>0.95),]$status='monoclonal'
    message("Percentage of monoclonal isolates: ",nrow(fws_overall[which(fws_overall$Fws>0.95),])/nrow(myFws))
    t = table(fws_overall$Location, fws_overall$status)
    t = as.matrix(t)
    pmono = round((t[,1]/rowSums(t))*100, 3); ppoly = round((t[,2]/rowSums(t))*100, 3)
    t = cbind(t, pmono, ppoly)
    colnames(t) = c('#monoclonal','#polyclonal','%monoclonal','%polyclonal')
    saveRDS(t, paste0(outputDir,'/','Populations_COI_Status.RDS'))
    sites = rep(rownames(t),2)
    ps = c(t[,3],t[,4])
    st = c(rep('monoclonal',nrow(t)), rep('polyclonal',nrow(t)))
    data = data.frame(cbind(sites,ps, st))
    names(data)=c('location','proportion','state')
    p=ggplot(data, aes(fill=state, y=proportion, x=location)) + 
        geom_bar(position="dodge", stat="identity")
    pdf(paste0(outputDir,'/','Populations_COI_Status.pdf'))
    print(p)
    dev.off()

    toc()
}














calculateCOI()