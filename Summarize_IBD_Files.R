

summarizeFiles = function()
{
    require(fst)
    require(doParallel)
    require(data.table)
    
    args = commandArgs(trailingOnly = TRUE)
    inputDir = as.character(args[1])
    outputDir = as.character(args[2])
    
    if(!dir.exists(inputDir)){
        stop(inputDir, " No such file or directory")
    }
    if(!dir.exists(outputDir)){
        system(sprintf("mkdir -p %s", outputDir))
    }
    directoriesList = list.dirs(inputDir, full.names = TRUE, recursive = TRUE)
    directoriesList = directoriesList[2:length(directoriesList)]
    # print(directoriesList)
    cl = makeCluster(10)
    registerDoParallel(cl)
    foreach(i = 1:length(directoriesList), .packages = c("data.table","fst")) %dopar%{
        prefix = unlist(strsplit(directoriesList[i],'/'))[length(unlist(strsplit(directoriesList[i],'/')))]
        maf0 = list.files(directoriesList[i], pattern = "*_0p0000.RDS", full.names = TRUE)
        if(length(maf0)==0){
          message("There was no file obtained for maf cutoff=0% for ", prefix)
          # next
        }else if(length(maf0)==1){
            summary_maf0 = readRDS(maf0[1])
            saveRDS(summary_maf0, paste0(outputDir,'/',prefix,'_all_maf>0.RDS'))
            # write_fst(summary_maf0, paste0(outputDir,'/',prefix,'_all_maf>0.fst'))
            summary_maf0 = as.data.table(summary_maf0)
            names(summary_maf0) = c('iid1','iid2','k','rMaf>0')
            summary_maf0 = subset(summary_maf0, select = -3)
            setkeyv(summary_maf0, c('iid1','iid2'))
        }else if(length(maf0)>1){
          # print(prefix)
          # print(length(maf0))
            summary_maf0 = readRDS(maf0[1])
            for(j in 2:length(maf0)){
                summary_maf0 = rbind(summary_maf0, readRDS(maf0[j]))
            }
            saveRDS(summary_maf0, paste0(outputDir,'/',prefix,'_all_maf>0.RDS'))
            # write_fst(summary_maf0, paste0(outputDir,'/',prefix,'_all_maf>0.fst'))
            summary_maf0 = as.data.table(summary_maf0)
            names(summary_maf0) = c('iid1','iid2','k','rMaf>0')
            summary_maf0 = subset(summary_maf0, select = -3)
            setkeyv(summary_maf0, c('iid1','iid2'))
        }
        
        maf10 = list.files(directoriesList[i], pattern = "*_0p1000.RDS", full.names = TRUE)
        if(length(maf10)==0){
            message("There was no file obtained for maf cutoff=10% for ", prefix)
            # next
        }else if(length(maf10)==1){
            summary_maf10 = readRDS(maf10[1])
            saveRDS(summary_maf10, paste0(outputDir,'/',prefix,'_all_maf>10.RDS'))
            # write_fst(summary_maf10, paste0(outputDir,'/',prefix,'_all_maf>10.fst'))
            summary_maf10 = as.data.table(summary_maf10)
            names(summary_maf10) = c('iid1','iid2','k','rMaf>10')
            summary_maf10 = subset(summary_maf10, select = -3)
            setkeyv(summary_maf10, c('iid1','iid2'))
            res1 = summary_maf0[summary_maf10, nomatch=0]
        }else if(length(maf10)>1){
            summary_maf10 = readRDS(maf10[1])
            for(j in 2:length(maf10))
                summary_maf10 = rbind(summary_maf10, readRDS(maf10[j]))
            saveRDS(summary_maf10, paste0(outputDir,'/',prefix,'_all_maf>10.RDS'))
            # write_fst(summary_maf10, paste0(outputDir,'/',prefix,'_all_maf>10.fst'))
            summary_maf10 = as.data.table(summary_maf10)
            names(summary_maf10) = c('iid1','iid2','k','rMaf>10')
            summary_maf10 = subset(summary_maf10, select = -3)
            setkeyv(summary_maf10, c('iid1','iid2'))
            res1 = summary_maf0[summary_maf10, nomatch=0]
        }
        
        maf20 = list.files(directoriesList[i], pattern = "*_0p2000.RDS", full.names = TRUE)
        if(length(maf20)==0){
            message("There was no file obtained for maf cutoff=20% for ", prefix)
            # next
        }else if(length(maf20)==1){
            summary_maf20 = readRDS(maf20[1])
            saveRDS(summary_maf20, paste0(outputDir,'/',prefix,'_all_maf>20.RDS'))
            # write_fst(summary_maf20, paste0(outputDir,'/',prefix,'_all_maf>20.fst'))
            summary_maf20 = as.data.table(summary_maf20)
            names(summary_maf20) = c('iid1','iid2','k','rMaf>20')
            summary_maf20 = subset(summary_maf20, select = -3)
            setkeyv(summary_maf20, c('iid1','iid2'))
            setkeyv(res1, c('iid1','iid2'))
            res2 = res1[summary_maf20, nomatch=0]
        }else if(length(maf20)>1){
            summary_maf20 = readRDS(maf20[1])
            for(j in 2:length(maf20))
                summary_maf20 = rbind(summary_maf20, readRDS(maf20[j]))
            saveRDS(summary_maf20, paste0(outputDir,'/',prefix,'_all_maf>20.RDS'))
            # write_fst(summary_maf20, paste0(outputDir,'/',prefix,'_all_maf>20.fst'))
            summary_maf20 = as.data.table(summary_maf20)
            names(summary_maf20) = c('iid1','iid2','k','rMaf>20')
            summary_maf20 = subset(summary_maf20, select = -3)
            setkeyv(summary_maf20, c('iid1','iid2'))
            setkeyv(res1, c('iid1','iid2'))
            res2 = res1[summary_maf20, nomatch=0]
        }
        
        maf30 = list.files(directoriesList[i], pattern = "*_0p3000.RDS", full.names = TRUE)
        if(length(maf30)==0){
            message("There was no file obtained for maf cutoff=30% for ", prefix)
            # next
        }else if(length(maf30)==1){
            summary_maf30 = readRDS(maf30[1])
            saveRDS(summary_maf30, paste0(outputDir,'/',prefix,'_all_maf>30.RDS'))
            # write_fst(summary_maf30, paste0(outputDir,'/',prefix,'_all_maf>30.fst'))
            summary_maf30 = as.data.table(summary_maf30)
            names(summary_maf30) = c('iid1','iid2','k','rMaf>30')
            summary_maf30 = subset(summary_maf30, select = -3)
            setkeyv(summary_maf30, c('iid1','iid2'))
            setkeyv(res2, c('iid1','iid2'))
            res3 = res2[summary_maf30, nomatch=0]
        }else{
            summary_maf30 = readRDS(maf30[1])
            for(j in 2:length(maf30))
                summary_maf30 = rbind(summary_maf30, readRDS(maf30[j]))
            saveRDS(summary_maf30, paste0(outputDir,'/',prefix,'_all_maf>30.RDS'))
            # write_fst(summary_maf30, paste0(outputDir,'/',prefix,'_all_maf>30.fst'))
            summary_maf30 = as.data.table(summary_maf30)
            names(summary_maf30) = c('iid1','iid2','k','rMaf>30')
            summary_maf30 = subset(summary_maf30, select = -3)
            setkeyv(summary_maf30, c('iid1','iid2'))
            setkeyv(res2, c('iid1','iid2'))
            res3 = res2[summary_maf30, nomatch=0]
        }
        
        maf40 = list.files(directoriesList[i], pattern = "*_0p4000.RDS", full.names = TRUE)
        if(length(maf40)==0){
            message("There was no file obtained for maf cutoff=40% for ", prefix)
            # next
        }else if(length(maf40)==1){
            summary_maf40 = readRDS(maf40[1])
            saveRDS(summary_maf40, paste0(outputDir,'/',prefix,'_all_maf>40.RDS'))
            # write_fst(summary_maf40, paste0(outputDir,'/',prefix,'_all_maf>40.fst'))
            summary_maf40 = as.data.table(summary_maf40)
            names(summary_maf40) = c('iid1','iid2','k','rMaf>40')
            summary_maf40 = subset(summary_maf40, select = -3)
            setkeyv(summary_maf40, c('iid1','iid2'))
            setkeyv(res3, c('iid1','iid2'))
            res4 = res3[summary_maf40, nomatch=0]
        }else if(length(maf40)>1){
            summary_maf40 = readRDS(maf40[1])
            for(j in 2:length(maf40))
                summary_maf40 = rbind(summary_maf40, readRDS(maf40[j]))
            saveRDS(summary_maf40, paste0(outputDir,'/',prefix,'_all_maf>40.RDS'))
            # write_fst(summary_maf40, paste0(outputDir,'/',prefix,'_all_maf>40.fst'))
            summary_maf40 = as.data.table(summary_maf40)
            names(summary_maf40) = c('iid1','iid2','k','rMaf>40')
            summary_maf40 = subset(summary_maf40, select = -3)
            setkeyv(summary_maf40, c('iid1','iid2'))
            setkeyv(res3, c('iid1','iid2'))
            res4 = res3[summary_maf40, nomatch=0]
        }
        # write_fst(res4, paste0(outputDir,'/',prefix,'_all_Relatedness.fst'))
        saveRDS(res4, paste0(outputDir,'/',prefix,'_all_Relatedness.RDS'))
        # system(sprintf("rm -f %s", paste0(directoriesList[i],'/*.RDS')))
    }
    stopCluster(cl)
}


summarizeFiles()