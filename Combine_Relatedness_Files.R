
combineRelatednessFiles = function()
{
    arguments = commandArgs(trailingOnly = TRUE)
    pathToOutput = as.character(arguments[1])
    pathToPhysicalDistanceFile = as.character(arguments[2])
    
    require(data.table)
    require(doParallel)
    
    tmpDir = pathToOutput
    hrpiDir = list.dirs(path = pathToOutput, full.names = TRUE)
    hrpiDir = hrpiDir[2:length(hrpiDir)]
    cl = makeCluster(5)
    registerDoParallel(cl)
    foreach(inc = 1:length(hrpiDir), .packages = "data.table") %dopar%{
        pathToOutput = hrpiDir[inc]
        detectUnique = function(x){unlist(strsplit(x,'_maf>'))[1]}
        physicalDistance = fread(pathToPhysicalDistanceFile)
        relatednessFiles = list.files(pathToOutput, pattern = "_All_Pairs.RDS", full.names = TRUE)
        # m = which(grepl('_Percentage_HRI',relatednessFiles)==TRUE)
        # relatednessFiles = relatednessFiles[-m]
        # cl = makeCluster(10)
        # registerDoParallel(cl)
        message('\ncombining the data together')
        for(i in 1:length(relatednessFiles)){
            if(i==1){
                comb = readRDS(relatednessFiles[i])
                sitePairs = unique(as.character(sapply(names(comb), detectUnique))) 
                if(length(sitePairs)==1){
                    relatednessMaf0 = comb[[1]]
                    relatednessMaf10 = comb[[2]]
                    relatednessMaf20 = comb[[3]]
                    relatednessMaf30 = comb[[4]]
                    relatednessMaf40 = comb[[5]]
                }else{
                    m = which(grepl(sitePairs[1],names(comb))==TRUE)
                    relatednessMaf0 = comb[[m[1]]]
                    relatednessMaf10 = comb[[m[2]]]
                    relatednessMaf20 = comb[[m[3]]]
                    relatednessMaf30 = comb[[m[4]]]
                    relatednessMaf40 = comb[[m[5]]]
                    for(j in 2:length(sitePairs)){
                        m = which(grepl(sitePairs[j],names(comb))==TRUE)
                        relatednessMaf0 = rbind(relatednessMaf0, comb[[m[1]]])
                        relatednessMaf10 = rbind(relatednessMaf10, comb[[m[2]]])
                        relatednessMaf20 = rbind(relatednessMaf20, comb[[m[3]]])
                        relatednessMaf30 = rbind(relatednessMaf30, comb[[m[4]]])
                        relatednessMaf40 = rbind(relatednessMaf40, comb[[m[5]]])
                    }
                }
            }else{
                comb = readRDS(relatednessFiles[i])
                sitePairs = unique(as.character(sapply(names(comb), detectUnique))) 
                if(length(sitePairs)==1){
                    relatednessMaf0 = rbind(relatednessMaf0, comb[[1]])
                    relatednessMaf10 = rbind(relatednessMaf10, comb[[2]])
                    relatednessMaf20 = rbind(relatednessMaf20, comb[[3]])
                    relatednessMaf30 = rbind(relatednessMaf30, comb[[4]])
                    relatednessMaf40 = rbind(relatednessMaf40, comb[[5]])
                }else{
                    m = which(grepl(sitePairs[1],names(comb))==TRUE)
                    relatednessMaf0 = rbind(relatednessMaf0, comb[[m[1]]])
                    relatednessMaf10 = rbind(relatednessMaf10, comb[[m[2]]])
                    relatednessMaf20 = rbind(relatednessMaf20, comb[[m[3]]])
                    relatednessMaf30 = rbind(relatednessMaf30, comb[[m[4]]])
                    relatednessMaf40 = rbind(relatednessMaf40, comb[[m[5]]])
                    for(j in 2:length(sitePairs)){
                        m = which(grepl(sitePairs[j],names(comb))==TRUE)
                        relatednessMaf0 = rbind(relatednessMaf0, comb[[m[1]]])
                        relatednessMaf10 = rbind(relatednessMaf10, comb[[m[2]]])
                        relatednessMaf20 = rbind(relatednessMaf20, comb[[m[3]]])
                        relatednessMaf30 = rbind(relatednessMaf30, comb[[m[4]]])
                        relatednessMaf40 = rbind(relatednessMaf40, comb[[m[5]]])
                    }
                }
            }
            # progress(i, length(relatednessFiles))
        }
        
        
        message('\nadding physical and temporal distance on rMAF>0')
        relatednessMaf0$physicalDistance=0
        relatednessMaf0$temporalDistance=0
        relatednessMaf0$temporalDistance = abs(relatednessMaf0$year1-relatednessMaf0$year2)
        unique1 = unique(relatednessMaf0$site1)
        unique2 = unique(relatednessMaf0$site2)
        for(l in 1:length(unique1)){
            site1 = unique1[l]
            for(j in 1:length(unique2)){
                site2 = unique2[j]
                pdistance = physicalDistance[which(physicalDistance$col2==site1 & physicalDistance$col3==site2),]$distance
                if(length(pdistance)==0){
                    pdistance = physicalDistance[which(physicalDistance$col2==site2 & physicalDistance$col3==site1),]$distance
                    relatednessMaf0[which(relatednessMaf0$site1==unique1[l] & relatednessMaf0$site2==unique2[j]),]$physicalDistance=pdistance
                    #relatednessMaf0$physicalDistance[i]=pdistance
                }else{
                    relatednessMaf0[which(relatednessMaf0$site1==unique1[l] & relatednessMaf0$site2==unique2[j]),]$physicalDistance=pdistance
                }
            }
            # progress(i, length(unique1))
        }
        #saveRDS(relatednessMaf0, "/media/Data/Data/Documents_Karim/DELGEME/Last_Try/Data/From_Raw/All_pairs.RDS")
        
        
        message('\nadding physical and temporal distance on rMAF>10')
        relatednessMaf10$physicalDistance=0
        relatednessMaf10$temporalDistance=0
        relatednessMaf10$temporalDistance = abs(relatednessMaf10$year1-relatednessMaf10$year2)
        unique1 = unique(relatednessMaf10$site1)
        unique2 = unique(relatednessMaf10$site2)
        for(i in 1:length(unique1)){
            site1 = unique1[i]
            for(j in 1:length(unique2)){
                site2 = unique2[j]
                pdistance = physicalDistance[which(physicalDistance$col2==site1 & physicalDistance$col3==site2),]$distance
                if(length(pdistance)==0){
                    pdistance = physicalDistance[which(physicalDistance$col2==site2 & physicalDistance$col3==site1),]$distance
                    relatednessMaf10[which(relatednessMaf10$site1==unique1[i] & relatednessMaf10$site2==unique2[j]),]$physicalDistance=pdistance
                    #relatednessMaf0$physicalDistance[i]=pdistance
                }else{
                    relatednessMaf10[which(relatednessMaf10$site1==unique1[i] & relatednessMaf10$site2==unique2[j]),]$physicalDistance=pdistance
                }
            }
            # progress(i, length(unique1))
        }
        
        
        
        message('\nadding physical and temporal distance on rMAF>20')
        relatednessMaf20$physicalDistance=0
        relatednessMaf20$temporalDistance=0
        relatednessMaf20$temporalDistance = abs(relatednessMaf20$year1-relatednessMaf20$year2)
        unique1 = unique(relatednessMaf20$site1)
        unique2 = unique(relatednessMaf20$site2)
        for(i in 1:length(unique1)){
            site1 = unique1[i]
            for(j in 1:length(unique2)){
                site2 = unique2[j]
                pdistance = physicalDistance[which(physicalDistance$col2==site1 & physicalDistance$col3==site2),]$distance
                if(length(pdistance)==0){
                    pdistance = physicalDistance[which(physicalDistance$col2==site2 & physicalDistance$col3==site1),]$distance
                    relatednessMaf20[which(relatednessMaf20$site1==unique1[i] & relatednessMaf20$site2==unique2[j]),]$physicalDistance=pdistance
                    #relatednessMaf0$physicalDistance[i]=pdistance
                }else{
                    relatednessMaf20[which(relatednessMaf20$site1==unique1[i] & relatednessMaf20$site2==unique2[j]),]$physicalDistance=pdistance
                }
            }
            # progress(i, length(unique1))
        }
        
        
        message('\nadding physical and temporal distance on rMAF>30')
        relatednessMaf30$physicalDistance=0
        relatednessMaf30$temporalDistance=0
        relatednessMaf30$temporalDistance = abs(relatednessMaf30$year1-relatednessMaf30$year2)
        unique1 = unique(relatednessMaf30$site1)
        unique2 = unique(relatednessMaf30$site2)
        for(i in 1:length(unique1)){
            site1 = unique1[i]
            for(j in 1:length(unique2)){
                site2 = unique2[j]
                pdistance = physicalDistance[which(physicalDistance$col2==site1 & physicalDistance$col3==site2),]$distance
                if(length(pdistance)==0){
                    pdistance = physicalDistance[which(physicalDistance$col2==site2 & physicalDistance$col3==site1),]$distance
                    relatednessMaf30[which(relatednessMaf30$site1==unique1[i] & relatednessMaf30$site2==unique2[j]),]$physicalDistance=pdistance
                    #relatednessMaf0$physicalDistance[i]=pdistance
                }else{
                    relatednessMaf30[which(relatednessMaf30$site1==unique1[i] & relatednessMaf30$site2==unique2[j]),]$physicalDistance=pdistance
                }
            }
            # progress(i, length(unique1))
        }
        
        
        message('\nadding physical and temporal distance on rMAF>40')
        relatednessMaf40$physicalDistance=0
        relatednessMaf40$temporalDistance=0
        relatednessMaf40$temporalDistance = abs(relatednessMaf40$year1-relatednessMaf40$year2)
        unique1 = unique(relatednessMaf40$site1)
        unique2 = unique(relatednessMaf40$site2)
        for(i in 1:length(unique1)){
            site1 = unique1[i]
            for(j in 1:length(unique2)){
                site2 = unique2[j]
                pdistance = physicalDistance[which(physicalDistance$col2==site1 & physicalDistance$col3==site2),]$distance
                if(length(pdistance)==0){
                    pdistance = physicalDistance[which(physicalDistance$col2==site2 & physicalDistance$col3==site1),]$distance
                    relatednessMaf40[which(relatednessMaf40$site1==unique1[i] & relatednessMaf40$site2==unique2[j]),]$physicalDistance=pdistance
                    #relatednessMaf0$physicalDistance[i]=pdistance
                }else{
                    relatednessMaf40[which(relatednessMaf40$site1==unique1[i] & relatednessMaf40$site2==unique2[j]),]$physicalDistance=pdistance
                }
            }
            # progress(i, length(unique1))
        }
        combineDir = paste0(tmpDir,"/Combined")
        system(sprintf("mkdir -p %s", combineDir))
        tDir = unlist(strsplit(pathToOutput,"/"))[length(unlist(strsplit(pathToOutput,"/")))]
        cutoff = unlist(strsplit(tDir,"_"))[2]
        saveRDS(list(relatednessMaf0, relatednessMaf10,relatednessMaf20,relatednessMaf30,relatednessMaf40), paste0(combineDir,"/Combined_",cutoff,".RDS"))
    }
    stopCluster(cl)
}


combineRelatednessFiles()

