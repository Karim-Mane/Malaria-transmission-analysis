

sumres_v2=function()
{
  require(data.table)
  # require(ggplot2)
  # require(yaml)
  require(officer)
  require(flextable)
  require(dplyr)
  require(doParallel)
  
  arguments = commandArgs(trailingOnly = TRUE)
  inputDir = as.character(arguments[1])
  # outputDir = as.character(arguments[2])
  metadataFile = as.character(arguments[2])
  # coff=as.numeric(as.character(arguments[3])) # for relatadeness
  
  meta = fread(metadataFile,header=TRUE)
  
  # allFiles = list.files(path = outputDir, pattern = '*.RDS', full.names = TRUE)
  cutoffs = c(0.2, 0.3, 0.4, 0.5, 0.6)
  
  for(inc in 1:length(cutoffs)){
    coff = cutoffs[inc]
    outputDir = paste0(inputDir,"/HRIP_",coff)  #HRIP stands for highly related infection pairs
    system(sprintf("mkdir -p %s", outputDir))
    allFiles = list.files(path = inputDir, pattern = '*_all_Relatedness.RDS', full.names = TRUE)
    cl = makeCluster(5)
    registerDoParallel(cl)
    foreach(i = 1:length(allFiles), .packages = c("dplyr","flextable","officer","data.table")) %dopar%{
          site1 = unlist(strsplit(basename(allFiles[i]),'_'))[1]
          site2 = unlist(strsplit(basename(allFiles[i]),'_'))[2]
          prefix = paste0(site1,'_',site2)
    
          target = allFiles[i]
          Xgam=readRDS(target)  #Xgam=readRDS(f)
          modifier = function(x){gsub('-','.',x)}
          meta$SampleID=lapply(meta$SampleID,modifier)
          m = match(Xgam$iid1, meta$SampleID)
          Xgam$site1=meta$Location[m]
          Xgam$year1=meta$Year[m]
          m = match(Xgam$iid2, meta$SampleID)
          Xgam$site2=meta$Location[m]
          Xgam$year2=meta$Year[m]
          
          g_gam=Xgam[,.(count = .N),by=list(site1,site2,year1,year2)]
          nams_gam=c(paste0(g_gam$site1,g_gam$year1),paste0(g_gam$site2,g_gam$year2))
          g_gam=as.matrix(g_gam)
          
          if(nrow(g_gam)==1)
          {
              Lgam=list()
              for ( j in 1 : (nrow(g_gam)+1)){
                  if(j==1)
                      Lgam[[j]]=Xgam[Xgam$site1==g_gam[1,j] & Xgam$year1==as.numeric(g_gam[1,(j+2)]),]$iid1
                  else
                      Lgam[[j]]=Xgam[Xgam$site2==g_gam[1,j] & Xgam$year2==as.numeric(g_gam[1,(j+2)]),]$iid2
              }
              
              relatedPairs=list()
              allPairs=list()
              ff=1
              nams_freq = c('maf>0','maf>10','maf>20','maf>30','maf>40')
              outlistName=NULL
              X = matrix(NA, 1, 6)   #X = matrix(NA, 1, 7)
              colnames(X)=c('sitePairs','maf>0%','maf>10%','maf>20%','maf>30%','maf>40%') #colnames(X)=c('sitePairs','maf>0%','maf>10%','maf>20%','maf>30%','maf>40%', '#Pairs')
              for (f in 1:5)
              {
                  targetRelatedness = as.matrix(subset(Xgam, select = (f+2)))
                  targetRelatedness=as.numeric(targetRelatedness[,1])
                  Yout_gam = numeric(1)
                  dum=character(1)
                  dum[1]=paste0(nams_gam[1],"_vs_",nams_gam[2])
                  r = match(Xgam$iid1,unlist(Lgam[[1]]))
                  j = match(Xgam$iid2,unlist(Lgam[[2]]))
                  k = !is.na(r) & !is.na(j)
                  Xgam_prime = Xgam[k,]
                  mm=which(targetRelatedness[k]>coff)
                  Yout_gam[1]=100*length(mm)/sum(k)
                  allPairs[[ff]] = Xgam_prime
                  relatedPairs[[ff]] = subset(Xgam_prime[mm,], select = c(1:2,(f+2),8:ncol(Xgam_prime)))
                  outlistName = c(outlistName, paste0(dum[1],'_',nams_freq[f]))
                  ff=ff+1
                  
                  X[1,(f+1)]=round(Yout_gam,2)
                  # X[1,7]=length(mm)
                  if (f==1){
                      X[1,f]=dum
                  }
              }
              names(relatedPairs) = outlistName
              names(allPairs) = outlistName
          }
          else
          {
              Lgam = list()
              comb1 = paste0(g_gam[,3],g_gam[,4])
              comb2 = paste0(g_gam[,4],g_gam[,3])
              sameComb = data.frame(index1=numeric(),index2=numeric())
              ooo=1
              saved=NULL
              for(o in 1:length(comb1))
              {
                  if(o %in% saved)
                      next
                  for(oo in 1:length(comb2))
                  {
                      if(comb1[o]==comb2[oo] & o!=oo)
                      {
                          sameComb[ooo,1]=o
                          sameComb[ooo,2]=oo
                          ooo=ooo+1
                          saved = c(saved, oo)
                      }
                  }
              }
              if(nrow(sameComb)==0)
              {
                  n=1
                  for ( j in 1 : nrow(g_gam)){
                      popTarget = Xgam[which((Xgam$site1==g_gam[j,1] & Xgam$year1==g_gam[j,3]) & (Xgam$site2==g_gam[j,2] & Xgam$year2==g_gam[j,4])),]
                      Lgam[[n]]=popTarget[popTarget$site1==g_gam[j,1] & popTarget$year1==as.numeric(g_gam[j,3]),]$iid1
                      n=n+1
                      Lgam[[n]]=popTarget[popTarget$site2==g_gam[j,2] & popTarget$year2==as.numeric(g_gam[j,4]),]$iid2
                      n=n+1
                  }
              }
              else
              {
                  saved=NULL
                  n=1
                  for ( j in 1 : nrow(g_gam)) #{  #(2*length(unique(g_gam[,3])))
                  {
                      if(j %in% saved)
                          next
                      else
                      {
                          v1 =match(j,sameComb$index1)   #which(!is.na(match(j,sameComb$index1)))
                          if(length(v1)!=0)
                          {
                              saved = c(saved, sameComb$index2[v1])
                              popTarget1 = Xgam[which((Xgam$site1==g_gam[j,1] & Xgam$year1==g_gam[j,3]) & (Xgam$site2==g_gam[j,2] & Xgam$year2==g_gam[j,4])),]
                              popTarget2 = Xgam[which((Xgam$site1==g_gam[sameComb$index2[v1],1] & Xgam$year1==g_gam[sameComb$index2[v1],3] ) & (Xgam$site2==g_gam[sameComb$index2[v1],2] & Xgam$year2==g_gam[sameComb$index2[v1],4])),]
                              popTarget = rbind(popTarget1, popTarget2)
                              lgam1=popTarget1[popTarget1$site1==g_gam[j,1] & popTarget1$year1==as.numeric(g_gam[j,3]),]$iid1
                              lgam2=popTarget2[popTarget2$site1==g_gam[sameComb$index2[v1],1] & popTarget2$year1==as.numeric(g_gam[sameComb$index2[v1],3]),]$iid1
                              Lgam[[n]]=c(lgam1, lgam2)#popTarget[popTarget$site1==g_gam[j,1] & popTarget$year1==as.numeric(g_gam[j,3]),]$iid1
                              n=n+1
                              lgam1=popTarget1[popTarget1$site2==g_gam[j,2] & popTarget1$year2==as.numeric(g_gam[j,4]),]$iid2
                              lgam2=popTarget2[popTarget2$site2==g_gam[sameComb$index2[v1],2] & popTarget2$year2==as.numeric(g_gam[sameComb$index2[v1],4]),]$iid2
                              Lgam[[n]]=c(lgam1, lgam2)#popTarget[popTarget$site2==g_gam[j,2] & popTarget$year2==as.numeric(g_gam[j,4]),]$iid2
                              n=n+1
                          }
                          else
                          {
                              popTarget = Xgam[which((Xgam$site1==g_gam[j,1] & Xgam$year1==g_gam[j,3]) & (Xgam$site2==g_gam[j,2] & Xgam$year2==g_gam[j,4])),]
                              Lgam[[n]]=popTarget[popTarget$site1==g_gam[j,1] & popTarget$year1==as.numeric(g_gam[j,3]),]$iid1
                              n=n+1
                              Lgam[[n]]=popTarget[popTarget$site2==g_gam[j,2] & popTarget$year2==as.numeric(g_gam[j,4]),]$iid2
                              n=n+1
                          }
                      }
                  }
                  
              }
              relatedPairs=list()
              allPairs=list()
              ff=1
              nams_freq = c('maf>0','maf>10','maf>20','maf>30','maf>40')
              outlistName=NULL
              
              X = matrix(NA, (nrow(g_gam)-nrow(sameComb)), 6) #X = matrix(-9, length(nams_gam_unique), 7)
              colnames(X)=c('sitePairs','maf>0%','maf>10%','maf>20%','maf>30%','maf>40%')  #colnames(X)=c('sitePairs','maf>0%','maf>10%','maf>20%','maf>30%','maf>40%', '#Pairs')
              # u_year2 = length(unique(g_gam[,4]))
              for (f in 1:5)
              {
                  targetRelatedness = as.matrix(subset(Xgam, select = (f+2)))
                  targetRelatedness=as.numeric(targetRelatedness[,1])
                  Yout_gam = numeric((nrow(g_gam)-nrow(sameComb))) #numeric(length(nams_gam_unique))
                  dum=character((nrow(g_gam)-nrow(sameComb))) #character(length(nams_gam_unique))
                  numPairs = numeric((nrow(g_gam)-nrow(sameComb)))
                  county = 0; m=1; q=1
                  while(m<=length(Lgam))    #while(m<=(length(Lgam)-u_year2))  #for ( m in 1 : nrow(g_gam))
                  {
                      if(!is.null(saved) & q %in% saved)
                          q=q+1
                      n=m+1  #n=(length(Lgam)-u_year2)+1
                      
                      county = county + 1
                      dum[county]=paste0(g_gam[q,1],g_gam[q,3],"_vs_",g_gam[q,2],g_gam[q,4]) #paste0(nams_gam_unique[m],"_vs_",nams_gam_unique[n])
                      r = match(Xgam$iid1,unlist(Lgam[[m]]))
                      j = match(Xgam$iid2,unlist(Lgam[[n]]))
                      k = !is.na(r) & !is.na(j)
                      Xgam_prime = Xgam[k,]
                      mm=which(targetRelatedness[k]>coff)
                      Yout_gam[county]=100*length(mm)/sum(k)
                      numPairs[county] = length(mm)
                      allPairs[[ff]] = Xgam_prime
                      relatedPairs[[ff]] = subset(Xgam_prime[mm,], select = c(1:2,(f+2),8:ncol(Xgam_prime)))
                      outlistName = c(outlistName, paste0(dum[county],'_',nams_freq[f]))
                      ff=ff+1; m=m+2; q=q+1
                      
                  }
                  X[1:length(Yout_gam),(f+1)]=as.numeric(round(Yout_gam,2))
                  # X[1:length(Yout_gam),7]=numPairs
                  if (f==1){
                      X[1:length(Yout_gam),f]=dum
                  }
              }
              names(relatedPairs) = outlistName
              names(allPairs) = outlistName
          }
          
          saveRDS(X, paste0(outputDir,'/',prefix,"_Percentage_HRI.RDS"))
          ft=flextable(as.data.frame(X)) %>% autofit() 
          colorToUse = c('#EFEFEF','yellow','cyan','steelblue','orange','purple','lightblue','lightgreen','pink','navy')
          for(o in 1:nrow(X)){
              ft = bg(ft, i = o, bg = colorToUse[o], part = "body")}
          ft = add_footer_row(ft, values = paste0(prefix," relatedness cut-off = ",coff), colwidths = 6, top = FALSE) 
          #ft = merge_at(ft, j = 1:5, part = "footer")# theme the tableft <- theme_box(ft)
          save_as_html(ft,path=paste0(outputDir,'/',prefix,"_MAFp2.html"),title="MAF")
          saveRDS(relatedPairs, paste0(outputDir,'/',prefix,"_HRI.RDS"))
          saveRDS(allPairs, paste0(outputDir,'/',prefix,"_All_Pairs.RDS"))
    }
    stopCluster(cl)
  }
    
}

sumres_v2()