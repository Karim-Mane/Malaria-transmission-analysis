

gen_mles=function()
{
    require(Rcpp)
    require(Rfast)
    require(yaml)
    
    args = commandArgs(trailingOnly = TRUE)
    yamlFile = as.character(args[1])
    iloop = as.integer(args[2])  #this is the SLURM array index
    workingDirectory = as.character(args[3])
    print(args)
    
    #config = yaml::yaml.load_file(paste0(workingDirectory,'/',yamlFile))
    config = yaml::yaml.load_file(yamlFile)
    
    #outFilePrefix=config$out_file_prefix #inputDir = config$inputDirectory;
    countryA=config$countryA; countryB=config$countryB
    #print(paste0('countryA=',countryA)); print(paste0('countryB=',countryB))
    
    outFilePrefix = paste0(countryA,'_',countryB)
    type = config$type
    freq_start = config$freq_start; freq_end = config$freq_end; freq_int = config$freq_int
    fileList = list.files(path = config$inputDirectory, pattern = '.RDS', full.names = TRUE)
    countryA = fileList[which(grepl(countryA, fileList)==TRUE)]
    countryB = fileList[which(grepl(countryB, fileList)==TRUE)]
    print(paste0("Country A=",countryA));print(paste0("Country B=",countryB))
    
    system(sprintf("mkdir -p %s", config$outputDirectory))
    out_file_root = paste0(config$outputDirectory,'/From_',type,'_Data')
    system(sprintf("mkdir -p %s", out_file_root))
    source(paste0(workingDirectory,"/simulate_data.R")) # Download this script from https://github.com/artaylor85/PlasmodiumRelatedness
    sourceCpp(paste0(workingDirectory,"/hmmloglikelihood.cpp")) # Download this script from https://github.com/artaylor85/PlasmodiumRelatedness
    
    nproc=700
    epsilon = 0.001 # Fix epsilon throughout
    nboot = 100 # For CIs
    Ps = c(0.025, 0.975) # CI quantiles
    
    outputDir = paste0(out_file_root,'/',outFilePrefix)
    system(sprintf("mkdir -p %s", outputDir))
    Ff=seq(freq_start,freq_end,freq_int)
    for (f in Ff){
        if (countryA==countryB){
            # within country comparison
            L=run_country(countryA,f)
            data_set=L$data_set
            individual_names=L$individual_names
            nindividuals=L$nindividuals
            chrom=L$chrom
            pos=L$pos
            name_combinations <- matrix(nrow = nindividuals*(nindividuals-1)/2, ncol = 2)
            Y=matrix(data=NA,nrow =length(name_combinations),ncol = 4 )
            k=0
            for ( i in 1 : (nindividuals-1)){
                j=(1+k):(k+nindividuals-i)
                name_combinations[j,1]=rep(individual_names[i],each=length(j))
                name_combinations[j,2]=individual_names[(i+1):nindividuals]
                k=k+length(j)      # within country comparison
            }
        }else{
            # within country comparison
            LA=run_country(countryA,f)
            individual_names_A=LA$individual_names
            nindividuals_A=LA$nindividuals
            # # within country comparison
            LB=run_country(countryB,f)
            individual_names_B=LB$individual_names
            nindividuals_B=LB$nindividuals
            c_offset = dim(LA$data_set)[2]
            L=run_2country(countryA,countryB,f)
            data_set=L$data_set
            individual_names=L$individual_names
            nindividuals=L$nindividuals
            chrom=L$chrom
            pos=L$pos
            name_combinations <- matrix(nrow = nindividuals_A*nindividuals_B, ncol = 2)
            Y=matrix(data=NA,nrow =length(name_combinations),ncol = 4 )
            name_combinations[,1]=rep(individual_names_A,each=nindividuals_B)
            name_combinations[,2]=rep(individual_names_B,nindividuals_A)
        } 
        X=as.matrix(data_set)
        data_set$fs = rowMeans(data_set, na.rm = TRUE) # Calculate frequencies
        data_set$pos =pos
        data_set$chrom=chrom
        data_set$dt <- c(diff(data_set$pos), Inf)
        pos_change_chrom <- 1 + which(diff(data_set$chrom) != 0) # find places where chromosome changes
        data_set$dt[pos_change_chrom-1] <- Inf
        # note NA result is undocumented - could change
        a0=Rfast::rowCountValues(X, rep(0,dim(X)[1]))  #getting the count of 0 on each row (SNPs)
        a1=Rfast::rowCountValues(X, rep(1,dim(X)[1]))  #getting the count of 1 on each row (SNPs)
        a2=Rfast::rowCountValues(X, rep(2,dim(X)[1]))  #getting the count of 2 on each row (SNPs)
        ana=rowSums(is.na(X))
        frequencies=cbind(a0,a1,a2)/(dim(X)[2]-ana)
        if (all.equal(rowSums(frequencies),rep(1,dim(X)[1]))!=TRUE){
            cat(paste0("frequency ERROR"))
            return()
        }
        if (iloop < nproc){
            N=floor(dim(name_combinations)[1]/nproc)
            starty=(iloop-1)*N+1
            endy=iloop*N
        }else{
            N=dim(name_combinations)[1]-(nproc-1)*floor(dim(name_combinations)[1]/nproc)
            starty = 1+(nproc-1)*floor(dim(name_combinations)[1]/nproc)
            endy=dim(name_combinations)[1]
        }
        
        X=matrix(data=NA,nrow =endy-starty+1,ncol = 4 )
        for (icombination in starty:endy){
            #cat(paste0("icombination=",icombination,"\n"))
            individual1 <- name_combinations[icombination,1]
            individual2 <- name_combinations[icombination,2]
            if (countryA==countryB){
                # Indices of pair
                i1 = which(individual1 == names(data_set))  #index of ind1 on the data frame
                i2 = which(individual2 == names(data_set)) #index of ind2 on the data frame
            }else{
                i1 = which(individual1 == individual_names_A) #index of ind1 on the data frame
                i2 = which(individual2 == individual_names_B) + c_offset #index of ind2 on the data frame
            }
            # Extract data
            subdata <- cbind(data_set[,c("fs","dt")],data_set[,c(i1,i2)])
            names(subdata) <- c("fs","dt","Yi","Yj") # note fs not used
            krhat_hmm <- compute_rhat_hmm(frequencies, subdata$dt, cbind(subdata$Yi, subdata$Yj), epsilon)
            X[icombination-starty+1,1]=individual1
            X[icombination-starty+1,2]=individual2
            X[icombination-starty+1,3:4]=krhat_hmm
            # Generate parametric bootstrap mles
            # krhats_hmm_boot = foreach(iboot = 1:nboot, .combine = rbind) %do% {
            #   cat(paste0("iboot=",iboot,"\r"))
            #   Ys_boot <- simulate_Ys_hmm(frequencies, dis tances = data_set$dt, k = krhat_hmm[1], r = krhat_hmm[2], epsilon)
            #   compute_rhat_hmm(frequencies, subdata$dt, Ys_boot, epsilon)
            # }"
            # CIs = apply(krhats_hmm_boot, 2, function(x)quantile(x, probs=Ps)) # Few seconds
            # X = data.frame('individual1' =0 individual1, 'individual2' = individual2,
            #                rhat = krhat_hmm[2], 'r2.5%' = CIs[1,2], 'r97.5%' = CIs[2,2],
            #                khat = krhat_hmm[1], 'k2.5%' = CIs[1,1], 'k97.5%' = CIs[2,1])
            #X= data.frame('individual1' = individual1, 'individual2' = individual2,khat = krhat_hmm[1],rhat = krhat_hmm[2])
        }
        #Y[starty:endy,]=X
        #saveRDS(Y,file =paste0(outputDir,'/',outFilePrefix,"_maf>",f,".rds"))
        saveRDS(X, file = paste0(outputDir,'/',outFilePrefix,"_",starty,"_",endy,"_",gsub("\\.","p",sprintf("%.4f",f)),".RDS"))
    } # maf loop
}


## Mechanism to compute MLE given fs, distances, Ys, epsilon
compute_rhat_hmm <- function(frequencies, distances, Ys, epsilon){
    #ndata <- nrow(frequencies)
    ll <- function(k, r) loglikelihood_cpp(k, r, Ys, frequencies, distances, epsilon, rho = 7.4 * 10^(-7))
    #ll = function(k, r) loglikelihood_R(k, r, Ys, frequencies, distances, epsilon, rho = 7.4 * 10^(-7))
    optimization <- optim(par = c(50, 0.5), fn = function(x) - ll(x[1], x[2]))
    rhat <- optimization$par
    return(rhat)
}



## Mechanism to generate Ys given fs, distances, k, r, epsilon
simulate_Ys_hmm <- function(frequencies, distances, k, r, epsilon){
    Ys <- simulate_data(frequencies, distances, k = k, r = r, epsilon, rho = 7.4 * 10^(-7))
    return(Ys)
}


run_country=function(countryA,f){
    matchy=function(x){regmatches(x, regexec('_(.*?)\\_', x))[[1]][2]}
    L=readRDS(countryA)
    q=data.frame(L$X)
    q=lapply(q,function(x) as.integer(x))
    Q=matrix(unlist(q), ncol = length(q[[1]]), byrow = TRUE)
    maf=colSums(Q==1,na.rm =TRUE)/colSums(Q<=1,na.rm=TRUE) #colSums(!is.na(Q))
    i=which(colSums(Q==0,na.rm =TRUE)<colSums(Q==1,na.rm =TRUE))
    maf[i]=colSums(Q[,i]==0,na.rm =TRUE)/colSums(Q[,i]<=1,na.rm=TRUE) #colSums(!is.na(Q[,i]))
    j=which(maf<=f)
    data_set = data.frame(q)[-j,]
    # Create indices for pairwise comparisons
    individual_names <- names(q)
    nindividuals <- length(individual_names)
    chrom=as.integer(unlist(lapply(L$chroms[-j],matchy)))
    pos=L$pos[-j]
    return(list(data_set=data_set,j=j,
                individual_names=individual_names,nindividuals=nindividuals,
                chrom=chrom,pos=pos))
}

run_2country=function(countryA,countryB,f){
    matchy=function(x){regmatches(x, regexec('_(.*?)\\_', x))[[1]][2]}
    Z=readRDS(countryA)
    M=readRDS(countryB)
    X=cbind(Z$X,M$X)
    samps=c(Z$samps,M$samps)
    chroms=c(Z$chroms)
    pos=c(Z$pos)#
    q=data.frame(X)
    q=lapply(q,function(x) as.integer(x))
    Q=matrix(unlist(q), ncol = length(q[[1]]), byrow = TRUE)
    maf=colSums(Q==1,na.rm =TRUE)/colSums(Q<=1,na.rm=TRUE) #colSums(!is.na(Q))
    i=which(colSums(Q==0,na.rm =TRUE)<colSums(Q==1,na.rm =TRUE))
    maf[i]=colSums(Q[,i]==0,na.rm =TRUE)/colSums(Q[,i]<=1,na.rm=TRUE) #colSums(!is.na(Q[,i]))
    j=which(maf<=f)
    data_set = data.frame(q)[-j,]
    # Create indices for pairwise comparisons
    individual_names <- names(q)
    nindividuals <- length(individual_names)
    chrom=as.integer(unlist(lapply(chroms[-j],matchy)))
    pos=pos[-j]
    return(list(data_set=data_set,j=j,
                individual_names=individual_names,nindividuals=nindividuals,
                chrom=chrom,pos=pos))
}









gen_mles()