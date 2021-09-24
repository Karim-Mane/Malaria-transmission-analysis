
run_relatedness = function(){
    args = commandArgs(trailingOnly = TRUE)
    inputDir = as.character(args[1])
    outputDir = as.character(args[2])
    popOrder = as.character(args[3])
    dataRecodingMethod = as.character(args[4])
    
    require(data.table)
    require(doParallel)
    
    if(!dir.exists(inputDir)){
        stop(inputDir, "does not exist!")
    }
    if(!dir.exists(outputDir)){
        system(sprintf("mkdir -p %s", outputDir))
    }
    if(!file.exists(popOrder)){
        stop(popOrder, "no such file or directory!")
    }else{
        pop = fread(popOrder, header = FALSE)
    }
    
    pop_pairs = makeCombination(pop$V1)
    #pop_pairs = pop_pairs[3:8,]
    #print(pop_pairs)
    pwd = getwd()
    cl = makeCluster(5)
    registerDoParallel(cl)
    foreach(i = 1:nrow(pop_pairs), .export = "createYaml") %dopar%{
        # print(paste0("i=",i))
        yamlFile = paste0(pwd,"/run_relatedness_",i,".yaml")
        createYaml(inputDir, outputDir, pwd, pop_pairs$pop1[i], pop_pairs$pop2[i], dataRecodingMethod, yamlFile)
        simFunction = paste0(pwd,'/runSim.sh')
        system(sprintf("%s %s", simFunction, yamlFile))
    }
    stopCluster(cl)
}





makeCombination = function(populations){
    pop_combinations = data.frame(pop1=character(),pop2=character(), stringsAsFactors = FALSE)
    k=1
    for(i in 1:length(populations)){   #for(i in 1:(length(populations)-1)){
        for(j in i:length(populations)){  #for(j in (i+1):length(populations)){
            pop_combinations[k,1] = populations[i]
            pop_combinations[k,2] = populations[j]
            k=k+1
        }
    }
    return(pop_combinations)
}


createYaml = function(indir, outdir, pwd, pop1, pop2, recodingMethod, yamlFile){
    myFile = yamlFile #paste0(pwd,"/run_relatedness.yaml")
    f = file(myFile, open = "w")
    cat("inputDirectory: ", file = f, sep = " ")
    cat(indir, file = f, sep = "\n")
    cat("outputDirectory: ", file = f, sep = " ")
    cat(outdir, file = f, sep = "\n")
    cat("countryA: ", file = f, sep = " ")
    cat(pop1, file = f, sep = "\n")
    cat("countryB: ", file = f, sep = " ")
    cat(pop2, file = f, sep = "\n")
    cat("type: ", file = f, sep = " ")
    cat(recodingMethod, file = f, sep = "\n")
    cat("freq_start: ", file = f, sep = " ")
    cat(0, file = f, sep = "\n")
    cat("freq_end: ", file = f, sep = " ")
    cat(0.4, file = f, sep = "\n")
    cat("freq_int: ", file = f, sep = " ")
    cat(0.1, file = f, sep = "\n")
}





run_relatedness()

