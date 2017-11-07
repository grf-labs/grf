load_datasets = function(dataname) {
  if (dataname == "winequality-white") {
    dataset = read.csv("datasets/WineQuality/winequality-white.csv",sep = ";");
    
    covariate.names <- names(dataset)[1:11]
    names(dataset)[12] <- "Y"
    
  } else if (dataname == "winequality-red") {
    dataset = read.csv("datasets/WineQuality/winequality-red.csv",sep = ";");
    
    covariate.names <- names(dataset)[1:11]
    names(dataset)[12] <- "Y"
    
  } else if (dataname == "OnlineNewsPopularity") {
    dataset = read.csv("datasets/OnlineNewsPopularity/OnlineNewsPopularity.csv");
    
    covariate.names <- names(dataset)[3:60]
    names(dataset)[names(dataset)=="shares"] <- "Y"
    
  } else if (dataname == "BeijingPM2.5") {
    dataset = read.csv("datasets/PRSA_data_2010.1.1-2014.12.31.csv");
    
    nanY = is.na(dataset$pm2.5)
    dataset = dataset[!nanY,]
    dataset$cbwd = as.numeric(factor(dataset$cbwd))
    covariate.names <- names(dataset)[-c(1,6)]
    names(dataset)[names(dataset)=="pm2.5"] <- "Y"
    
  } else if (dataname == "EnergyData") {
    dataset = read.csv("datasets/energydata_complete.csv");
    
    covariate.names <- names(dataset)[3:29]
    print(covariate.names)
    names(dataset)[names(dataset)=="Appliances"] <- "Y"
    
  } else if (dataname == "CommunityCrimes") {
    # http://archive.ics.uci.edu/ml/datasets/Communities+and+Crime
    dataset = read.csv("datasets/communities.txt", header=FALSE);
    
    is.na(dataset) = dataset=="?"
    dataset[6:127] = lapply(dataset[6:127], as.numeric)
    covariate.names <- names(dataset)[6:127]
    names(dataset)[128] <- "Y"
    
  } else if (dataname == "AirFoil") {
    # http://archive.ics.uci.edu/ml/datasets/Airfoil+Self-Noise
    dataset = read.csv("datasets/airfoil_self_noise.txt", header=FALSE, sep="");
    
    covariate.names <- names(dataset)[1:5]
    names(dataset)[6] <- "Y"
    
  } else if (dataname == "BikeSharing_day") {
    dataset = read.csv("datasets/Bike-Sharing-Dataset/day.csv");
    
    is.na(dataset) = dataset=="?"
    dataset[2:15] = lapply(dataset[2:15], as.numeric)
    covariate.names <- names(dataset)[2:15]
    names(dataset)[names(dataset)=="cnt"] <- "Y"
    
  } else if (dataname == "BikeSharing_hour") {
    dataset = read.csv("datasets/Bike-Sharing-Dataset/hour.csv");
    
    is.na(dataset) = dataset=="?"
    dataset[2:16] = lapply(dataset[2:16], as.numeric)
    covariate.names <- names(dataset)[2:16]
    names(dataset)[names(dataset)=="cnt"] <- "Y"
  } else if (dataname == "TomsHardware") {
    dataset = read.csv("datasets/TomsHardware/TomsHardware.csv",header=FALSE);
    
    covariate.names <- names(dataset)[1:96]
    names(dataset)[97] <- "Y"
  } else if (dataname == "FBComments") {
    dataset = read.csv("datasets/FB_Variant_5.csv",header=FALSE);
    
    covariate.names <- names(dataset)[1:53]
    names(dataset)[54] <- "Y"
  } 
  
  # Extract the dependent variable
  Y <- dataset[["Y"]]
  covariates <- dataset[covariate.names]
  covariates.scaled <- scale(covariates)
  new_dataset = cbind(Y, covariates.scaled);
  return(new_dataset)
}

dataset = load_datasets("AirFoil")
Y = dataset[,1]
X = dataset[,-1]

dataset = read.csv("AirFoil.txt", header=FALSE, sep="\t");
Y = dataset[,1]
X = dataset[,-1]
