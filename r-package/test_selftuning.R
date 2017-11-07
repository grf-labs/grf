
# http://archive.ics.uci.edu/ml/datasets/Airfoil+Self-Noise
dataset = read.csv("datasets/AirFoil.txt", header=FALSE, sep="\t");
Y = dataset[,1]
X = dataset[,-1]

grf::selftuning(X,Y)