# Three-dimensional system of lennard-JOnes Particles

#Cleaning up global variables

rm(list = ls())
setwd("D:/r/")

#Build the cluster of the particles

n = 1000; #The number of molecules in this simulation
Dimension = 3;
Cluster <- matrix(nrow = n, ncol = Dimension);

#Declare box size and put molecules into the box

Boxsize = 10; #The unit is Amstrong
XValue = abs(runif(n)) * Boxsize;
YValue = abs(runif(n)) * Boxsize;
ZValue = abs(runif(n)) * Boxsize;

Cluster[,1] = XValue;
Cluster[,2] = YValue;
Cluster[,3] = ZValue;

#Function of the calculation of the distance

Distance = function(x,y){
  if(length(x) == 3 & length(y) == 3){
    distance = ((x[1]-y[1])^2 + (x[2]-y[2])^2 + (x[3] - y[3])^2)^0.5
    return(distance)
  }
  if(length(x) == 2 & length(y) == 2){
    distance = ((x[1]-y[1])^2 + (x[2]-y[2])^2)^0.5
    return(distance)
  }
}

#Function of the Calculation of the Lennard-Jones potential

LJpotential = function(x){
  ep = 35.6;
  Sigma = 2.75;
  Kb = 1.38064852 * 10^(-23) #Boltzmann constant
  Potential = 4*ep*((Sigma/x)^12 - (Sigma/x)^6) * Kb
  return(Potential)
}

#The calculation of energy from one molecule

DistanceMatrix = matrix(nrow = n, ncol = n)

for(ii in 1:n)
for(i in 1:n){
  DistanceMatrix[i,ii] = Distance(as.vector(Cluster[i,]), as.vector(Cluster[ii,]))
}

DistanceVector = as.vector(DistanceMatrix)
DistanceVector = DistanceVector[DistanceVector != 0]

PotentialEnergy = LJpotential(DistanceVector)
TotalPotential = sum(PotentialEnergy)/2

#Graphics

split.screen(c(1,2))

screen(1)
plot(Cluster[,1],Cluster[,2], cex = 0.5, pch = 16,
     xlab = "The Box", ylab = "The Box",
     main = "The projection onto the Z and Plane")
screen(2)
plot(Cluster[,1],Cluster[,3], cex = 0.5, pch = 16,
     xlab = "The Box", ylab = "The Box",
     main = "The projection onto the Y Plane")

#Result

print(paste(TotalPotential, "J"))
rm(list = ls())





