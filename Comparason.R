library("microbenchmark")
#This shows a comparason of two algorithms

Particles = 40;
BoxSize = 1; #The Boxsize has the unit of Amstrong

LJPotential = function(Particles,BoxSize){
  
  #Generate Random Distribution of Particles.
  Initialize = function(particles,boxsize){
    Location = matrix(nrow = Particles, ncol = 3)
    Location[,1] = abs(runif(Particles)) * BoxSize;
    Location[,2] = abs(runif(Particles)) * BoxSize;
    Location[,3] = abs(runif(Particles)) * BoxSize;
    return(Location)
  }
  
  Location = Initialize(Particles,BoxSize);
  
  #The Function that Calculate the Distance.
  Distance = function(x,y){
    if(length(x) == length(y) & length(x) == 3){
      distance = sqrt((x[1] - y[1])^2 + (x[2] - y[2])^2 + (x[3] - y[3])^2);
      return(distance);
    }
    else{
      print("This is not the right input.")
    }
  }
  
  #Generate a Matrix With Regard of Distance
  print("Calculating the Distance Between Particles...")
  Dismatrix = matrix(ncol = Particles,nrow = Particles)
  for(ii in 1: Particles){
    for(i in 1:Particles){
      Dismatrix[i,ii] = Distance(Location[i,],Location[ii,])
    }
  }
  
  #Get Rid of half of the Redundent Component of the Matrix
  print("Optimizing...")
  for(ii in 1:Particles){
    for(i in 1:Particles){
      if(i > ii){
        Dismatrix[i,ii] = Dismatrix[i,ii];
      }
      if(i < ii){
        Dismatrix[i,ii] = 0;
      }
    }
  }
  
  #Generate the Distance Vector
  Disvec = as.vector(Dismatrix);
  Disvec = Disvec[Disvec != 0];
  
  #Function to Calculate Lennard-Jones Potential Energy
  Lennard = function(x){
    Kb = 1.38064852 * 10^(-23)
    Ipsilon = 35.61;
    Sigma = 2.75;
    4 * Ipsilon * Kb * ((Sigma/x)^12 - (Sigma/x)^6)
  }
  
  Potential = sum(Lennard(Disvec));
  
  Char = paste("The Lennard-Jones Potential of the Given System is: ", 
               round(Potential,5), " J")
  return(Char)
}
LJPotentialPeriodic = function(Particles,BoxSize){
  
  #Generate Random Distribution of Particles.
  Initialize = function(particles,boxsize){
    Location = matrix(nrow = Particles, ncol = 3)
    Location[,1] = abs(runif(Particles)) * BoxSize;
    Location[,2] = abs(runif(Particles)) * BoxSize;
    Location[,3] = abs(runif(Particles)) * BoxSize;
    return(Location)
  }
  
  Location = Initialize(Particles,BoxSize);
  
  #Generate Boundary
  Boundary = function(locavec){
    
    border = vector();
    index = c(-1,0,1);
    indexmat = vector();
    for(i in 1:3){
      for(j in 1:3){
        for(k in 1:3){
          indexmat = rbind(indexmat, c(index[i],index[j],index[k]))
        }
      }
    }
    rm(i,j,k);
    
    for(i in 1:dim(indexmat)[1]){
      temlocation = locavec;
      for(j in 1:3){
        temlocation[j] = temlocation[j] + indexmat[i,j] * BoxSize;
        border = rbind(border,temlocation);
      }
    }
    return(border)
  }
  
  Surrounding = vector();
  print("Generating Boundary Condition...")
  for(i in 1:Particles){
    Surrounding = rbind(Surrounding, Boundary(Location[i,]))
  }
  
  #The Function that Calculate the Distance.
  Distance = function(x,y){
    if(length(x) == length(y) & length(x) == 3){
      distance = sqrt((x[1] - y[1])^2 + (x[2] - y[2])^2 + (x[3] - y[3])^2);
      return(distance);
    }
    else{
      print("This is not the right input.")
    }
  }
  
  #Generate a Matrix of Distance
  print("Calculating the Distance Between the Molecules...")
  DimSur = dim(Surrounding)[1];
  Dismat = matrix(nrow = DimSur,ncol = Particles)
  for(i in 1:Particles){
    for(j in 1:DimSur){
      Dismat[j,i] = Distance(Location[i,],Surrounding[j,]);
    }
  }
  
  #Get Rid of Redundant Particles and Invalid Values
  print("Optimizing...")
  LongDis = sqrt(3 * BoxSize^2);
  Disvec = Dismat[Dismat <= LongDis];
  Disvec = Disvec[Disvec != 0];
  
  #Function to Calculate Lennard-Jones Potential Energy
  Lennard = function(x){
    Kb = 1.38064852 * 10^(-23)
    Ipsilon = 35.61;
    Sigma = 2.75;
    4 * Ipsilon * Kb * ((Sigma/x)^12 - (Sigma/x)^6)
  }
  
  #Calculate the Lennard-Jones Potential Energy of the System
  print("Calculating the Lennard-Jones Potential Energy...")
  Potential = sum(Lennard(Disvec));
  
  Char = paste("The Lennard-Jones Potential of the Given System is: ", 
               round(Potential,5), " J")
  print(Char);
  return(Potential)
}

LJPotential(Particles,BoxSize)
LJPotentialPeriodic(Particles,BoxSize)

microbenchmark::microbenchmark(LJPotential(Particles,BoxSize),
                               LJPotentialPeriodic(Particles,BoxSize), times = 500L)