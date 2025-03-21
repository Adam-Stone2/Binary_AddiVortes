library('truncnorm')
library('FNN')

Sampling_Z<-function(SumofAllTess,y){
  active_indexes<-which(y==1)
  non_active_indexes<-1:length(y)
  non_active_indexes<-non_active_indexes[-active_indexes]
  
  Z_values<-vector(length= length(y))
  
  for(i in active_indexes){
    Z_values[i]<-rtruncnorm(1, a=0,b=Inf, mean = SumofAllTess[i] ,sd = 1)
  }
  for(i in non_active_indexes){
    Z_values[i]<-rtruncnorm(1, a=-Inf,b=0, mean = SumofAllTess[i] ,sd = 1)
  }
  
  return(Z_values)
}

NewTess<-function(x,j,Tess,Dim,var){ #Propose a new tessellation
  
  p=runif(1,0,1) #Randomly sample p to decide the proposed modification to tessellation.
  
  DimStar=Dim # Let proposed dimension matrix equal original dimension matrix.
  TessStar=Tess #Similar for the tessellation matrix.
  
  if (p<0.2 & length(Dim[[j]])!=length(x[1,]) | length(Dim[[j]])==1 & p<0.4){ #Add a dimension if p is less then 0.2 or if p is less then 0.4 when there is only one dimension in the Tessellation due to adjustments (Supplementary Material).
    NumberOfCovariates=1:length(x[1,]) #Let NumberOfCovariates be a vector from 1 to the number of covariates considered.
    NumberOfCovariates=NumberOfCovariates[-Dim[[j]]] #Remove all values that for the covariates that are already in the tessellation.
    DimStar[[j]]<-c(Dim[[j]],sample(NumberOfCovariates,1)) # Uniformly sample a new covariate and add it to the dimension matrix.
    TessStar[[j]]=cbind(Tess[[j]],rnorm(length(Tess[[j]][,1]),0,var)) # Sample new coordinates from Normal distribution for the new dimension and add it to the Tessellation matrix.
    Modification="AD"}
  else if (p<0.4){ #Remove a dimension if p is less then 0.4.
    RemovedDim=sample(1:length(Dim[[j]]),1) #Uniformly sample the dimension to be removed.
    DimStar[[j]]=DimStar[[j]][-RemovedDim] #Remove the dimension from the dimesion Matrix.
    TessStar[[j]]=matrix(TessStar[[j]][,-RemovedDim],ncol=length(DimStar[[j]])) #Remove the coordinates in the Tessellation matrix corresponding to the dimension removed.
    Modification="RD"}
  else if (p<0.6 || p<0.8 & length(Tess[[j]][,1])==1){ #Add a centre if p is less then 0.6 or if p is less then 0.4 when there is only one center in the Tessellation due to adjustments (Supplementary Material).
    TessStar[[j]]=rbind(Tess[[j]],rnorm(length(Dim[[j]]),0,var)) #Add a new row of coordinates, sampled from a normal distribution, to the Tessellation matrix to add a center.
    Modification="AC"}
  else if (p<0.8){ #Add a centre if p is less then 0.8. 
    CenterRemoved=sample(1:length(TessStar[[j]][,1]),1) #Sample a row.
    TessStar[[j]]=matrix(TessStar[[j]][-CenterRemoved,],ncol=length(Dim[[j]])) #Remove row sampled.
    Modification="RC"}
  else if (p<0.9 || length(Dim[[j]])==length(x[1,])){ #Change a center if p is less then 0.9 or if the all the covariates are in the tessellation.
    TessStar[[j]][sample(1:length(TessStar[[j]][,1]),1),]=rnorm(length(Dim[[j]]),0,var) # Sample a row in the tessellaion matrix and change the coordinates of the centre by sampling from a normal distribution.
    Modification="Change"}
  else{ #Swop a dimension.
    NumberOfCovariates=1:length(x[1,])  #Let NumberOfCovariates be a vector from 1 to the number of covariates considered.
    NumberOfCovariates=NumberOfCovariates[-Dim[[j]]]  #Remove all values that for the covariates that are already in the tessellation.
    DimToChange=sample(1:length(Dim[[j]]),1) #Uniformly sample a dimension to change.
    DimStar[[j]][DimToChange]=sample(NumberOfCovariates,1) #Replace the Dimension to a new uniforly sampled covariate that is not already in the tessellaion.
    TessStar[[j]][,DimToChange]=rnorm(length(Tess[[j]][,1]),0,var) #Add new normally sampled coordinates new dimension added.
    Modification="Swop"}
  
  TessStar[[j]]<-matrix(TessStar[[j]],ncol=length(DimStar[[j]])) #Ensure the the Tessellation matrix is a "matrix" type.
  
  return(list(TessStar,DimStar,Modification)) #Return new proposed tessellation.
}

Indexes<-function(x,Tess,Dim){ #Gives the row (the center) of the tessellation that each obseravtion falls within.
  if (length(Tess[,1])==1){ #only 1 centre
    CellsForGivenTess=rep(1,length(x[,1]))
  }
  else{ #multiple
    CellsForGivenTess=knnx.index(Tess,matrix(x[,Dim],ncol = length(Dim)),k=1)
  }
  return(CellsForGivenTess)
}

AlphaCalculation<-function(x,Tess,Dim,j,R_ijOld,n_ijOld,R_ijNew,n_ijNew,Modification,SigmaSquaredMu,Omega,lambda_rate){ #Calculates the acceptence rate of the proposed tessellation.
  
  d=length(Dim[[j]]);
  NumCovariates=length(x[1,]);
  cStar=length(Tess[[j]][,1]);
  
  #The Log Likelihood Ratio in the acceptance ratio
  LOGlikelihoodRatio=0.5*(log(prod(n_ijOld*SigmaSquaredMu+1))-log(prod(n_ijNew*SigmaSquaredMu+1)))+((SigmaSquaredMu/(2))*(-sum((R_ijOld^2)/(n_ijOld*SigmaSquaredMu+1))+sum((R_ijNew^2)/(n_ijNew*SigmaSquaredMu+1))))
  
  #Calculating the acceptance probability for "AD"=Adding a dimension, "RD"=Removing a dimension, "AC"=Adding a center, "RC"=Removing a center, "Change"=Changing the coordinates of a center and Swopping a dimension.
  if (Modification == "AD"){ 
    TessStructure=(dbinom(d-1,NumCovariates-1,Omega/NumCovariates))/(dbinom(d-2,NumCovariates-1,Omega/NumCovariates)*(NumCovariates-d+1))
    TransitionRatio=(NumCovariates-d+1)/d;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)
    
    #Adjustments.
    if (length(Dim[[j]])==1){
      AcceptenceProb=AcceptenceProb+log(1/2)
    }
    else if (length(Dim[[j]])==NumCovariates-1){
      AcceptenceProb=AcceptenceProb+log(2)}
  }
  else if (Modification == "RD"){
    TessStructure=(dbinom(d-1,NumCovariates,Omega/NumCovariates)*(NumCovariates-d))/(dbinom(d,NumCovariates,Omega/NumCovariates))
    TransitionRatio=(d+1)/(NumCovariates-d)
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)
    
    #Adjustments.
    if (length(Dim[[j]])==NumCovariates){
      AcceptenceProb=AcceptenceProb+log(1/2)
    }
    else if (length(Dim[[j]])==2){
      AcceptenceProb=AcceptenceProb+log(2);
    }
  }
  else if(Modification == "AC"){
    TessStructure=dpois(cStar-1,lambda_rate)/dpois(cStar-2,lambda_rate)
    TransitionRatio=1/cStar;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)
    
    #Adjustments.
    if (cStar==1){
      AcceptenceProb=AcceptenceProb+log(1/2);
    }
  }
  else if (Modification == "RC"){
    TessStructure=dpois(cStar-1,lambda_rate)/dpois(cStar,lambda_rate);
    TransitionRatio=cStar+1;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)
    
    #Adjustments.
    if (cStar==2){
      AcceptenceProb=AcceptenceProb+log(2);
    }
  }
  else if (Modification == "Change"){
    TessStructure=1;
    TransitionRatio=1;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio);
  }
  else {
    TessStructure=1;
    TransitionRatio=1;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio);
  }
  return(AcceptenceProb)
}

CalculateResiduals<-function(y,x,j,SumOfAllTess,Tess,Dim,Pred,TessStar,DimStar,LastTessPred){ #A function that calculates the n-vector of partial residuals derived from a fitting process that excludes the jth tessellation.
  if (j==1){
    indexes=Indexes(x,Tess[[j]],Dim[[j]]);
    CurrentTessPred<-Pred[[j]][indexes]
    SumOfAllTess=SumOfAllTess-CurrentTessPred}
  else{
    indexes=Indexes(x,Tess[[j]],Dim[[j]]);
    CurrentTessPred<-Pred[[j]][indexes]
    SumOfAllTess=SumOfAllTess+LastTessPred-CurrentTessPred;
  }
  
  IndexesStar=Indexes(x,TessStar[[j]],DimStar[[j]]);
  R_j<-y-SumOfAllTess
  
  #Initializing Sizes
  
  R_ijOld=rep(0,length(Pred[[j]]))
  n_ijOld=rep(0,length(Pred[[j]]))
  
  for (i in 1:length(Pred[[j]])){
    R_ijOld[i]<-sum(R_j[indexes==i])
    n_ijOld[i]<-sum(indexes==i)
  }
  
  R_ijNew=rep(0,length(TessStar[[j]][,1]))
  n_ijNew=rep(0,length(TessStar[[j]][,1]))
  
  for (i in 1:length(TessStar[[j]][,1])){
    R_ijNew[i]<-sum(R_j[IndexesStar==i])
    n_ijNew[i]<-sum(IndexesStar==i)}
  
  return(list(R_ijOld,n_ijOld,R_ijNew,n_ijNew,SumOfAllTess,IndexesStar,indexes))
}

NewPredSet<-function(j,Tess,R_ijNew,n_ijNew,sigmaSquaredMu){ #Sampling the new output values for the new tessellation.
  PredSet=rep(0,length(Tess[[j]][,1]))
  for (i in 1:length(Tess[[j]][,1])){
    PredSet[i]=rnorm(1,sigmaSquaredMu*R_ijNew[i]/(sigmaSquaredMu*n_ijNew[i]+1),(sigmaSquaredMu/(n_ijNew[i]*sigmaSquaredMu+1))^0.5);
  }
  return(PredSet)
}

TestPrediction<-function(x,m,Tess,Dim,Pred){ #A function that derives a prediction for the output variable give an observation, for one iteration in AddiVortes.
  Prediction=rep(0,length(x[,1]));
  for (j in 1:m){
    NewTessIndexes=Indexes(x,Tess[[j]],Dim[[j]]);
    Prediction=Prediction+Pred[[j]][NewTessIndexes]
  }
  return(Prediction)
}

AddiVortes_Algorithm<-function(y,x,m,max_iter,burn_in,k,var,Omega,lambda_rate,YTest,XTest){
  
  xScaled=x;
  for (i in 1:length(x[1,])){
    xScaled[,i]=(x[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
  }
  
  
  for (i in 1:length(XTest[1,])){
    XTest[,i]=(XTest[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
  }
  
  Pred<-rep(list(matrix(mean(y)/m)),m)
  Dim=vector(length = m)
  Tess=vector(length = m)
  for (i in 1:m){
    Dim[i]<-list(sample(1:length(x[1,]), 1))
    Tess[i]<-(list(matrix(rnorm(1,0,var))))
  }
  
  #Prepare some variables used in the backfitting algorithm
  SumOfAllTess=rep(mean(y),length(y))
  SigmaMu=3/(k*sqrt(m))
  LastTessPred=matrix
  variable_selection_counter<-numeric(ncol(x))
  
  #Matrices that will hold the samples from the poseterior distribution for the training samples and test samples.
  PredictionMatrix<-array(dim=c(length(y),(max_iter-burn_in)))
  TestMatrix<-array(dim=c(length(YTest),(max_iter-burn_in)))
  
  for (i in 1:max_iter){
    Z_values<-Sampling_Z(SumOfAllTess,y)
    
    for (j in 1:m){
      NewTessOutput<-NewTess(xScaled,j,Tess,Dim,var) #Propose new Tessellation 
      TessStar<-NewTessOutput[[1]]  
      DimStar<-NewTessOutput[[2]]
      Modification<-NewTessOutput[[3]]
      
      ResidualsOutput<-CalculateResiduals(Z_values,xScaled,j,SumOfAllTess,Tess,Dim,Pred,TessStar,DimStar,LastTessPred) #Calculate the n-vector of partial residuals derived from a fitting process that excludes the jth tessellation and the number of observations in each cell.
      R_ijOld<-ResidualsOutput[[1]]   #Old and New refer to the original and proposed tessellations
      n_ijOld<-ResidualsOutput[[2]]
      R_ijNew<-ResidualsOutput[[3]]
      n_ijNew<-ResidualsOutput[[4]]
      SumOfAllTess<-ResidualsOutput[[5]] #Keeps track of the prediction for all tessellations to help sample sigma squared.
      IndexesStar<-ResidualsOutput[[6]] #Gives the row of each observation for the cell it falls in for the proposed tessellation.
      Indexes<-ResidualsOutput[[7]]  #Gives the row of each observation for the cell it falls in for the original tessellation.
      
      if (!any(n_ijNew==0)){ #automatically rejects proposed tessellation if there exists a cell with no observations in.
        
        LOGAcceptenceProb=AlphaCalculation(xScaled,TessStar,DimStar,j,R_ijOld,n_ijOld,R_ijNew,n_ijNew,Modification,SigmaMu^2,Omega,lambda_rate) #Gives the log of the acceptence probability.
        
        if (log(runif(n=1, min=0, max=1))<LOGAcceptenceProb){ #Accepts the proposed tessellation is accepted then calculates the new output values for the new tessellation. 
          Tess=TessStar
          Dim=DimStar
          Pred[[j]]=NewPredSet(j,TessStar,R_ijNew,n_ijNew,SigmaMu^2)
          LastTessPred=Pred[[j]][IndexesStar]
        }
        else { #Rejects the proposed tesellation then calculates new output values for the original tessellation.
          Pred[[j]]=NewPredSet(j,Tess,R_ijOld,n_ijOld,SigmaMu^2);
          LastTessPred=Pred[[j]][Indexes];
        }
      }
      else{ #Rejects the proposed tesellation then calculates new output values for the original tessellation.
        Pred[[j]]=NewPredSet(j,Tess,R_ijOld,n_ijOld,SigmaMu^2);
        LastTessPred=Pred[[j]][Indexes];
      }
      if (j==m){ #If j equals m then adds the last tessellation output values to give a prediction.
        SumOfAllTess=SumOfAllTess+LastTessPred;
      }
      if (i>burn_in){ 
        variable_selection_counter[Dim[[j]]]=variable_selection_counter[Dim[[j]]]+1
      }
    }
    
    if (i>burn_in){ #vectors that hold the predictions for each iteration after burn in.
      PredictionMatrix[,(i-burn_in)]=SumOfAllTess;
      TestMatrix[,(i-burn_in)]=TestPrediction(XTest,m,Tess,Dim,Pred);
    }
    print(i)
  }
  
  #finding the mean of the predition over the iterations and then unscaling the predictions.
  
  probabilities_test<-pnorm(TestMatrix)
  mean_yhat=(rowSums(pnorm(PredictionMatrix))/(max_iter-burn_in))
  mean_yhat_Test=(rowSums(probabilities_test)/(max_iter-burn_in))
  
  lower_quantiles<-apply(probabilities_test,1, function(row) quantile(row,probs=0.25))
  upper_quantiles<-apply(probabilities_test,1, function(row) quantile(row,probs=0.75))
  
  predicted_values<-ifelse(mean_yhat_Test>0.5,1,0)
  TP<-sum(YTest==1 & predicted_values==1)
  Precision<-TP/(TP+sum(YTest==0 & predicted_values==1))
 
  
  return( #Returns the RMSE value for the test samples.
    list(
      ProbTrain = mean_yhat,
      ProbTest =mean_yhat_Test,
      Accuracy_train = length(which(y==round(mean_yhat)))/length(y),
      Accuracy_test = length(which(YTest==round(mean_yhat_Test)))/length(YTest),
      RMSE_Train = sqrt(mean((y-round(mean_yhat))^2)),
      RMSE_Test = sqrt(mean((YTest-round(mean_yhat_Test))^2)),
      lower_quant_test=lower_quantiles,
      upper_quant_test=upper_quantiles,
      Precision= Precision,
      variable_selection_counter=variable_selection_counter
    )
  )
}


# Rotated Axes Function
rotated_axes <- function(x, theta, delta) {
  u <- c(x[1] * cos(theta) - x[2] * sin(theta),
         x[1] * sin(theta) + x[2] * cos(theta))
  return(delta * (2 * (u[1] * u[2] > 0) - 1))
}

# Sinusoid Function
sinusoid <- function(x, theta, delta) {
  return(delta * (2 * (x[2] > theta * sin(10 * x[1])) - 1))
}

set.seed(123)  # For reproducibility
n <- 10000       # Number of points
x <- matrix(runif(2 * n, -1, 1), ncol = 2)  # Uniform points in [0, 1]^2
theta_rot <- pi / 6  # Example rotation angle for rotated axes
theta_sin <- 0.5    # Example amplitude for sinusoid
delta <- 4           # Jump size for the true function

# Generate responses
y_rotated <- ifelse((apply(x, 1, rotated_axes, theta = theta_rot, delta = delta) + rnorm(n, 0, 1))>0,1,0)
y_sinusoidal <- ifelse((apply(x, 1, sinusoid, theta = theta_sin, delta = delta) + rnorm(n, 0, 1))>0,1,0)

y_rotated_true <- apply(x, 1, rotated_axes, theta = theta_rot, delta = delta)
y_sinusoidal_true <- apply(x, 1, sinusoid, theta = theta_sin, delta = delta) 

library(BayesTree)  # For axis-aligned BART
#devtools::install_github(repo = "paulhnguyen/obliqueBART", subdir = "obliqueBART")
#library(obliqueBART)  # For obliqueBART (from GitHub)

# Axis-aligned BART
#bart_rotated <- pbart(x, y_rotated)
#bart_sinusoidal <- pbart(x, y_sinusoidal)


TrainSet<-1:(n/5)
TestSet<-(n/5):n
Addi_results<-AddiVortes_Algorithm(y_rotated[TrainSet],x[TrainSet,],200,500,100,3,0.8,1,25,y_rotated[TestSet],x[TestSet,])
bart_results<-bart(x[TrainSet,],y_rotated[TrainSet],x[TestSet,])

Addi_results_sin<-AddiVortes_Algorithm(y_sinusoidal[TrainSet],x[TrainSet,],200,500,100,3,0.8,1,25,y_sinusoidal[TestSet],x[TestSet,])
bart_results_sin<-bart(x[TrainSet,],y_sinusoidal[TrainSet],x[TestSet,])

# Create a gradient palette
gradient_palette <- colorRampPalette(c("red", "yellow", "blue"))

# Map probabilities to colours
point_colours <- gradient_palette(100)[as.numeric(cut(Addi_results$ProbTest, breaks = 100))]
point_colours_bart <- gradient_palette(100)[as.numeric(cut(colMeans(bart_results$yhat.test), breaks = 100))]

par(mfrow=c(1,3),mar = c(5, 2, 6, 5))
# Plot with gradient colours
plot(x, col = ifelse(y_rotated_true > 0, "blue", "red"), main = "Rotated Axes (True)")
plot(x[TestSet,], col = point_colours, pch = 16, main = "Probability Gradient Plot", xlab = "x1", ylab = "x2")
par(mar = c(5, 2, 6, 7))
plot(x[TestSet,], col = point_colours_bart, pch = 16, main = "Probability Gradient Plot", xlab = "x1", ylab = "x2")


image.plot(
  zlim = c(0, 1),  # The range of probabilities
  legend.only = TRUE,  # Only add the legend (not the full image)
  col = gradient_palette(100),  # Use the gradient palette
  legend.lab = "Probability Of Being Active",  # Label for the legend # Distance of the label from the bar
  legend.line = 3,                 # Distance of the label from the bar
  legend.width = 1.5,              # Width of the colour bar
  legend.mar = 3
)

                         
set.seed(14)  # For reproducibility
n <- 2000       # Number of points
x <- matrix(runif(2 * n, -1, 1), ncol = 2)  # Uniform points in [0, 1]^2
theta_rot <- c(0,pi/20,pi/14,pi/10,pi/6,pi/4)  # Example rotation angle for rotated axes
theta_sin <- c(0,0.1,0.3,0.5,0.75,1)    # Example amplitude for sinusoid
delta <- 4           # Jump size for the true function

Accuracy_values<-c("Addi_rotated","BART_rotated","Addi_sinusoidal","BART_sinusoidal")
Accuracy_values_rf<-c("Random_Forest_rotated","Random_Forest_rotated","Random_Forest_sinusoidal","Random_Forest_sinusoidal")

num_cores <- 10 # Specify the number of cores you want to use
cl <- makeCluster(num_cores)
registerDoParallel(cl)

Accuracy_results<-foreach(i = 1:6, .combine = rbind) %dopar%{ 
  
  library(randomForest)
  library(xgboost)
  #library(e1071)
  library(BayesTree)
  library('truncnorm')
  library('FNN')
  
  # Generate responses
  y_rotated <- ifelse((apply(x, 1, rotated_axes, theta = theta_rot[i], delta = delta) + rnorm(n, 0, 1))>0,1,0)
  #y_rotated_true <- apply(x, 1, rotated_axes, theta = theta_rot[i], delta = delta)

  
  TrainSet<-1:(4*n/8)
  TestSet<-(4*n/8):n
  Addi_results<-AddiVortes_Algorithm(y_rotated[TrainSet],x[TrainSet,],200,1200,400,1,0.4,1,25,y_rotated[TestSet],x[TestSet,])
  bart_results<-bart(x[TrainSet,],y_rotated[TrainSet],x[TestSet,])
  rf<-randomForest(x[TrainSet,],as.factor(y_rotated[TrainSet]),ntree=1000)
  xgb<-xgboost(data = xgb.DMatrix(data = x[TrainSet,], label = y_rotated[TrainSet]), nrounds = 1000,print_every_n = 100,objective = "binary:logistic")
  #svm<-svm(x[TrainSet,],y_rotated[TrainSet],probability = TRUE)  

  xgb_accuracy_rot<-sum(y_rotated[TestSet]==ifelse(predict(xgb,xgb.DMatrix(data = x[TestSet,])) > 0.5, 1, 0))/length(TestSet)
  bart_accuracy_rot<-sum(y_rotated[TestSet]==ifelse(colMeans(bart_results$yhat.test) > 0.5, 1, 0))/length(TestSet)
  random_forest_accuracy_rot<-sum(y_rotated[TestSet]==predict(rf,x[TestSet,]))/length(TestSet)
  
  
  return(c(Addi_results$Accuracy_test,bart_accuracy_rot,random_forest_accuracy_rot,xgb_accuracy_rot))
}

Accuracy_results_sin<-foreach(i = 1:6, .combine = rbind) %dopar%{ 
  
  library(randomForest)
  library(xgboost)
  #library(e1071)
  library(BayesTree)
  library('truncnorm')
  library('FNN')

  #y_sinusoidal_true <- apply(x, 1, sinusoid, theta = theta_sin[i], delta = delta) 
  
  y_sinusoidal <- ifelse(apply(x, 1, sinusoid, theta = theta_sin[i], delta = delta)>0,1,0)
  
  TrainSet<-1:(4*n/8)
  TestSet<-(4*n/8):n
  Addi_results_sin<-AddiVortes_Algorithm(y_sinusoidal[TrainSet],x[TrainSet,],200,1200,300,1,0.4,1,25,y_sinusoidal[TestSet],x[TestSet,])
  bart_results_sin<-bart(x[TrainSet,],as.factor(y_sinusoidal[TrainSet]),x[TestSet,])
  bart_accuracy_sin<-sum(y_sinusoidal[TestSet]==ifelse(colMeans(bart_results_sin$yhat.test) > 0.5, 1, 0))/length(TestSet)
  rf<-randomForest(x[TrainSet,],as.factor(y_sinusoidal[TrainSet]),ntree=1000)
  randome_forest_accuracy_sin<-sum(y_sinusoidal[TestSet]==predict(rf,x[TestSet,]))/length(TestSet)
  xgb<-xgboost(data = xgb.DMatrix(data = x[TrainSet,], label = y_sinusoidal[TrainSet]), nrounds = 1000,print_every_n = 100,objective = "binary:logistic")
  xgb_accuracy_sin<-sum(y_sinusoidal[TestSet]==ifelse(predict(xgb,xgb.DMatrix(data = x[TestSet,])) > 0.5, 1, 0))/length(TestSet)
  
  return(c(Addi_results_sin$Accuracy_test,bart_accuracy_sin,randome_forest_accuracy_sin, xgb_accuracy_sin))
}

stopCluster(cl)

par(mfrow=c(1,2))
par(mar = c(4, 4, 4, 4))

plot(theta_rot,as.numeric(Accuracy_results[,1]),ylim = c(0.87,1),pch=16,type="b",xlab="Rotation Angle",ylab="Accuracy")
points(theta_rot,as.numeric(Accuracy_results[,2]),col="red",pch=16,type = "b")
points(theta_rot,as.numeric(Accuracy_results[,3]),col="green",pch=16,type = "b")
points(theta_rot,as.numeric(Accuracy_results[,4]),col="blue",pch=16,type = "b")


#par(mfrow=c(1,2))
plot(theta_sin,as.numeric(Accuracy_results_sin[,1]),ylim = c(0.87,1),pch=16,type="b",xlab="Amplitude",ylab="Accuracy")
points(theta_sin,as.numeric(Accuracy_results_sin[,2]),col="red",pch=16,type = "b")
points(theta_sin,as.numeric(Accuracy_results_sin[,3]),col="green",pch=16,type = "b")
points(theta_sin,as.numeric(Accuracy_results_sin[,4]),col="blue",pch=16,type = "b")
