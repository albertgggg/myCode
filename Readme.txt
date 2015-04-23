This is a Matlab code of dependent multi-output Gaussian process modelling. 
The maths is explained in the paper 'Dependent Gaussian Process' by Boyle and Frean.
This folder contains the following matlab code files:
1. mainprg.m  generates training points from the models, trains a GP model which is then used for prediction.
2. trnGP.m trains a GP model with the training points from the mainprg file. Hyper paremters of GP is evaluated
   using maximum likelihood estimation (MLE).
3. predGP.m take the trained model for predictions.
4. Csub.m computes the correlation coefficient matrix for the GP model.
5. Direct.m contains the optimization algorithm 'DIviding RECTangle' method, which is used for MLE.
