%% Based on the paper 'Dependent Gaussian Process' -- Boyle & Frean
% Programmed by Chen Liang: Jan 2014
% When multiple outputs are tightly coupled, it is desired to build a 
% multi-output surrogate model that can take advantage of this dependency. 
% When one model have lost training points, the missed information can be
% inferred from other models according to the dependence.This code adopts 
% a Gaussian process surrogate modeling technique, which incorporate the 
% covariance between different outputs into the estimation of GP hyper 
% parameter. 

close all;
clear;clc
global counter
counter=0
% Define one dimensional input
xlo = 1; xhi = 4;
X = linspace(xlo,xhi,101).';

%% Define functions
F1= @(x) (sin(3*x-2)+(x-1).^3/11);
F2 = @(x) -(cos(4*x)+exp(x/50)/18);
F3 = @(x) tan(x/3)+sin(5*x);
F4 = @(x) (x-2.5).^2;
%% True function values
Y1 = F1(X); Y2 = F2(X);Y3 = F3(X);Y4=F4(X);
Yori = [Y1, Y2, Y3, Y4];

%% Actual training points
Ns(1)=6; Ns(2)=10;Ns(3)=12;Ns(4)=8;  % Number of training points
sigma_n1 =.1;sigma_n2 =.1;sigma_n3 =.1;sigma_n4 =.1;  % noise

% Training points
xs1 = [linspace(1,1.6,Ns(1)/2),linspace(3,4,Ns(1)/2)]'; 
xs2 = linspace(xlo,xhi,Ns(2))';
xs3 = linspace(xlo,xhi,Ns(3))';
xs4 = xs1;

ys1 = F1(xs1) + normrnd(0, sigma_n1, size(xs1));
ys2 = F2(xs2) + normrnd(0, sigma_n2, size(xs2));
ys3 = F3(xs3) + normrnd(0, sigma_n3, size(xs3));
ys4 = F4(xs4) + normrnd(0, sigma_n4, size(xs4));

xtrn = [{xs1},{xs2},{xs3},{xs4}];
ytrn = [{ys1},{ys2},{ys3},{ys4}];

%% Define hyper parameters of dependent GP
% Define no. of input and no. of output
Nin = size(xs1,2);
Nout = size(ytrn,2);
%% Define upper& lower bounds of the hyper parameters
v_lolim = -1.5*ones(Nout,1)';
v_hilim =  1.5*ones(Nout,1)';
w_lolim = -1.5*ones(Nout,1)';
w_hilim =  1.5*ones(Nout,1)';
f_lolim =  0*ones(Nout,1)';
f_hilim =  3*ones(Nout,1)';
g_lolim =  0*ones(Nout,1)';
g_hilim =  3*ones(Nout,1)';
Beta_lolim =  -6*ones(Nout,1)';
Beta_hilim =  -2*ones(Nout,1)';
mu_lolim =  -1.5*ones(1,1)';
mu_hilim =  1.5*ones(1,1)';

mgp.lb = [v_lolim,w_lolim,f_lolim,g_lolim,Beta_lolim,mu_lolim]';
mgp.ub = [v_hilim,w_hilim,f_hilim,g_hilim,Beta_hilim,mu_hilim]';

% Training GP by MLE
[Model,OptHist] = trnGP(mgp,xtrn,ytrn);

% Making predictions
[Ypred,S2pred] = predGP(Model,X);

% Making plots
for i = 1:Nout
        figure(i)
        plot(X,Yori(:,i),'b--','linewidth',4);
        hold on
        plot(xtrn{i},ytrn{i},'bo','markerfacecolor','w','MarkerSize',12,'linewidth',4);
        plot(X,Ypred(:,i),':m','linewidth',4)
        plot(X,Ypred(:,i)+1.96*sqrt(S2pred(:,i)),'g:','linewidth',2);
        legend('Original function','Training points','GP prediction','95% confidence interval');
        plot(X,Ypred(:,i)-1.96*sqrt(S2pred(:,i)),'g:','linewidth',2);
        legend boxoff
        set(gca,'FontSize',25)
        box off
end