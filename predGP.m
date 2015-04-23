function [muYpred, varYpred] = predGP(Model, Xp)

N = Model.Ns;
xtrn_cell = Model.x;
ytrn_cell = Model.y;
Nin = size(cell2mat(xtrn_cell(1)),2);  % Number of input variables
Nout = size(ytrn_cell,2);              % Number of output variables

ytrn_num = cat(1,ytrn_cell(:));
y = cell2mat(ytrn_num);

%% Obtain the hyper parameters
C = Model.Cov;
Nin = Model.Nin;
Nout = Model.Nout;
hyp_est = Model.optsolu;
LowC = chol(C,'lower');
InvLowC = inv(LowC);
invC=InvLowC'*InvLowC;

for i = 1:Nout
        paras.v = hyp_est(1:Nout);
        paras.w = hyp_est(Nout+1:2*Nout);
        f = hyp_est(2*Nout+1:3*Nout);
        g = hyp_est(3*Nout+1:4*Nout);
        Beta=hyp_est(4*Nout+1:5*Nout);
        paras.mu = hyp_est(5*Nout+1);
end
paras.A = exp(f); paras.B = exp(g);
paras.sig = exp(Beta);
paras.p = Model.p;

%% Prediction
Npred = length(Xp);
for i=1: Npred  % At different training points
        for m = 1 : Nout         % For each output
            index  = 0;
            for n = 1 : Nout      % Cov with each variable 
                xs =  xtrn_cell{n};    
                    for j = 1: N(n) 
                         if n == 1
                             K(j,m) = Csub(paras,Xp(i),xs(j),m,n);
                         else
                             K(j+index,m) = Csub(paras,Xp(i),xs(j),m,n);
                         end
                    end
                          index = index+N(n);
            end
        end
        % Prediction
                muYpred(i,:) = (K'*invC*y)';
                kappa =paras.v.^2+paras.w.^2+paras.sig.^2;
                varYpred(i,:) = kappa-diag(K'*invC*K);       
end



