function [Model,hist] = trnGP(bounds,xtrn_cell,ytrn_cell)

problem.f = @LikFunc; % Likelihood function for the hyper parameters
mgp = bounds;

%% Setup parameters for DIRECT.
opts.ep        = 1e-4; %Jones factor                      (default is 1e-4)
opts.maxevals  =0.5e3; %max. number of function evals     (default is 20)
opts.maxits    = 5e1; %max. number of iterations         (default is 10)
opts.maxdeep   = 1e3;  %max. number of rect. divisions    (default is 100)
opts.testflag  = 0;    %if globalmin known, 0 otherwise (default is 0)
opts.showits   = 0;    %1 if disp. stats shown, 0 oth.
[fval,hyp_est,hist] = Direct(problem,[mgp.lb mgp.ub],opts,xtrn_cell,ytrn_cell); % Maximum likelihood estimation of the parameters

Nin = size(cell2mat(xtrn_cell(1)),2);  % Number of input variables
Nout = size(ytrn_cell,2);              % Number of output variables

for i = 1:Nout
        paras.v = hyp_est(1:Nout);
        paras.w = hyp_est(Nout+1:2*Nout);
        f = hyp_est(2*Nout+1:3*Nout);
        g = hyp_est(3*Nout+1:4*Nout);
        Beta=hyp_est(4*Nout+1:5*Nout);
        paras.mu = hyp_est(5*Nout+1);
        paras.A(i) = exp(f(i)); paras.B(i) = exp(g(i));
        paras.sig(i) = exp(Beta(i));
end

paras.p = Nin;

for i = 1 : Nout
     N(i) = length(ytrn_cell{i});
end

%% Calculate the correlation coefficient matrix with optimal hyper parameters
%% Vectorized programming for Csub
for m = 1:Nout
    %% Construct diagonal entries
    xs = xtrn_cell{m};
    for i = 1:N(m)
        for j = 1:N(m)
            Cii(i,j) = Csub(paras,xs(i),xs(j),m,m);
        end
    end
    
    if m>1  % If not the first row/column
       rownum = 0;
      
       for i = 1:(m-1)
          rownum = rownum + N(i);
       end
       
    else
        rownum = 0; %First row/column
    end
    
        C(rownum+1:rownum+N(m),rownum+1:rownum+N(m))=Cii;
              clearvars Cii
    
%% Construct non-diagonal entries
    for n = 1:Nout
        
        if m~=n
        
            xsm = xtrn_cell{m};
            xsn = xtrn_cell{n};
            
            for i = 1 : N(m)
                for j = 1 : N(n)
                     Cmn(i,j) = Csub(paras,xsm(i),xsn(j),m,n);
                     Cnm(j,i) = Csub(paras,xsn(j),xsm(i),n,m);
                end
            end
   %% starting point
   
   if m>1
   
       rownum = 0;
       
       for i = 1:(m-1)
           rownum = rownum + N(i);
       end
       
   else
       
       rownum = 0;
   
   end
   
   if n>1
   
       colnum = 0;
    
       for i = 1:(n-1)
        colnum = colnum + N(i);
       end
       
   else
       
       colnum = 0;
   
   end

    C(rownum+1:rownum+N(m),colnum+1:colnum+N(n))=Cmn;
    C(colnum+1:colnum+N(n),rownum+1:rownum+N(m))=Cnm;
    clearvars Cmn Cnm
        end
    end
end


%% Assembly the parameters in 'Model' class
Model.Cov = C;
Model.Ns = N;
Model.Nin = Nin; Model.Nout = Nout;
Model.x = xtrn_cell;
Model.y = ytrn_cell;
Model.p = paras.p;
Model.optsolu = hyp_est;

