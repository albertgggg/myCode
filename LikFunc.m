function L = LikFunc(x,varargin)
global counter
counter = counter+1;
if mod(counter,100) == 0
    counter
end
% 
xtrn_cell = varargin{:,1};
ytrn_cell = varargin{:,2};
Nin = size(cell2mat(xtrn_cell(1)),2);  % Number of input variables
Nout = size(ytrn_cell,2);              % Number of output variables

% Parameter assignment
paras.v = x(1:Nout);
paras.w= x(Nout+1:2*Nout);
f = x(2*Nout+1:3*Nout);
g = x(3*Nout+1:4*Nout);
Beta = x(4*Nout+1:5*Nout);
paras.mu = x(5*Nout+1);

paras.A = exp(f);
paras.B = exp(g);
paras.p = 1;
paras.sig = exp(Beta);

for i = 1 : Nout
     N(i) = length(ytrn_cell{i});
end
%% C (correlation coefficient matrix)
%% Try vectorized programming for Csub
for m = 1:Nout
    %% Construct diagonal entries
    xs = xtrn_cell{m};
    for i = 1:N(m)
        for j = 1:N(m)
            Cii(i,j) = Csub(paras,xs(i),xs(j),m,m);
        end
    end
    
    if m>1  % Not First row/column
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


%% Likelihood function
index = 0;
for i = 1:Nout
    y(index+1:index + N(i),1) = ytrn_cell{i};
    index = index + N(i);
end

LowC = chol(C,'lower');
InvLowC = inv(LowC);
invC=InvLowC'*InvLowC;
invC=inv(C);
L =-(-1/2*log(det(C)) - 1/2*y'*invC*y - sum(N)/2*log(2*pi));

