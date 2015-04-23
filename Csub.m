%% Try vectorized programming
function Cyij = Csub(paras,si,sj,i,j)
    A = paras.A;
    B = paras.B;
    p = paras.p;
    v = paras.v;
    w = paras.w;
    Sig = paras.sig;
    Sigma = A(i)*A(j)/(A(i)+A(j)); % For 1 D input case
    mu = paras.mu;
 d = si - sj;
 if d == 0
     delta = 1;
 else
     delta = 0;
 end

if i == j
    Cuij =  pi^(p/2)*v(i)^2/sqrt(det(A(i)))*exp(-1/4*d*A(i)*d);
    Cvij =  pi^(p/2)*w(i)^2/sqrt(det(B(i)))*exp(-1/4*d*B(i)*d);
    Cyij = Cuij + Cvij + delta*Sig(i)^2;
else
    Cuij =  (2*pi)^(p/2)*v(i)*v(j)/sqrt(det(A(i)+A(j)))*exp(-1/2*(d-mu)'*Sigma*(d-mu));
    Cyij = Cuij;
end
