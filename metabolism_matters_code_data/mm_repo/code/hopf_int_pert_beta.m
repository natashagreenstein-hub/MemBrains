function [FC,CV,Cvth,A]=hopf_int_pert_beta(gC,f_diff,sigmafinal)
N=size(gC,1);
a=-0.01;
wo = f_diff'*(2*pi);

Cvth = zeros(2*N);

% Jacobian:

s = sum(gC,2);
B = diag(s);

aa=a*eye(N);

Axx = aa - B + gC;
Ayy = Axx;
Axy = -diag(wo);
Ayx = diag(wo);

A = [Axx Axy; Ayx Ayy];
sigmaf=sigmafinal;

sigmaf=[sigmaf sigmaf];
Qn = diag(sigmaf.^2);

Cvth=sylvester(A,A',-Qn);
Cvth=0.5*(Cvth+Cvth');
FCth=corrcov(Cvth);
FC=FCth(1:N,1:N);
CV=Cvth(1:N,1:N);
