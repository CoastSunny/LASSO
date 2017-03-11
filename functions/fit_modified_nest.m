function [bnew,obj_nest] = fit_modified_nest(D,Dtrans,Lipsh,Y,tau,TOL,bold)
% Solves argmin_b 0.5*||Y-Xb||^2 subject to ||b|| <= tau
%X = UDV' is a SVD of X
%Y is the target, tau is constraint parameter

%compute relevant quantities

% X=U*D*V';

% Lipsh=D(1,1)^2;

iter=1;

if isempty(bold)
bold=zeros(size(Dtrans,2),1);
end

bcouple=bold; bnew=bold;
tk=1;
Max_nest=1000;

obj_nest=zeros(Max_nest,1);

obj_nest(iter)=norm(Y-Dtrans*bnew,2)^2;
tol=1000;

while ((tol > TOL)&&(iter<Max_nest))

%   lambda = lambda - f(Y,R,tau)/f_prime(Y,R);
    iter=iter+1;

grad_vec= -D*(Y-Dtrans*bcouple);

bnew= bcouple - grad_vec/Lipsh;
nrm_bnew=norm(bnew,2);

if (nrm_bnew > tau)
bnew=tau*bnew/nrm_bnew;
end

tk1= 1 + sqrt(1 + 4*tk^2); tk1=tk1/2;

bcouple= bnew + ((tk- 1)*(bnew-bold)/tk1) ;

tol = Lipsh*norm(bold- bnew,2);

bold=bnew;

obj_nest(iter)=norm(Y-Dtrans*bnew,2)^2;

tk=tk1;


end


obj_nest=obj_nest(1:iter);



