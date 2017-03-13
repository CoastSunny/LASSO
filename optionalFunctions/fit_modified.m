function output = fit_modified(U,D,V,Y,tau,TOL)
% Solves argmin_b ||Y-Xb||^2 subject to ||b|| <= tau
%X = UDV' is a SVD of X
%Y is the target, tau is constraint parameter

%compute relevant quantities

D=sparse(D);
Y = D'*U'*Y;
p = size(D,2);
D_sq = D'*D;

%initialize lambda

% search over eta=exp(lambda);


if f(Y,diag(D_sq),tau,1) <= 0

output = V*(Y./diag(D_sq));    
 
else
    
lambda=10;

R=diag(D_sq) + lambda;

%newton phase
tol = 1;
iter=0;Max_newton=30;  back_track=0.8;

while ((tol > TOL)&&(iter<Max_newton))

%    lambda = lambda - f(Y,R,tau)/f_prime(Y,R);
    iter=iter+1;
    R=diag(D_sq) + lambda;
    change=f(Y,R,tau,lambda)/f_prime(Y,R,lambda);

    lambdaN = lambda - change;    R=diag(D_sq) + lambdaN;
    
    while (abs(f(Y,R,tau,lambdaN)) > abs(f(Y,R,tau,lambda))) || (lambdaN<0)
    change=change*back_track;  
    lambdaN = lambda - change;    R=diag(D_sq) + lambdaN;
    end
    
    lambda=lambdaN;   
   
    obj(iter)=f(Y,R,tau,lambda);
    %%obj(iter)
end

%output = V*((D_sq+lambda*eye(p))\Y);
output = V*(Y./diag(D_sq+ lambda));


end

%required functions------------------------------------
function value = f(Y,R,tau,lambda)
%value = Y'*(R^2\Y) - tau^2;
R=max(R,10^-8).^2;
value = Y'*(Y./R) - tau^2;

function value = f_prime(Y,R,lambda)
%value = -2*Y'*(R^3\Y);
R=max(R,10^-8).^3;
value = -2*Y'*(Y./R);





