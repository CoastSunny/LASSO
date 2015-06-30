function output = fit_modified_newton(U,D,V,Y,tau,TOL)
% Solves argmin_b ||Y-Xb||^2 subject to ||b|| <= tau
%X = UDV' is a SVD of X
%Y is the target, tau is constraint parameter

%compute relevant quantities

Y = D'*U'*Y;
D_sq = D'*D;

%initialize lambda

% search over eta=exp(lambda);


if f(Y,diag(D_sq),tau) <= 0
    output = V*(Y./diag(D_sq));    
else
    
lambda=10;
R = diag(D_sq) + lambda;
%newton phase
tol = 1;
iter=0;Max_newton=30;  back_track=0.8;
obj = zeros(1,Max_newton);

while ((tol > TOL)&&(iter<Max_newton))
    iter=iter+1;
    while lambda - f(Y,R,tau)/f_prime(Y,R) < 0
        lambda = lambda/2;
        R = diag(D_sq)+lambda;
    end
    change=f(Y,R,tau)/f_prime(Y,R);
    lambda = lambda - change;
    R = diag(D_sq) + lambda;
    
%     while (abs(f(Y,R_new,tau)) > abs(f(Y,R,tau))) || (lambdaN<0)
%     change=change*back_track;  
%     lambdaN = lambda - change;    R_new=diag(D_sq) + lambdaN;
%     end
     
    obj(iter)=f(Y,R,tau)/tau^2;
    tol = obj(iter);
end

%output = V*((D_sq+lambda*eye(p))\Y);
%output = V*(Y./diag(D_sq+ lambda));
output = V*(Y./R);


end

%required functions------------------------------------
function value = f(Y,R,tau)
%R=max(R,10^-8).^2;
R = R.^2;
value = Y'*(Y./R) - tau^2;

function value = f_prime(Y,R)
% R=max(R,10^-8).^3;
R = R.^3;
value = -2*Y'*(Y./R);





