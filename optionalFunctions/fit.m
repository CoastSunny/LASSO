function output = fit(U,D,V,Y,tau,TOL)
% Solves argmin_b ||Y-Xb||^2 subject to ||b|| <= tau
%X = UDV' is a SVD of X
%Y is the target, tau is constraint parameter

%compute relevant quantities
D=sparse(D);
Y = D'*U'*Y;
p = size(D,2);
D_sq = D'*D;

%initialize lambda
lambda = 1;
R = (D_sq+lambda*eye(p));
%%R = (D_sq+lambda*speye(p));
%Rv=diag(D_sq) + lambda;

%feasibility check
while lambda - f(Y,R,tau)/f_prime(Y,R) < 0
    lambda = lambda/2;
    R = (D_sq+lambda*speye(p));
    %R = (D_sq+lambda*eye(p));
end
    

%newton phase
tol = 1;
while tol > TOL
    lambda = lambda - f(Y,R,tau)/f_prime(Y,R);
    %R = (D_sq+lambda*eye(p));
     R = (D_sq+lambda*speye(p));
    tol = abs(f(Y,R,tau));
    
end

%output = V*((D_sq+lambda*eye(p))\Y);
output = V*(Y./diag(D_sq+ lambda));




%required functions------------------------------------
function value = f(Y,R,tau)
%value = Y'*(R^2\Y) - tau^2;
r=diag(R); r=max(r,10^-8).^2;
value = Y'*(Y./r) - tau^2;

function value = f_prime(Y,R)
%value = -2*Y'*(R^3\Y);
r=diag(R); r=max(r,10^-8).^3;
value = -2*Y'*(Y./r);





