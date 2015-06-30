function output = fit_mike(U,D,V,Y,tau,TOL)

%X = UDV' is a SVD of X
%Y is the target, tau is constraint parameter

%compute relevant quantities
Y = D'*U'*Y;
p = size(D,2);
D_sq = D'*D;

%initialize lambda
lambda = 10;
R = (D_sq+lambda*eye(p));

%feasibility check
while lambda - f(Y,R,tau)/f_prime(Y,R) < 0
    lambda = lambda/2;
    R = (D_sq+lambda*eye(p));
    lambda - f(Y,R,tau)/f_prime(Y,R) 
end
    

%newton phase
tol = 1;
while tol > TOL
    lambda = lambda - f(Y,R,tau)/f_prime(Y,R);
    R = (D_sq+lambda*eye(p));
    tol = abs(f(Y,R,tau))
end

output = V*((D_sq+lambda*eye(p))\Y);

%required functions------------------------------------
function value = f(Y,R,tau)
value = Y'*(R^2\Y) - tau^2;

function value = f_prime(Y,R)
value = -2*Y'*(R^3\Y);