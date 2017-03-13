function output = fit_modified_newton_large(D_sq,U,sigma_sq,Y,tau,TOL)
% Solves argmin_b ||Y-Xb||^2 subject to ||b|| <= tau
%R = U*D_sq*U' is the eigendecomposition of the connectivity matrix
%sigma_sq is nonzero for the large groups
%Y is the target, tau is constraint parameter
if sigma_sq <= 0
    fprintf('%s\n','error');
    hey
end
%compute relevant quantities
normZ = norm(Y,2)^2 - norm(U'*Y,2)^2;
Y_comp = sqrt(D_sq).*(U'*Y);

if f(Y_comp,D_sq,tau,0,sigma_sq,normZ) <= 0
    output = Y;

else

%{
if f(Y_comp,D_sq,tau,1,sigma_sq,normZ) <= 0
    output = U'*(Y_comp./D_sq);    

%}

lambda=10;

%newton phase
tol = 1;
iter=0;Max_newton=30;  back_track=0.8;
obj = zeros(1,Max_newton);
while ((tol > TOL)&&(iter<Max_newton))
    iter=iter+1;
    while lambda-f(Y_comp,D_sq,tau,lambda,sigma_sq,normZ)/f_prime(Y_comp,D_sq,lambda,sigma_sq,normZ)<0
    lambda = lambda/2;
    end
    change=f(Y_comp,D_sq,tau,lambda,sigma_sq,normZ)/f_prime(Y_comp,D_sq,lambda,sigma_sq,normZ);

    lambdaN = lambda - change;
    
%     while (abs(f(Y_comp,D_sq,tau,lambdaN,sigma_sq,normZ)) > abs(f(Y_comp,D_sq,tau,lambda,sigma_sq,normZ))) || (lambdaN<0) 
%     change=change*back_track;  
%     lambdaN = lambda - change;
%     end

      
    obj(iter)=f(Y_comp,D_sq,tau,lambdaN,sigma_sq,normZ)/tau^2;
    %%obj(iter)
    tol = obj(iter);
    lambda = lambdaN;
end

output = U*((U'*Y)./(1+lambda*D_sq));
if sigma_sq > 0
    output = output + 1/(1+lambda*sigma_sq)*(Y-U*(U'*Y));
end


end

%required functions------------------------------------
function value = f(Y,D_sq,tau,lambda,sigma_sq,normZ)
value = Y'*(Y./(1+lambda*D_sq).^2) - tau^2;
if sigma_sq > 0
    value = value + sigma_sq/(1+lambda*sigma_sq)^2*normZ;
end

function value = f_prime(Y,D_sq,lambda,sigma_sq,normZ)
value = -2 * Y'*(D_sq.*Y./(1+lambda*D_sq).^3);
if sigma_sq > 0
    value = value - 2*normZ*sigma_sq^2/(1+lambda*sigma_sq)^3;
end





