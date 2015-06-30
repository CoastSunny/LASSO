function [beta, objValues] = get_solution(X, Y, betaInit, lambda, tol, maxIter)

    p = size(betaInit, 1);
    objValues = zeros(1, maxIter);
    activeIndex = ones(1, p);
    indicator = 0;
    iter = 0;

    beta = betaInit;
    res = Y - X*beta;
    pen = get_penalty(beta, lambda);
    obj = 0.5*norm(res, 'fro')^2 + pen;
    
    while (indicator == 0)
        tolerance = 1;
        while tolerance > tol && iter < maxIter
            objOld = obj;
            for i = find(activeIndex > 0)
                r = res + X(:,i) * beta(i,:);
                pen = pen - lambda * norm(beta(i,:));
                u = X(:,i)'*r;
                if norm(u) > lambda
                    activeIndex(i) = 1;
                    beta(i,:) = 1/norm(X(:,i))^2 * (1-lambda/norm(u)) * u;
                    res = r - X(:,i) * beta(i,:);
                    pen = pen + lambda * norm(beta(i,:));                
                else 
                    activeIndex(i) = 0;
                    beta(i,:) = 0 * beta(i,:);
                    res = r;
                end
                obj = 0.5*norm(res,'fro')^2 + pen;
                tolerance = (objOld - obj)/obj;
                iter = iter + 1;
                objValues(iter) = obj;
            end
        end
        %converged, so check if active set needs updating
        idx = zeros(1, p);
        for i = 1:p
            r = res + X(:,i) * beta(i,:);
            if norm(X(:,i)'*r) > lambda
                idx(i) = 1;
            end
        end
        if any(idx ~= activeIndex)
            activeIndex = idx;
        else 
            indicator = 1;
        end
    end
    objValues = objValues(objValues>0);
end

function penalty = get_penalty(beta, lambda)
    penalty = 0;
    p = size(beta, 1);
    for i = 1:p
        penalty = penalty + norm(beta(i,:));
    end
    penalty = penalty * lambda;
end