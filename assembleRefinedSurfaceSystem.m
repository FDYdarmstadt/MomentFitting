function [ A, b, X ] = assembleRefinedSurfaceSystem(boundarySegments, gradPhi, order, safetyFactor)
    [ basis, ~ ] = getBasis(order);
    M = size(basis, 2);
    
    % Select nodes
    if (~exist('safetyFactor', 'var'))
        safetyFactor = 1.0;
    end
    X = getNodes(ceil(sqrt(safetyFactor * M)));
    N = size(X, 1); 
    
    % Assemble system matrix 'A'
    A = zeros(M, N);
    for i = 1:N
        x = X(i, 1);
        y = X(i, 2);
        
        for j = 1:M
            A(j,i) = basis{j}(x, y);
        end
    end
    
    % Assemble right-hand side
    b = zeros(M, 1);
    [ADivFree, bDivFree, nodesDivFree] = assembleSurfaceSystem(boundarySegments, gradPhi, order);
    origWarning = warning('query', 'MATLAB:rankDeficientMatrix');
    warning('off', 'MATLAB:rankDeficientMatrix');
    wDivFree = ADivFree\bDivFree;
    warning(origWarning);
    for i = 1:size(nodesDivFree, 1)
        x = nodesDivFree(i, 1);
        y = nodesDivFree(i, 2);

        for j = 1:M
            b(j) = b(j) + basis{j}(x, y) * wDivFree(i);
        end
    end
end

