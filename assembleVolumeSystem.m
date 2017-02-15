function [ A, b, X ] = assembleVolumeSystem(boundarySegments, gradPhi, order)
%ASSEMBLEVOLUMESYSTEM Constructs the moment-fitting system for the volume
%case.
%   Constructs the matrix 'A' (M rows, N columns), the vector 'b' (M rows)
%   and the list of nodes 'X' (N rows, 2 columns) which define the
%   under-determined system 'A w = b' that can be used to determine the
%   quadrature weights 'w' at the nodes 'X' for evaluating volume
%   integrals over the zero set of the level set

    [ basis, basisAntiderivative ] = getBasis(order);
    M = size(basis, 2);
    
    % Select nodes
    safetyFactor = 1.0;
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
    
    %
    % Assemble right-hand side
    %
    
    % Integral over sub-domain boundary without zero iso-contour
    intBoundary = zeros(M, 1);
    for j = 1:M
        % Top
        segment = boundarySegments(1,:);
        funTop = @(x)basisAntiderivative{2, j}(x, 1);
        intTop = integral(@(X)cellfun(funTop, num2cell(X)), segment(1), segment(2));
        
        % Bottom
        segment = boundarySegments(2,:);
        funBottom = @(x)-1*basisAntiderivative{2, j}(x, -1);
        intBottom = integral(@(X)cellfun(funBottom, num2cell(X)), segment(1), segment(2));
        
        % Right
        segment = boundarySegments(3,:);
        funRight = @(y)basisAntiderivative{1, j}(1, y);
        intRight = integral(@(X)cellfun(funRight, num2cell(X)), segment(1), segment(2));
        
        % Left
        segment = boundarySegments(4,:);
        funLeft = @(y)-1*basisAntiderivative{1, j}(-1, y);
        intLeft = integral(@(X)cellfun(funLeft, num2cell(X)), segment(1), segment(2));
        
        % Sum up
        intBoundary(j) = intTop + intBottom + intRight + intLeft;
    end
    
    % Integral over zero iso-contour by making use of the surface rule
    % Note that 'order + 1' is required to integrate the antiderivatives
    % of the basis functions with sufficient accuracy
    intSurface = zeros(M, 1);
    [ASurface, bSurface, nodesSurface] = assembleSurfaceSystem(boundarySegments, gradPhi, order + 1);
    origWarning = warning('query', 'MATLAB:rankDeficientMatrix');
    warning('off', 'MATLAB:rankDeficientMatrix');
    wSurface = ASurface\bSurface;
    warning(origWarning);
    for i = 1:size(nodesSurface, 1)
        x = nodesSurface(i, 1);
        y = nodesSurface(i, 2);

        [ normalX, normalY ] = getLevelSetNormal(gradPhi{1}(x,y), gradPhi{2}(x,y));
            
        for j = 1:M
            basisX = basisAntiderivative{1, j}(x, y);
            basisY = basisAntiderivative{2, j}(x, y);
            intSurface(j) = intSurface(j) + ...
                (basisX * normalX + basisY * normalY) * wSurface(i);
        end
    end
    
    b = -intSurface + intBoundary;
end

