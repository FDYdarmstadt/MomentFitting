function [ A, b, X ] = assembleSurfaceSystem(boundarySegments, gradPhi, order)
%ASSEMBLESURFACESYSTEM Constructs the moment-fitting system for the volume
%case.
%   Constructs the matrix 'A' (M rows, N columns), the vector 'b' (M rows)
%   and the list of nodes 'X' (N rows, 2 columns) which define the
%   under-determined system 'A w = b' that can be used to determine the
%   quadrature weights 'w' at the nodes 'X' for evaluating surface
%   integrals over the zero set of the level set

    divBasis = getDivergenceFreeBasis(order);
    M = size(divBasis, 2);
    
    % Select nodes
    safetyFactor = 1.6;
    X = getNodes(ceil(sqrt(safetyFactor * M)));
    N = size(X, 1);
    
    % Assemble system matrix
    A = zeros(M, N);
    for i = 1:N
        x = X(i, 1);
        y = X(i, 2);
        
        [ normalX, normalY ] = getLevelSetNormal(gradPhi{1}(x,y), gradPhi{2}(x,y));
        
        for j = 1:M
            A(j,i) = divBasis{1, j}(x, y) * normalX + ...
                divBasis{2, j}(x, y) * normalY;
        end
    end
    
    % Assemble right-hand side
    b = zeros(M, 1);
    for j = 1:M
        % Top
        segment = boundarySegments(1,:);
        funTop = @(x)divBasis{2, j}(x, 1);
        intTop = integral(@(X)cellfun(funTop, num2cell(X)), segment(1), segment(2));
        
        % Bottom
        segment = boundarySegments(2,:);
        funBottom = @(x)-1*divBasis{2, j}(x, -1);
        intBottom = integral(@(X)cellfun(funBottom, num2cell(X)), segment(1), segment(2));
        
        % Right
        segment = boundarySegments(3,:);
        funRight = @(y)divBasis{1, j}(1, y);
        intRight = integral(@(X)cellfun(funRight, num2cell(X)), segment(1), segment(2));
        
        % Left
        segment = boundarySegments(4,:);
        funLeft = @(y)-1*divBasis{1, j}(-1, y);
        intLeft = integral(@(X)cellfun(funLeft, num2cell(X)), segment(1), segment(2));
        
        % Sum up
        b(j) = intTop + intBottom + intRight + intLeft;
    end
end

