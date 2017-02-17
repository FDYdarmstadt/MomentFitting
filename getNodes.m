function [ X ] = getNodes(nodesPerDirection)
%GETNODES Nodes for a tensor product quadrate rule
%   Creates a set of nodes from the tensor product of one-dimensional node
%   distributions with 'nodesPerDirection' nodes, i.e. the result contains
%   'nodesPerDirection'^2 nodes. For simplicity, equidistant nodes are used
%   here, even though Gauss nodes typically deliver better results.

    %nodes1D = linspace(-1, 1, nodesPerDirection);
    nodes1D = gaussNodes(nodesPerDirection);
    N = nodesPerDirection^2;
    
    k = 1;
    X = zeros(N, 2);
    for i = 1:nodesPerDirection
        for j = 1:nodesPerDirection
            X(k, 1) = nodes1D(i);
            X(k, 2) = nodes1D(j);
            k = k + 1;
        end
    end
end