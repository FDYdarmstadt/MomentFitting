function [ nodes, weights ] = getVolumeRule(phi, gradPhi, order)
%GETVOLUMERULE Constructs a volume quadrature rule
%   Constructs a quadrature rule to integrate over the sub-domain of the
%   cell the cell [-1,1]x[-1,1] where the level set function 'phi' is
%   positive. The rule is constructed using moments up to the specified
%   'order'.

    boundarySegments = getBoundarySegments(phi);
    [A, b, nodes] = assembleVolumeSystem(boundarySegments, gradPhi, order);
    
    % Solves the under-determined system
    % Note that 'A' may be severely ill-conditioned for higher orders and
    % challenging zero-contours, so the backslash operator might not be the
    % best choice in call cases
    origWarning = warning('query', 'MATLAB:rankDeficientMatrix');
    warning('off', 'MATLAB:rankDeficientMatrix');
    weights = A\b;
    warning(origWarning);
end