function [ nodes, weights ] = getSurfaceRule(phi, gradPhi, order, safetyFactor)
%GETSURFACERULE Constructs a surface quadrature rule
%   Constructs a quadrature rule to integrate over the zero-set of the
%   level set function 'phi' within the cell [-1,1]x[-1,1] using moments up
%   the specified 'order' within the moment-fitting approach.
    
    boundarySegments = getBoundarySegments(phi);
    if (exist('safetyFactor', 'var'))
        [A, b, nodes] = assembleSurfaceSystem(boundarySegments, gradPhi, order, safetyFactor);
    else
        [A, b, nodes] = assembleSurfaceSystem(boundarySegments, gradPhi, order);
    end
    
    % Solves the under-determined system
    % Note that 'A' may be severely ill-conditioned for higher orders and
    % challenging zero-contours, so the backslash operator might not be the
    % best choice in call cases
    origWarning = warning('query', 'MATLAB:rankDeficientMatrix');
    warning('off', 'MATLAB:rankDeficientMatrix');
    weights = A\b;
    warning(origWarning);
end