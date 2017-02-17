function [ nodes, weights ] = getRefinedSurfaceRule(phi, gradPhi, order, safetyFactor)
    boundarySegments = getBoundarySegments(phi);
    if (exist('safetyFactor', 'var'))
        [A, b, nodes] = assembleRefinedSurfaceSystem(boundarySegments, gradPhi, order, safetyFactor);
    else
        [A, b, nodes] = assembleRefinedSurfaceSystem(boundarySegments, gradPhi, order);
    end
    
    origWarning = warning('query', 'MATLAB:rankDeficientMatrix');
    warning('off', 'MATLAB:rankDeficientMatrix');
    weights = A\b;
    warning(origWarning);
end