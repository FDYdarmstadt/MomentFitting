function [ normalX, normalY ] = getLevelSetNormal(gradPhiX, gradPhiY)
%GETLEVELSETNORMAL Safely normalizes the given gradient
    gradNorm = sqrt(gradPhiX * gradPhiX + gradPhiY * gradPhiY);
        
    if (gradNorm < 1e-10)
        warningMessage = 'Norm of level set gradient is almost zero within domain of integration. Don''t expect nice convergence rates';
        warning(warningMessage);
    end

    % Make sure zero gradient doesn't cause NaNs
    if (gradNorm < 1e-10)
        gradNorm = 1; 
    end

    normalX = gradPhiX / gradNorm;
    normalY = gradPhiY / gradNorm;
end

