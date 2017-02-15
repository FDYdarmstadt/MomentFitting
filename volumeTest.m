function tests = volumeTest
    tests = functiontests(localfunctions);
end

function testStraightVertical(~)
    phi = @(x,y)x;
    gradPhi = {@(x,y)1 @(x,y)0};
    exactVolume = 2;

    [~, w] = getVolumeRule(phi, gradPhi, 2);

    error = abs(sum(w) - exactVolume);
    assert(error < 1e-14);
end

function testStraightHorizontal(~)
    phi = @(x,y)y;
    gradPhi = {@(x,y)0 @(x,y)1};
    exactVolume = 2;

    [~, w] = getVolumeRule(phi, gradPhi, 2);

    error = abs(sum(w) - exactVolume);
    assert(error < 1e-14);
end

function testStraightVerticalShifted(~)
    phi = @(x,y)x - 0.5;
    gradPhi = {@(x,y)1 @(x,y)0};
    exactVolume = 1;

    [~, w] = getVolumeRule(phi, gradPhi, 2);

    error = abs(sum(w) - exactVolume);
    assert(error < 1e-14);
end

function testStraightDiagonal(~)
    phi = @(x,y)x + y;
    gradPhi = {@(x,y)1 @(x,y)1};
    exactVolume = 2;

    [~, w] = getVolumeRule(phi, gradPhi, 2);
    
    error = abs(sum(w) - exactVolume);
    assert(error < 1e-14);
end

function testQuarterCircle(~)
    phi = @(x,y)(x-1.5)^2 + (y-1.5)^2 - 5/2;
    gradPhi = {@(x,y)2*(x-1.5) @(x,y)2*(y-1.5)};
    exactVolume = 3.34088097749782;
    
    errors = zeros(1, 4);
    for i=1:3
        [~, w] = getVolumeRule(phi, gradPhi, i);
        errors(i) = abs(sum(w) - exactVolume);
    end
    
    % Order 4 requires higher safetyFactor (for equidistant nodes)
    [~, w] = getVolumeRule(phi, gradPhi, i, 2);
    errors(i) = abs(sum(w) - exactVolume);
    
    thresholds = [ 0.17 0.06 0.003 0.002 ];
    assert(all(errors < thresholds));
    
    % For reference: Errors obtained with an orthonormal basis, a
    % divergence-free basis based on this orthonormal basis, Gauss nodes
    % and complete orthogonal factorization:
    % optimizedResults = [ 0.00520 0.00114 0.00040 0.00016 ];
    % Note that these results where obtained _without_ using the increased
    % safety factor used in line 62
end