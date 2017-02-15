function [ boundarySegments ] = getBoundarySegments(phi)
%GETBOUNDARYSEGMENTS Returns the line segments the boundary of the cell
%[-1,1]x[-1,1] where 'phi' is positive
%   Determines the line segments that bound the sub-cell created by the
%   intersection of the reference cell [-1,1]x[-1,1] and the positive level
%   set region. The current implementation assumes that each edge of
%   [-1,1]x[-1,1] cannot be cut more than once, even though this is not a
%   requirement of the moment-fitting approach.
%   The result is returned as matrix (rows: edge index; columns; spatial
%   coordinate index), where the edges are ordered as follows:
%   - Row index 1: Top edge (i.e., point (-1, 1) to (1, 1))
%   - Row index 2: Bottom edge
%   - Row index 3: Right edge
%   - Row index 4: Left edge

    signTopLeft = sign(phi(-1, 1));
    signTopRight = sign(phi(1, 1));
    signBottomLeft = sign(phi(-1, -1));
    signBottomRight = sign(phi(1, -1));
    
    boundarySegments = zeros(4, 2);

    % Top
    if signTopLeft * signTopRight > 0
        % Uncut
        if signTopLeft > 0
            boundarySegments(1, :) = [-1 1];
        else
            boundarySegments(1, :) = [0 0];
        end
    else
        % Also catches case where zero is at corner
        xRoot = fzero(@(x)phi(x,1), [-1 1]);
        if signTopLeft > 0
            boundarySegments(1, :) = [-1 xRoot];
        else
            boundarySegments(1, :) = [xRoot 1];
        end
    end

    % Bottom
    if signBottomLeft * signBottomRight > 0
        % Uncut
        if signBottomLeft > 0
            boundarySegments(2, :) = [-1 1];
        else
            boundarySegments(2, :) = [0 0];
        end
    else
        % Also catches case where zero is at corner
        xRoot = fzero(@(x)phi(x,-1), [-1 1]);
        if signBottomLeft > 0
            boundarySegments(2, :) = [-1 xRoot];
        else
            boundarySegments(2, :) = [xRoot 1];
        end
    end

    % Right
    if signTopRight * signBottomRight > 0
        % Uncut
        if signTopRight > 0
            boundarySegments(3, :) = [-1 1];
        else
            boundarySegments(3, :) = [0 0];
        end
    else
        % Also catches case where zero is at corner
        yRoot = fzero(@(y)phi(1,y), [-1 1]);
        if signBottomRight > 0
            boundarySegments(3, :) = [-1 yRoot];
        else
            boundarySegments(3, :) = [yRoot 1];
        end
    end

    % Left
    if signTopLeft * signBottomLeft > 0
        % Uncut
        if signTopLeft > 0
            boundarySegments(4, :) = [-1 1];
        else
            boundarySegments(4, :) = [0 0];
        end
    else
        % Also catches case where zero is at corner
        yRoot = fzero(@(y)phi(-1,y), [-1 1]);
        if signBottomLeft > 0
            boundarySegments(4, :) = [-1 yRoot];
        else
            boundarySegments(4, :) = [yRoot 1];
        end
    end
end

