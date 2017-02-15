function [ divBasis ] = getDivergenceFreeBasis(order)
%GETDIVERGENCEFREEBASIS Returns a divergence-free vector basis of the
%specified 'order'
%   Returns a divergence-free vector basis of the specified 'order'. The
%   current implementation is based on monomials for simplicity, even
%   though this is a very bad choice in practice, especially for higher
%   orders. 
%   The result is given as a matrix (2 rows, M columns), where the first
%   index corresponds to the spatial coordinate.
    divBasis = {
        @(x,y)1         @(x,y)0
        @(x,y)0         @(x,y)1
        %
        @(x,y)y         @(x,y)0
        @(x,y)x         @(x,y)-y
        @(x,y)0         @(x,y)x
        %
        @(x,y)y^2       @(x,y)0
        @(x,y)x*y       @(x,y)-1/2*y^2
        @(x,y)x^2       @(x,y)-2*x*y
        @(x,y)0         @(x,y)x^2
        %
        @(x,y)y^3       @(x,y)0
        @(x,y)x*y^2     @(x,y)-1/3*y^3
        @(x,y)x^2*y     @(x,y)-x*y^2
        @(x,y)x^3       @(x,y)-3*x^2*y
        @(x,y)0         @(x,y)x^3
        %
        @(x,y)y^4       @(x,y)0
        @(x,y)x*y^3     @(x,y)-1/4*y^4
        @(x,y)x^2*y^2   @(x,y)-2/3*x*y^3
        @(x,y)x^3*y     @(x,y)-3/2*x^2*y^2
        @(x,y)x^4       @(x,y)-4*x^3*y
        @(x,y)0         @(x,y)x^4
        %
        @(x,y)y^5       @(x,y)0
        @(x,y)x*y^4     @(x,y)-1/5*y^5
        @(x,y)x^2*y^3   @(x,y)-1/2*x*y^4
        @(x,y)x^3*y^2   @(x,y)-x^2*y^3
        @(x,y)x^4*y     @(x,y)-2*x^3*y^2
        @(x,y)x^5       @(x,y)-5*x^4*y
        @(x,y)0         @(x,y)x^5
    }';

    M = (order + 1) * (order + 4) / 2;
    if size(divBasis, 2) < M
        error('Requested order not implemented');
    end
    
    divBasis = divBasis(:, 1:M);
end

