function [ basis, basisAntiderivative ] = getBasis(order)
%GETBASIS Returns a polynomial basis of the specified 'order'
%   Returns a basis of the polynomial space of the specified 'order' and
%   the corresponding 'antiderivatives' in the sense that the divergence of
%   the antiderivative yields the corresponding basis function. The current
%   implementation uses monomials for simplicity, even though this is a
%   very bad choice in practice, especially for higher orders
%   The result consists of a row-vector 'basis' (1 row, M columns) and the
%   matrix 'basisAntiderivative' (2 rows, M columns). For the latter, the
%   column index corresponds to the spatial coordinate.
    basis = {
        @(x,y)1
        %
        @(x,y)y
        @(x,y)x
        %
        @(x,y)y^2
        @(x,y)x*y
        @(x,y)x^2
        %
        @(x,y)y^3
        @(x,y)x*y^2
        @(x,y)x^2*y
        @(x,y)x^3
        %
        @(x,y)y^4
        @(x,y)x*y^3
        @(x,y)x^2*y^2
        @(x,y)x^3*y
        @(x,y)x^4
    }';
    basisAntiderivative = {
        @(x,y)1/2*x         @(x,y)1/2*y
        %
        @(x,y)1/2*x*y       @(x,y)1/4*y^2
        @(x,y)1/4*x^2       @(x,y)1/2*x*y
        %
        @(x,y)1/2*x*y^2     @(x,y)1/6*y^3
        @(x,y)1/2*x^2*y     @(x,y)1/2*x*y^2
        @(x,y)1/6*x^3       @(x,y)1/2*x^2*y
        %
        @(x,y)1/2*x*y^3     @(x,y)1/8*y^4
        @(x,y)1/4*x^2*y^2   @(x,y)1/6*x*y^3
        @(x,y)1/6*x^3*y     @(x,y)1/4*x^2*y^2
        @(x,y)1/8*x^4       @(x,y)1/2*x^3*y
        %
        @(x,y)1/2*x*y^4     @(x,y)1/10*y^5
        @(x,y)1/4*x^2*y^3   @(x,y)1/8*x*y^4
        @(x,y)1/6*x^3*y^2   @(x,y)1/6*x^2*y^3
        @(x,y)1/8*x^4*y     @(x,y)1/4*x^3*y^2
        @(x,y)1/10*x^5      @(x,y)1/2*x^4*y
    }';
    
    M = (order + 1) * (order + 2) / 2;
    if size(basis, 2) < M
        error('Requested order not implemented');
    end
    
    basis = basis(1, 1:M);
    basisAntiderivative = basisAntiderivative(:, 1:M);
end

