function [p,S] = F2DNR(xbar,V,dx,x,i,k,r) 
% F2DNR.m
% Benjamin Hanson, 2024
%
% Following the framework from "Essentially Non-Oscillatory and Weighted
% Essentially  Non-Oscillatory Schemes for Hyperbolic Conservation Laws'' 
% by Chi-Wang Shu (1997), F1DUR stands for Fixed, 2D Newtonian Reconstruction. 
% This function performs a kth order accurate interpolation of v(x) using 
% the Newtonian Difference interpolation given by V, with a fixed stencil 
% determined by r
%
% Inputs:
%   xbar -- cell centers I_ij, i=1,...,Nx, j=1,...,Ny
%   V    -- Divided difference matrix of size N x k, where N is the length 
%           of xbar/vbar. The value at V(i,m) represents the m-th degree
%           divided different of the left-most point x_i
%   dx   -- grid widths Delta x_ij (or Delta x if uniform), i=1,...,Nx, 
%           j=1,...,Ny
%   x    -- (x-value,y-value) of interpolation
%   i    -- cell of interest where interpolation is being performed, i.e.
%           p_ij(x)
%   k    -- order of accuracy of interpolation, generating a k-1th
%           degree polynomial p(x)
%   r    -- number of cells to the left used in stencil S(i)
%
% Outputs:
%   p    -- interpolated value of v at x
%   S    -- stencil used at I_i

d = 'n'; % Assume nonuniform divided difference
N = length(xbar); 
if length(dx)==1
    dx = dx.*ones(N,1); 
    d = 'u'; % Uniform divided difference
end

n = length(x); 
p = NaN(n,1); S = {};
for q=1:n
    iq = i(q); rq = r(q); 
    Sq = [];
    for j = iq-rq:iq-rq+k-1
        Sq(end+1) = j; 
    end
    S{q} = Sq; 

    pq = 0; 
    for j=1:k
        if strcmp(d,'u')
            v = V(iq-rq,j)/(factorial(j)*dx(1)^(j-1));
        else
            v = V(iq-rq,j);
        end
        temp_sum = 0;
        for m=0:j-1
            temp_prod = 1;
            for l=0:j-1
                if(l~=m)
                    temp_prod = temp_prod*(x(q)-(xbar(iq-rq+l)-(dx(iq-rq+l)/2)));
                end
            end
            temp_sum = temp_sum + temp_prod;
        end
        pq = pq + v*temp_sum;
    end
    p(q)=pq;
end
if n==1, S = S{1}; end
end 