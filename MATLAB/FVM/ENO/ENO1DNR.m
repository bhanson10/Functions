function [p,S] = ENO1DNR(xbar,V,dx,x,i,k)
% ENO1DNR.m
% Benjamin Hanson, 2024
%
% Following the framework from "Essentially Non-Oscillatory and Weighted
% Essentially  Non-Oscillatory Schemes for Hyperbolic Conservation Laws'' 
% by Chi-Wang Shu (1997), ENO1DNR stands for ENO, 1D Newtonian Reconstruction. 
% This function performs a kth order accurate interpolation of v(x) using 
% the Newtonian Difference interpolation given by V, with an adaptive stencil 
% determined by V
%
% Inputs:
%   xbar -- cell centers I_i, i=1,...,N
%   V    -- Divided difference matrix of size N x k, where N is the length 
%           of xbar/vbar. The value at V(i,m) represents the m-th degree
%           divided different of the left-most point x_i
%   dx   -- grid widths Delta x_i (or Delta x if uniform), i=1,...,N
%   x    -- x-value of interpolation
%   i    -- cell of interest where interpolation is being performed, i.e.
%           p_i(x)
%   k    -- order of accuracy of interpolation, generating a k-1th
%           degree polynomial p(x)
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
    iq = i(q);
    Sq = [iq]; 
    r = 0; s = 0; %r: # of cells left, s: # of cells right
    for j=1:k-1
        if iq-r-1 < 1 
            s = s + 1; Sq = [Sq iq+s]; 
        elseif isnan(V(iq-r,j+1))
            r = r + 1; Sq = [iq-r Sq];
        else
            if abs(V(iq-r-1,j+1)) < abs(V(iq-r,j+1))
                r = r + 1; Sq = [iq-r Sq]; 
            else
                s = s + 1; Sq = [Sq iq+s]; 
            end
        end
    end
    S{q} = Sq;

    pq = 0; 
    for j=1:k
        if strcmp(d,'u')
            v = V(iq-r,j)/(factorial(j)*dx(1)^(j-1));
        else
            v = V(iq-r,j);
        end
        temp_sum = 0;
        for m=0:j-1
            temp_prod = 1;
            for l=0:j-1
                if(l~=m)
                    temp_prod = temp_prod*(x(q)-(xbar(iq-r+l)-(dx(iq-r+l)/2)));
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