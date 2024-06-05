function V = ND1D(xbar,vbar,dx,k)
% ND1D.m
% Benjamin Hanson, 2024
%
% Following the framework from "Essentially Non-Oscillatory and Weighted
% Essentially  Non-Oscillatory Schemes for Hyperbolic Conservation Laws'' 
% by Chi-Wang vShu (1997), ND1D stands for Newton Differences, 1D. 
% This function solves for the undivided/divided difference matrix V(N,k)
%
% Inputs:
%   xbar -- cell centers I_i, i=1,...,N
%   vbar -- cell average values of v(x) at cell I_i, i=1,...,N
%   dx   -- grid widths Delta x_i (or Delta x if uniform), i=1,...,N
%   k    -- order of accuracy of interpolation
%   d    -- Difference type
%           * 'u' for uniform divided difference
%           * 'n' for nonuniform divided difference
%
% Outputs:
%   V    -- Divided difference matrix of size N x k, where N is the length 
%           of xbar/vbar. The value at V(i,m) represents the m-th degree
%           divided different of the left-most point x_i

d = 'n'; % Assume nonuniform divided difference
N = length(xbar); 
V = NaN(N,k);
if length(dx)==1
    dx = dx.*ones(N,1); 
    d = 'u'; % Uniform divided difference
end

for i=1:N
    x = [xbar(i)-(dx(i)/2) xbar(i)+(dx(i)/2)]; 
    v = [vbar(i)]; 

    V(i,1) = vbar(i); 
    for m=1:k-1
        if i+m <= N
            x(end+1) = xbar(i+m)+(dx(i+m)/2);
            v(end+1) = vbar(i+m);
            if strcmp(d,'u')
                V(i,m+1) = UDD1D(x, v); 
            else
                V(i,m+1) = DD1D(x, v); 
            end
        end
    end
end
end

function V = DD1D(x,v) % Divided Difference, 1D
    if length(v)==1
        V = v(1);
    else
        V = (DD1D(x(2:end),v(2:end))-DD1D(x(1:end-1),v(1:end-1)))/(x(end)-x(1));
    end
end

function V = UDD1D(x,v) % Undivided Difference, 1D
    if length(v)==1
        V = v(1);
    else
        V = (UDD1D(x(2:end),v(2:end))-UDD1D(x(1:end-1),v(1:end-1)));
    end
end

