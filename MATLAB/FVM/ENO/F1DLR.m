function [p,S] = F1DLR(xbar,vbar,dx,x,i,k,r) 
% F1DLR.m
% Benjamin Hanson, 2024
%
% Following the framework from "Essentially Non-Oscillatory and Weighted
% Essentially  Non-Oscillatory Schemes for Hyperbolic Conservation Laws'' 
% by Chi-Wang Shu (1997), F1DUR stands for Fixed, 1D Lagrange Reconstruction. 
% This function performs a kth order accurate interpolation of v(x) using 
% the Lagrange form interpolation given by vbar, with a fixed stencil 
% determined by r
%
% Inputs:
%   xbar -- cell centers I_i, i=1,...,N
%   vbar -- cell average values of v(x) at cell I_i, i=1,...,N
%   dx   -- grid widths Delta x_i (or Delta x if uniform), i=1,...,N
%   x    -- x-value of interpolation
%   i    -- cell of interest where interpolation is being performed, i.e.
%           p_i(x)
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

S = [];
for j = i-r:i-r+k-1
    S(end+1) = j; 
end

p = 0;
for m=0:k
    temp_sum=0;
    for j=0:m-1
        temp_temp_sum=0; 
        for l=0:k
            if(l~=m)
                temp_prod=1; 
                for q=0:k
                    if(q~=m)&&(q~=l)
                        if i-r+q > N
                            xq = xbar(i-r+q-1)+(dx(i-r+q-1)/2);
                        else
                            xq = xbar(i-r+q)-(dx(i-r+q)/2);
                        end
                        temp_prod=temp_prod*(x-xq);
                    end
                end
                temp_temp_sum=temp_temp_sum+temp_prod;
            end
        end
        temp_prod=1;
        for l=0:k
            if(l~=m)
                if i-r+m > N
                    xm = xbar(i-r+m-1)+(dx(i-r+m-1)/2);
                else
                    xm = xbar(i-r+m)-(dx(i-r+m)/2);
                end
                if i-r+l > N
                    xl = xbar(i-r+l-1)+(dx(i-r+l-1)/2);
                else
                    xl = xbar(i-r+l)-(dx(i-r+l)/2);
                end
                temp_prod=temp_prod*(xm-xl);
            end
        end
        temp_sum=temp_sum+(vbar(i-r+j)*dx(i-r+j)*(temp_temp_sum/temp_prod));
    end
    p=p+temp_sum;
end
end 