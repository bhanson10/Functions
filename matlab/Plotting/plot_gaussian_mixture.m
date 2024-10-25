function plot_gaussian_mixture(X,Ps,W,isovalue,p)
% plot_gaussian_mixture.m
% Benjamin Hanson, 2024
% 
% Given a set of 2D/3D state vectors X, ensemble covariances, and
% associated weights, generate an isosurface representing a curve of
% isovalues
% 
% Inputs:
%          X -- set of 2D/3D state vectors, MxN
%         Ps -- ensemble covariances of 2D/3D state vectorsn MxMxN
%          W -- weights of ensemble members, Nx1
%   isovalue -- isosurface value(s) to plot
%          p -- plotting parameters (optional)
%               *   color -- isosurface color
%               * display -- handle visibility
%               *    name -- display name, if display==1
%               *   means -- plot weighted mean of point mass PDF
%               *     axh -- figure axis
%               *   alpha -- surface visibility
%               *    type -- distribution type

% Checks and Balances
if length(X)~=length(Ps)
    error("Incongruous state vector/weight sets.")
end
if (max(isovalue)>1)||(min(isovalue)<0)
    error("Isovalue is outside of probability bounds [0,1].")
end
if ~exist('p','var')
    for i = 1:numel(isovalue)
        p.color{i}=[1 0 0];
    end
    p.display=0; 
    p.means=0;
    p.axh=gca; 
    p.alpha=flip(logspace(log(0.3),log(0.6),numel(isovalue)));
    p.type = "grid";
else
    if ~isfield(p,'color')
        for i = 1:numel(isovalue)
            p.color{i}=[1 0 0];
        end
    else
        if (isstring(p.color))||(ischar(p.color))||((all(size(p.color) == [1,3]))&&(~iscell(p.color)))
            col = p.color; p.color = {}; 
            for i = 1:numel(isovalue)
                p.color{i}=col;
            end
        end 
    end

    if ~isfield(p,'display')
        p.display=0;
    end
    if(p.display == 1)
        if ~isfield(p,'name')
            p.display=0;
        end
    end
    if isfield(p,'name')
        if ~isfield(p,'display')
            p.display=1;
        end
    end
    if ~isfield(p,'means')
        p.means=0;
    end
    if ~isfield(p,'axh')
        p.axh=gca;
    end
    if ~isfield(p,'alpha')
        p.alpha=flip(logspace(log(0.3),log(0.6),numel(isovalue)));
    else
        p.alpha = p.alpha.*ones(1,numel(isovalue));
    end

    if ~isfield(Ps, "type")
        p.type = "grid";
    end
end
p.alpha = sort(p.alpha); 
isovalue = sort(isovalue); 

% Getting number of state vectors
M=size(X, 1); 

switch M
    case 2, plot_gaussian_mixture2D(X,Ps,W,isovalue,p);
    case 3, plot_gaussian_mixture3D(X,Ps,W,isovalue,p);
   otherwise
      error('Unsupported dimensionality');
end

function plot_gaussian_mixture2D(X,Ps,W,isovalue,p)

% Get boundaries of domain 
M = 2; 
N = size(X, 2); 
Npts = 100; 
minmax = [inf.*ones([M,1]), -inf.*ones([M,1])];

for n = 1:N
    x = X(:, n); ps = Ps(:,:,n);
    bounds = bound_ellipse(x, ps, 3);
    for m = 1:M
        minmax(m, 1) = min(minmax(m, 1), min(bounds(m, 1))); 
        minmax(m, 2) = max(minmax(m, 2), max(bounds(m, 2)));
    end
end

% Generate linearly spaced points for each dimension
linspaces = cell(M, 1);
for m = 1:M
    linspaces{m} = linspace(minmax(m, 1), minmax(m, 2), Npts);
end
[mesh{1:M}] = ndgrid(linspaces{:});
chi = reshape(cat(M + 1, mesh{:}), [], M);

% Calculate probability at points
% P_full = zeros(size(chi, 1), 1); 
% K = N; 
% idx = knnsearch(X', chi, 'K', K);
% for i = 1:size(chi, 1)
%     flag_not_pd = 0; 
%     x = chi(i, :); 
%     as = zeros(size(idx(i,:),2), 1);
%     count = 1; 
%     for nm = idx(i,:)
%         mu = X(:, nm); ps = Ps(:, :, nm); w = W(nm);
%         try chol(ps);
%             R = chol(ps);
%             as(count) = (-0.5*sum((x' - mu).*(ps\(x' - mu)), 1) ...
%                 - sum(log(diag(R))) - (M/2)*log(2*pi) + log(w))';
%         catch ME
%             flag_not_pd = 1; 
%         end
%         count = count + 1; 
%     end
%     if flag_not_pd
%         P_full(i) = 0; 
%     else
%         ma = max(as);
%         lp = ma + log(sum(exp(as - ma),1));
%         P_full(i) = exp(lp);
%     end
% end
% P_full = P_full./max(P_full); 

% Calculate probability at points
flag_not_pd = 0; 
as = zeros(size(chi',2), N);
for n = 1:N
    x = X(:, n); ps = Ps(:, :, n); w = W(n); 
    try chol(ps);
        R = chol(ps);
        as(:, n) = (-0.5*sum((chi' - x).*(ps\(chi' - x)), 1) ...
            - sum(log(diag(R))) - (M/2)*log(2*pi) + log(w))';
    catch ME
        flag_not_pd = 1; 
    end
end
if flag_not_pd
    P_full = zeros(size(chi',2), 1); 
else
    ma = max(as, [], 2);
    lp = ma + log(sum(exp(as - ma), 2));
    lp = lp.';
    P_full = exp(lp)';
    P_full = P_full./max(P_full); 
end

P_full = reshape(P_full, Npts, Npts);
count = 1; 

for i=isovalue
    if (count == 1)&&p.display
        contour(p.axh, mesh{1}, mesh{2}, P_full, [i i], '-', 'EdgeAlpha', p.alpha(count), 'EdgeColor', 'none', 'FaceAlpha', p.alpha(count), 'FaceColor', p.color{count}, 'Fill', 'on', 'DisplayName', p.name);
    else
        contour(p.axh, mesh{1}, mesh{2}, P_full, [i i], '-', 'EdgeAlpha', p.alpha(count), 'EdgeColor', 'none', 'FaceAlpha', p.alpha(count), 'FaceColor', p.color{count}, 'Fill', 'on', 'HandleVisibility', 'off');
    end
    count = count + 1; 
end

if p.means, mean_X=X.*W; scatter(p.axh, mean(mean_X(:,1)), mean(mean_X(:,2)), 100, p.color{1}, "pentagram", "filled", 'HandleVisibility', 'off'); end

function plot_gaussian_mixture3D(X,Ps,W,isovalue,p)

if p.means, mean_X=sum(X.*W); scatter3(p.axh, mean_X(1), mean_X(2), mean_X(3), 100, p.color{1}, "pentagram", "filled", 'HandleVisibility', 'off'); end

% Get boundaries of domain 
M = 3; 
N = size(X, 2); 
Npts = 100; 
minmax = [inf.*ones([M,1]), -inf.*ones([M,1])];

for n = 1:N
    x = X(:, n); ps = Ps(:,:,n);
    bounds = bound_ellipse(x, ps, 3);
    for m = 1:M
        minmax(m, 1) = min(minmax(m, 1), min(bounds(m, 1))); 
        minmax(m, 2) = max(minmax(m, 2), max(bounds(m, 2)));
    end
end

% Generate linearly spaced points for each dimension
linspaces = cell(M, 1);
for m = 1:M
    linspaces{m} = linspace(minmax(m, 1), minmax(m, 2), Npts);
end
[mesh{1:M}] = ndgrid(linspaces{:});
chi = reshape(cat(M + 1, mesh{:}), [], M);

% Calculate probability at points
% P_full = zeros(size(chi, 1), 1); 
% K = N; 
% idx = knnsearch(X', chi, 'K', K);
% for i = 1:size(chi, 1)
%     flag_not_pd = 0; 
%     x = chi(i, :); 
%     as = zeros(size(idx(i,:),2), 1);
%     count = 1; 
%     for nm = idx(i,:)
%         mu = X(:, nm); ps = Ps(:, :, nm); w = W(nm);
%         try chol(ps);
%             R = chol(ps);
%             as(count) = (-0.5*sum((x' - mu).*(ps\(x' - mu)), 1) ...
%                 - sum(log(diag(R))) - (M/2)*log(2*pi) + log(w))';
%         catch ME
%             flag_not_pd = 1; 
%         end
%         count = count + 1; 
%     end
%     if flag_not_pd
%         P_full(i) = 0; 
%     else
%         ma = max(as);
%         lp = ma + log(sum(exp(as - ma),1));
%         P_full(i) = exp(lp);
%     end
% end
% P_full = P_full./max(P_full); 

as = zeros(size(chi',2), N);
for n = 1:N
    x = X(:, n); ps = Ps(:, :, n); w = W(n); 
    try chol(ps);
        R = chol(ps);
        as(:, n) = (-0.5*sum((chi' - x).*(ps\(chi' - x)), 1) ...
            - sum(log(diag(R))) - (M/2)*log(2*pi) + log(w))';
    catch ME
        flag_not_pd = 1; 
    end
end
if flag_not_pd
    P_full = zeros(size(chi',2), 1); 
else
    ma = max(as, [], 2);
    lp = ma + log(sum(exp(as - ma), 2));
    lp = lp.';
    P_full = exp(lp)';
    P_full = P_full./max(P_full); 
end

P_full = reshape(P_full, Npts, Npts, Npts);
count = 1; 
for i=isovalue
    if (count == 1)&&p.display
        patch(isosurface(mesh{1}, mesh{2}, mesh{3}, P_full, i), 'EdgeColor', 'none', 'FaceAlpha', p.alpha(count), 'FaceColor', p.color{count},'DisplayName', p.name); 
    else
        patch(isosurface(mesh{1}, mesh{2}, mesh{3}, P_full, i), 'EdgeColor', 'none', 'FaceAlpha', p.alpha(count), 'FaceColor', p.color{count},'HandleVisibility', 'off');
    end
    count = count + 1; 
end

function bounds = bound_ellipse(mu, sigma, n)
L = chol(sigma, 'lower');

d = size(sigma, 1);
num_points = 20^d; 
points = randn(d, num_points);
points = points ./ vecnorm(points);

scaled_points = n * L * points;
ellipse_points = mu + scaled_points;

bounds = zeros(d, 2);
for i=1:d
    bounds(i,1) = min(ellipse_points(i,:));
    bounds(i,2) = max(ellipse_points(i,:));
end