function [plots, shps] = plot_gaussian_ellipsoid(X, S, sd, p)
% plot_gaussian_ellipsoid.m
% Benjamin Hanson, 2024
% 
% Given a set of 2D/3D state vectors X with associated weights P,
% generate an isosurface representing a curve of isovalue
% 
% Inputs:
%          X -- 2D/3D mean vector
%          S -- 2x2/3x3 covariance matrix
%         sd -- standard deviation of ellipse(s)
%          p -- plotting parameters (optional)
%               *   color -- isosurface color
%               * display -- handle visibility
%               *    name -- display name, if display==1
%               *   means -- plot weighted mean of point mass PDF
%               *     axh -- figure axis
%               *   alpha -- surface visibility
%               *     plt -- plotting boolean
% 
% Outputs:
%   plots -- plots for legends, with first being the mean

if ~exist('sd', 'var'), sd = 1; end
if ~exist('p','var')
    p.type = 'fill';
    p.display = 0; 
    for i = 1:numel(isovalue)
        p.color{i}=[1 0 0];
    end
    p.name = num2str(sd) + "\sigma covariance"; 
    p.means = 0; 
    p.axh=gca; 
    p.alpha=flip(logspace(log(0.3),log(0.6),numel(sd)));
    p.plt = 1; 
else
    if ~isfield(p,'type')
        p.type = 'fill';
    end
    if ~isfield(p,'display')
        p.display = 0;
    end
    if(p.display == 1)
        if ~isfield(p,'name')
            p.name = num2str(sd) + "\sigma covariance";
        end
    end
    if strcmp(p.type, 'line')
        if ~isfield(p,'color')
            p.color = 'r';
        end
        if ~isfield(p,'lw')
            p.lw = 1;
        end
        if ~isfield(p,'ls')
            p.ls = '-';
        end
    end
    if ~isfield(p,'means')
        p.means = 0;
    end
    if ~isfield(p,'color')
        for i = 1:numel(isovalue)
            p.color{i}=[1 0 0];
        end
    else
        if (isstring(p.color))||(ischar(p.color))||((all(size(p.color) == [1,3]))&&(~iscell(p.color)))
            col = p.color; p.color = {}; 
            for i = 1:numel(sd)
                p.color{i}=col;
            end
        end 
    end

    if ~isfield(p,'alpha')
        p.alpha=flip(logspace(log(0.3),log(0.6),numel(sd)));
    else
       p.alpha = p.alpha.*ones(1,numel(sd));
    end

    if ~isfield(p, 'axh'), p.axh = gca; end

    if ~isfield(p, 'plt'), p.plt = 1; end
end

if numel(X) ~= length(X), 
    error('M must be a vector'); 
end

if ~isfield(p, 'hgt')
    if ~( all(numel(X) == size(S)) )
        error('Dimensionality of M and S must match');
    end
end

if ~(isscalar(p.axh) && ishandle(p.axh) && strcmp(get(p.axh,'type'), 'axes'))
    error('Invalid axes handle');
end
p.alpha = sort(p.alpha, 'descend'); 
sd = sort(sd); 

set(p.axh, 'nextplot', 'add');

switch numel(X)
    case 2, [plots,shps] = plot_gaussian_ellipsoid2D(X,S,sd,p);
    case 3, [plots,shps] = plot_gaussian_ellipsoid3D(X,S,sd,p);
   otherwise
      error('Unsupported dimensionality');
end

if nargout==0,
    clear plots;
end

%-----------------------------
function [plots,shps] = plot_gaussian_ellipsoid2D(X, S, sd, p)

count = 1; 
for i=sd
    npts=50; 
    % plot the gaussian fits
    tt = linspace(0, 2*pi, npts);
    x = cos(tt); y=sin(tt);
    ap = [x(:) y(:)]';
    [v,d]=eig(S); 
    d = i * sqrt(d); % convert variance to sdwidth*sd
    bp = (v*d*ap) + repmat(X, 1, size(ap,2)); 
    if p.plt
        if p.means, scatter(p.axh, X(1), X(2), 100, p.color, "pentagram", 'filled','HandleVisibility','off'); end
        if(p.display==1)
            if(p.type=="fill")
                plots{count} = fill(p.axh, bp(1,:), bp(2,:), p.color{count}, 'EdgeColor', 'none', 'FaceAlpha', p.alpha(count),'DisplayName',p.name);
            elseif(p.type=="line")
                plots{count} = plot(p.axh, bp(1,:), bp(2,:), 'LineStyle', p.ls, 'Color', p.color{count}, 'LineWidth',1, 'DisplayName',p.name);
            end
        else
            if(p.type=="fill")
                plots{count} = fill(p.axh, bp(1,:), bp(2,:), p.color{count}, 'EdgeColor', 'none', 'FaceAlpha', p.alpha(count),'HandleVisibility','off');
            elseif(p.type=="line")
                plots{count} = plot(p.axh, bp(1,:), bp(2,:), 'LineStyle', p.ls, 'Color', p.color{count}, 'LineWidth',1, 'HandleVisibility','off');
            end
        end
    else
        plots{count} = NaN;
    end
    shps{count} = bp;
    p.display = 0; p.means = 0; count = count + 1; 
end

%-----------------------------
function [plots, shps] = plot_gaussian_ellipsoid3D(X, S, sd, p)

count = 1; 
for i=sd
    npts = 50; 
    if isfield(p, 'hgt')
        % plot the gaussian fits
        tt = linspace(0, 2*pi, npts);
        x = cos(tt); y=sin(tt);
        ap = [x(:) y(:)]';
        [v,d]=eig(S); 
        d = i * sqrt(d); % convert variance to sdwidth*sd
        bp = (v*d*ap) + repmat(X(2:3), 1, size(ap,2)); 
        xp = [X(1)*ones(size(bp(1,:))) + (length(sd)-count+1)*p.hgt X(1)*ones(size(bp(1,:))) + (length(sd)-count)*p.hgt]';
        yp = [bp(1,:) bp(1,:)]'; 
        zp = [bp(2,:) bp(2,:)]'; 
        [k, ~] = boundary(xp, yp, zp, 0.1); 

        if p.plt
            if p.means, scatter3(p.axh, X(1), X(2), X(3), 100, p.color, "pentagram", 'filled','HandleVisibility','off'); end
            if p.display
                plots{count} = trisurf(k, xp, yp, zp, 'Parent', p.axh, 'Edgecolor', 'none', 'FaceColor', p.color{count}, 'FaceAlpha', p.alpha(count), 'DisplayName', p.name);  
            else
                plots{count} = trisurf(k, xp, yp, zp, 'Parent', p.axh, 'Edgecolor', 'none', 'FaceColor', p.color{count}, 'FaceAlpha', p.alpha(count), 'HandleVisibility', 'off');
            end
        else
            plots{count} = NaN;
        end
        shps{count} = [reshape(xp,1,[])', reshape(yp,1,[])', reshape(zp,1,[])']; 
        p.display = 0; p.means = 0; count = count + 1; 
    else
        [x,y,z] = sphere(npts);
        ap = [x(:) y(:) z(:)]';
        [v,d]=eig(S); 
        if any(d(:) < 0)
           fprintf('warning: negative eigenvalues\n');
           d = max(d,0);
        end
        d = i * sqrt(d); % convert variance to sdwidth*sd
        bp = (v*d*ap) + repmat(X, 1, size(ap,2)); 
        xp = reshape(bp(1,:), size(x));
        yp = reshape(bp(2,:), size(y));
        zp = reshape(bp(3,:), size(z));
        if p.plt
            if p.means, scatter3(p.axh, X(1), X(2), X(3), 100, p.color, "pentagram", 'filled','HandleVisibility','off'); end
            if(p.display==1)
                plots{count} = surf(p.axh,xp,yp,zp,"FaceColor",p.color{count},"EdgeColor","none","FaceAlpha",p.alpha(count),'DisplayName',p.name);
            else
                plots{count} = surf(p.axh,xp,yp,zp,"FaceColor",p.color{count},"EdgeColor","none","FaceAlpha",p.alpha(count),'HandleVisibility','off');
            end
        else
            plots{count} = NaN;
        end
        shps{count} = [reshape(xp,1,[])', reshape(yp,1,[])', reshape(zp,1,[])']; 
        p.display = 0; p.means = 0; count = count + 1; 
    end
end