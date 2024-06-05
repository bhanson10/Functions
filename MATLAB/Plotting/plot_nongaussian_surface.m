function [plots, shps] = plot_nongaussian_surface(X,P,isovalue,alpha,p)
% plot_nongaussian_surface.m
% Benjamin Hanson, 2024
% 
% Given a set of 2D/3D state vectors X with associated weights P,
% generate an isosurface representing a curve of isovalue
% 
% Inputs:
%          X -- set of 2D/3D state vectors
%          P -- weights of 2D/3D state vectors
%   isovalue -- isosurface value(s) to plot
%      alpha -- shrinkage parameter from [0 (convex hull), 1] (optional)
%               (default = 0.5)
%          p -- plotting parameters (optional)
%               *   color -- isosurface color
%               * display -- handle visibility
%               *    name -- display name, if display==1
%               *   means -- plot weighted mean of point mass PDF
%               *     axh -- figure axis
%               *   alpha -- surface visibility
%               *     plt -- plotting boolean
%               *    type -- Contour plot type
%                            +  'sharp' -- plots true alpha convex hull
%                                          (default)
%                            + 'smooth' -- uses 3D smoothing technique to
%                                          smooth alpha convex hull 
% 
% Outputs:
%   plots -- plots for legends, with first being the mean
%    shps -- alphaShapes at isovalues

% Checks and Balances
if length(X)~=length(P)
    error("Incongruous state vector/weight sets.")
end
if (max(isovalue)>1)||(min(isovalue)<0)
    error("Isovalue is outside of probability bounds [0,1].")
end
if ~exist('alpha','var')
    alpha=0.5;
end
if ~exist('p','var')
    for i = 1:numel(isovalue)
        p.color{i}=[1 0 0];
    end
    p.display=0; 
    p.means=0;
    p.axh=gca; 
    p.alpha=flip(logspace(log(0.3),log(0.6),numel(isovalue)));
    p.plt = 1; 
    p.type='sharp'; 
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
    if ~isfield(p,'means')
        p.means=0;
    end
    if ~isfield(p,'axh')
        p.axh=gca;
    end
    if ~isfield(p,'plt')
        p.plt=1;
    end
    if ~isfield(p,'type')
        p.type='sharp';
    end

    if ~isfield(p,'alpha')
        p.alpha=flip(logspace(log(0.3),log(0.6),numel(isovalue)));
    else
        p.alpha = p.alpha.*ones(1,numel(isovalue));
    end
end

% Setting sum of weights equal to 1
P=P./sum(P);

% Getting number of state vectors
N=size(X); 

switch N(2)
    case 2, [plots,shps]=plot_nongaussian_surface2D(X,P,isovalue,alpha,p);
    case 3, [plots,shps]=plot_nongaussian_surface3D(X,P,isovalue,alpha,p);
   otherwise
      error('Unsupported dimensionality');
end

function [plots,shps] = plot_nongaussian_surface2D(X,P,isovalue,alpha,p)
N=size(X); N=N(1); 

[P, sort_idx]=sort(P, "descend"); X=X(sort_idx,:); 

count_i=1; 
for i=isovalue
    sum_P=0; count_j=1;
    while (sum_P < i)
        Xs(count_j,:)=X(count_j,:);
        sum_P=sum_P+P(count_j);
        count_j=count_j+1;
        if count_j>N
            break
        end
    end 

    [k, shp] = boundary(Xs(:,1), Xs(:,2),alpha); 
    shps{count_i} = shp; 
    
    if p.plt
        if strcmp(p.type, 'sharp')
            if (count_i==1)&&p.display
                ploti = plot(shp, 'Parent', p.axh, 'Edgecolor', 'none', 'FaceColor', p.color{count_i}, 'FaceAlpha', p.alpha(count_i), 'DisplayName', p.name);
            else
                ploti = plot(shp, 'Parent', p.axh, 'Edgecolor', 'none', 'FaceColor', p.color{count_i}, 'FaceAlpha', p.alpha(count_i), 'HandleVisibility', 'off');
            end

            x_range = max(X(:,1)) - min(X(:,1));
            y_range = max(X(:,2)) - min(X(:,2));
            
            % Find maximum range among the axes
            max_range = max([x_range, y_range]);
            
            % Calculate scaling factors
            x_scale = (max_range / x_range)^-1;
            y_scale = (max_range / y_range)^-1;
            
            % Set the daspect ratio
            daspect([x_scale, y_scale, 1]);
        elseif strcmp(p.type, 'smooth')
            % Interpolating at discretized grid
            nbins = 200; 
            dx=[(max(X(:,1))-min(X(:,1)))/nbins; (max(X(:,2))-min(X(:,2)))/nbins];
            bin_centers=[min(X(:,1))+dx(1)/2:dx(1):max(X(:,1))-dx(1)/2; min(X(:,2))+dx(2)/2:dx(2):max(X(:,2))-dx(2)/2];
            [xq, yq]=meshgrid(bin_centers(1,:), bin_centers(2,:));
        
            xqf=reshape(xq,1,[]); yqf=reshape(yq,1,[]); Pqf=zeros(size(xqf));
        
            for j=1:length(xqf)
                Pqf(j)=inShape(shp,xqf(j),yqf(j));
            end
        
            Pq=reshape(Pqf,size(xq)); Pq= conv2(Pq,(1/16)*ones(4),'same');
            if (count_i==1)&&p.display
                [~,ploti] = contour(p.axh,xq,yq,Pq,1,'Edgecolor', 'none', 'FaceColor', p.color{count_i}, 'FaceAlpha', p.alpha(count_i), 'DisplayName', p.name, 'Fill', 'on');
            else
                [~,ploti] = contour(p.axh,xq,yq,Pq,1,'Edgecolor', 'none', 'FaceColor', p.color{count_i}, 'FaceAlpha', p.alpha(count_i), 'HandleVisibility', 'off', 'Fill', 'on');
            end
        else
            error('Unsupported plot type.')
        end
        plots{count_i} = ploti; 
    else
        plots{count_i} = NaN;
    end
    count_i=count_i+1; 
end

if p.plt, if p.means, mean_X=sum(X.*P); plots{1}=scatter(p.axh, mean_X(1), mean_X(2), 100, p.color{1}, "pentagram", "filled", 'HandleVisibility', 'off'); end, end

function [plots,shps] = plot_nongaussian_surface3D(X,P,isovalue,alpha,p)
N=size(X); N=N(1); 

if p.plt, if p.means, mean_X=sum(X.*P); plots{1}=scatter3(p.axh, mean_X(1), mean_X(2), mean_X(3), 100, p.color{1}, "pentagram", "filled", 'HandleVisibility', 'off'); end, end

[P, sort_idx]=sort(P, "descend"); X=X(sort_idx,:); 

count_i=1; 
for i=isovalue
    sum_P=0; count_j=1;
    while (sum_P < i)
        Xs(count_j,:)=X(count_j,:);
        sum_P=sum_P+P(count_j);
        count_j=count_j+1;
        if count_j>N
            break
        end
    end

    if isfield(p,'hgt')
        Xt(:,1) = [Xs(:,1) + (length(isovalue)-count_i+1)*p.hgt; Xs(:,1) + (length(isovalue)-count_i)*p.hgt];
        Xt(:,2) = [Xs(:,2); Xs(:,2)];
        Xt(:,3) = [Xs(:,3); Xs(:,3)];
        Xs = Xt; clear Xt; 
    end
    
    [k, shp] = boundary(Xs(:,1), Xs(:,2), Xs(:,3),alpha); 
    shps{count_i} = shp;

    if p.plt
        if strcmp(p.type, 'sharp')
            if (count_i==1)&&p.display
                ploti = trisurf(k, Xs(:,1), Xs(:,2), Xs(:,3), 'Parent', p.axh, 'Edgecolor', 'none', 'FaceColor', p.color{count_i}, 'FaceAlpha', p.alpha(count_i), 'DisplayName', p.name);  
            else
                ploti = trisurf(k, Xs(:,1), Xs(:,2), Xs(:,3), 'Parent', p.axh, 'Edgecolor', 'none', 'FaceColor', p.color{count_i}, 'FaceAlpha', p.alpha(count_i), 'HandleVisibility', 'off');
            end
        elseif strcmp(p.type, 'smooth')
            % Interpolating at discretized grid
            nbins = 200; 
            dx=[(max(X(:,1))-min(X(:,1)))/nbins; (max(X(:,2))-min(X(:,2)))/nbins; (max(X(:,3))-min(X(:,3)))/nbins];
            bin_centers=[min(X(:,1))+dx(1)/2:dx(1):max(X(:,1))-dx(1)/2; min(X(:,2))+dx(2)/2:dx(2):max(X(:,2))-dx(2)/2; min(X(:,3))+dx(3)/2:dx(3):max(X(:,3))-dx(3)/2];
            [xq, yq, zq]=meshgrid(bin_centers(1,:), bin_centers(2,:), bin_centers(3,:));
        
            xqf=reshape(xq,1,[]); yqf=reshape(yq,1,[]); zqf=reshape(zq,1,[]); Pqf=zeros(size(xqf));
        
            for j=1:length(xqf)
                Pqf(j)=inShape(shp,xqf(j),yqf(j),zqf(j));
            end
        
            Pq=reshape(Pqf,size(xq)); Pq= convn(Pq,(1/16)*ones(4),'same'); 
            if (count_i==1)&&p.display
                ploti = patch(p.axh,isosurface(xq,yq,zq,Pq,0.01),'Edgecolor', 'none', 'FaceColor', p.color{count_i}, 'FaceAlpha', p.alpha(count_i), 'DisplayName', p.name);  
            else
                ploti = patch(p.axh,isosurface(xq,yq,zq,Pq,0.01),'Edgecolor', 'none', 'FaceColor', p.color{count_i}, 'FaceAlpha', p.alpha(count_i), 'HandleVisibility', 'off'); 
            end
        else
            error('Unsupported plot type.')
        end
        plots{count_i}=ploti; 
    else
        plots{count_i}=NaN;
    end
    count_i=count_i+1;
end