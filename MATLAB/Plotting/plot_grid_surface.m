function [plots, shps] = plot_grid_surface(X,P,isovalue,alpha,p)
% plot_grid_surface.m
% Benjamin Hanson, 2024
% 
% Given a set of 2D/3D grid centers X with associated weights P,
% generate an isosurface representing a curve of isovalue
% 
% Inputs:
%          X -- set of 2D/3D state grid centers
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
    p.color=[1 0 0];
    p.display=0; 
    p.means=0;
    p.axh=gca; 
    p.plt = 1; 
else
    if ~isfield(p,'color')
        p.color=[1 0 0];
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

    if ~isfield(p,'alpha')
        p.alpha=logspace(log(0.3),log(0.6),numel(isovalue));
    else
        p.alpha = p.alpha.*ones(1,numel(isovalue));
    end
end

% Setting sum of weights equal to 1
P=P./sum(P);

% Getting number of state vectors
N=size(X); 

switch N(2)
    case 2, [plots,shps]=plot_grid_surface2D(X,P,isovalue,alpha,p);
    case 3, [plots,shps]=plot_grid_surface3D(X,P,isovalue,alpha,p);
    otherwise
      error('Unsupported dimensionality');
end

function [plots,shps] = plot_grid_surface2D(X,P,isovalue,alpha,p)
x_list = unique(X(:,1));
y_list = unique(X(:,2));

Pfull=zeros(length(y_list),length(x_list));

mean_X = 0; mean_Y = 0; weight_sum = 0; 
for m=1:length(P)
    x_val=X(m,1); y_val=X(m,2); 
    i=find(x_list==x_val); j=find(y_list==y_val); 
    Pfull(j,i)=Pfull(j,i)+P(m);

    mean_X = mean_X + x_val*P(m); 
    mean_Y = mean_Y + y_val*P(m); 
    weight_sum = weight_sum + P(m); 
end
mean_X = mean_X/weight_sum; mean_Y = mean_Y/weight_sum; 

Pfull = Pfull./max(Pfull,[],'all'); 
c = contourc(x_list,y_list,Pfull,isovalue);

contourSets = struct(); i = 1; 
while i < size(c, 2)
    level = c(1, i); level = find(isovalue==level); 
    numPoints = c(2, i);
    points = c(:, i + 1 : i + numPoints);
    fieldName = ['level' strrep(num2str(level), '.', '_')];
    if isfield(contourSets, fieldName)
        contourSets.(fieldName) = [contourSets.(fieldName), points];
    else
        contourSets.(fieldName) = points;
    end
    i = i + numPoints + 1;
end

count_i=1; levelNames = fieldnames(contourSets);
for i=isovalue
    pts = contourSets.(levelNames{count_i}); pts = pts';

    [k, shp] = boundary(pts(:,1), pts(:,2), alpha); 
    shps{count_i} = shp; 

    if p.plt
        if (count_i==1)&&p.display
            ploti = plot(shp, 'Parent', p.axh, 'Edgecolor', 'none', 'FaceColor', p.color, 'FaceAlpha', p.alpha(count_i), 'DisplayName', p.name);
        else
            ploti = plot(shp, 'Parent', p.axh, 'Edgecolor', 'none', 'FaceColor', p.color, 'FaceAlpha', p.alpha(count_i), 'HandleVisibility', 'off');
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
      
        plots{count_i} = ploti; 
    else
        plots{count_i} = NaN;
    end
    count_i=count_i+1; 
end

if p.plt, if p.means, plots{1}=scatter(p.axh, mean_X, mean_Y, 100, p.color, "pentagram", "filled", 'HandleVisibility', 'off'); end, end

function [plots,shps] = plot_grid_surface3D(X,P,isovalue,alpha,p)
if isfield(p,'hgt')
    x_list = unique(X(:,2));
    y_list = unique(X(:,3));
    
    Pfull=zeros(length(y_list),length(x_list));
    
    mean_X = 0; mean_Y = 0; weight_sum = 0; 
    for m=1:length(P)
        x_val=X(m,2); y_val=X(m,3); 
        i=find(x_list==x_val); j=find(y_list==y_val); 
        Pfull(j,i)=Pfull(j,i)+P(m);
    
        mean_X = mean_X + x_val*P(m); 
        mean_Y = mean_Y + y_val*P(m); 
        weight_sum = weight_sum + P(m); 
    end
    mean_X = mean_X/weight_sum; mean_Y = mean_Y/weight_sum; 
    
    Pfull = Pfull./max(Pfull,[],'all'); 
    c = contourc(x_list,y_list,Pfull,isovalue);
    
    contourSets = struct(); i = 1; 
    while i < size(c, 2)
        level = c(1, i); level = find(isovalue==level); 
        numPoints = c(2, i);
        points = c(:, i + 1 : i + numPoints);
        fieldName = ['level' strrep(num2str(level), '.', '_')];
        if isfield(contourSets, fieldName)
            contourSets.(fieldName) = [contourSets.(fieldName), points];
        else
            contourSets.(fieldName) = points;
        end
        i = i + numPoints + 1;
    end
    
    count_i=1; levelNames = fieldnames(contourSets);
    for i=isovalue
        pts = contourSets.(levelNames{count_i}); pts = pts';
        Xt(:,1) = [X(1,1).*ones(size(pts(:,1))) + (count_i+1)*p.hgt; X(1,1).*ones(size(pts(:,1))) + (count_i)*p.hgt];
        Xt(:,2) = [pts(:,1) - mean_X + 0.03; pts(:,1) - mean_X + 0.03];
        Xt(:,3) = [pts(:,2) - mean_Y; pts(:,2) - mean_Y];
        pts = Xt; clear Xt; 
    
        [k, shp] = boundary(pts(:,1), pts(:,2), pts(:,3), alpha); 
        shps{count_i} = shp; 
    
       if p.plt
            if (count_i==1)&&p.display
                ploti = trisurf(k, pts(:,1), pts(:,2), pts(:,3), 'Parent', p.axh, 'Edgecolor', 'none', 'FaceColor', p.color{count_i}, 'FaceAlpha', p.alpha(count_i), 'DisplayName', p.name);  
            else
                ploti = trisurf(k, pts(:,1), pts(:,2), pts(:,3), 'Parent', p.axh, 'Edgecolor', 'none', 'FaceColor', p.color{count_i}, 'FaceAlpha', p.alpha(count_i), 'HandleVisibility', 'off');
            end
            
            %{
            x_range = max(X(:,1)) - min(X(:,1));
            y_range = max(X(:,2)) - min(X(:,2));
            
            % Find maximum range among the axes
            max_range = max([x_range, y_range]);
            
            % Calculate scaling factors
            x_scale = (max_range / x_range)^-1;
            y_scale = (max_range / y_range)^-1;
            
            % Set the daspect ratio
            daspect([x_scale, y_scale, 1]);
            %}
          
            plots{count_i} = ploti; 
        else
            plots{count_i} = NaN;
        end
        count_i=count_i+1; 
    end
else
    x_list = unique(X(:,1));
    y_list = unique(X(:,2));
    z_list = unique(X(:,3));
    
    Pfull=zeros(length(y_list),length(x_list), length(z_list));
    
    mean_X = 0; mean_Y = 0; mean_Z = 0; weight_sum = 0; 
    for m=1:length(P)
        x_val=X(m,1); y_val=X(m,2); z_val=X(m,3); 
        i=find(x_list==x_val); j=find(y_list==y_val); k=find(z_list==z_val); 
        Pfull(j,i,k)=Pfull(j,i,k)+P(m);
    
        mean_X = mean_X + x_val*P(m); 
        mean_Y = mean_Y + y_val*P(m); 
        mean_Z = mean_Z + z_val*P(m); 
        weight_sum = weight_sum + P(m); 
    end
    mean_X = mean_X/weight_sum; mean_Y = mean_Y/weight_sum; mean_Z = mean_Z/weight_sum; 
    
    if p.plt, if p.means, plots{1}=scatter3(p.axh, mean_X, mean_Y, mean_Z, 100, p.color, "pentagram", "filled", 'HandleVisibility', 'off'); end, end
    
    max_P = max(Pfull,[],'all'); Pfull = Pfull.*(1/max_P); 
    
    count_i=1; 
    for i=isovalue
        pa = patch(isosurface(x_list,y_list,z_list,Pfull,i));
    
        pts = pa.Vertices; 
    
        [k, shp] = boundary(pts(:,1), pts(:,2), pts(:,3), alpha); 
        shps{count_i} = shp; 
    
        if p.plt
            if (count_i==1)&&p.display
                ploti = trisurf(k, pts(:,1), pts(:,2), pts(:,3), 'Parent', p.axh, 'Edgecolor', 'none', 'FaceColor', p.color, 'FaceAlpha', p.alpha(count_i), 'DisplayName', p.name);  
            else
                ploti = trisurf(k, pts(:,1), pts(:,2), pts(:,3), 'Parent', p.axh, 'Edgecolor', 'none', 'FaceColor', p.color, 'FaceAlpha', p.alpha(count_i), 'HandleVisibility', 'off');
            end
                plots{count_i}=ploti; 
        else
            plots{count_i}=NaN;
        end
        count_i=count_i+1;
    end
end