function initialize_figures(varargin)
% initialize_figures.M
% Benjamin Hanson, 2024
% 
% Initialize figures for standard images
% 
% Inputs:
% n       -- Figures numbers 
%            * [n1]        -- Single figure
%            * [n1,...,nN] -- List of figures
% spacing -- Size of figures
%            * 'even'          -- Evenly spaced in screen
%            * 'equal'         -- Equal sized filling the screen
%            * Position vector -- {Pos1,...,PosN} specify [left, bottom,
% bg      -- Background color of figure
%            * 'white'
%            * 'black'
%            * 'clear'
%            width, height]
% lgd     -- Legend object
%            * Location
%            * Orientation
%            * Interpreter
%            * FontSize
%            * FontName
% ttl     -- Title object
%            * String
%            * Interpreter
%            * FontSize
%            * FontName
% lbl     -- Label object
%            * XString 
%            * YString
%            * ZString
%            * XPosition
%            * YPosition
%            * ZPosition
%            * Interpreter
%            * FontSize
%            * FontName
% axs     -- Axis type
% margin  -- Image margins
%            * [margin_x margin_y]
% vw      -- View object
% lmt     -- Limits object
%            * XLimit 
%            * YLimit 
%            * ZLimit
% lght     -- Lighting
%            * Position
%            * Type

% Default
n          = 1; 
screen_size = get(0, 'ScreenSize');
w = screen_size(3);
h = screen_size(4);
spacing    = {[100 100 w-200 h-200]};
bg_exist  = 0; 
lgd_exist  = 0; 
ttl_exist  = 0;
lbl_exist  = 0;
axs_exist  = 0;
lmts_exist = 0; lmt.XLimit = [-inf inf]; lmt.YLimit = [-inf inf]; lmt.ZLimit = [-inf inf];  
lght_exist = 0; 
margin     = [20 20]; 
vw_exist = 0; 

% Optionals
for i=1:2:nargin
    if strcmp('n',varargin{i})
        n = varargin{i+1};
    elseif strcmp('spacing',varargin{i})
        spacing = varargin{i+1};
    elseif strcmp('bg',varargin{i})
        bg = varargin{i+1};
        bg_exist = 1; 
    elseif strcmp('lgd',varargin{i})
        lgd = varargin{i+1};
        if lgd==1
            clear lgd; 
            lgd.Location = {'northeast'};
        end
        lgd_exist = 1;
    elseif strcmp('ttl',varargin{i})
        ttl = varargin{i+1};
        ttl_exist = 1; 
    elseif strcmp('lbl',varargin{i})
        lbl = varargin{i+1};
        lbl_exist = 1; 
    elseif strcmp('axs',varargin{i})
        axs = varargin{i+1};
        axs_exist = 1; 
    elseif strcmp('margin',varargin{i})
        margin = varargin{i+1};
    elseif strcmp('vw',varargin{i})
        vw = varargin{i+1};
        vw_exist = 1;
    elseif strcmp('lmt',varargin{i})
        lmt = varargin{i+1};
        lmts_exist = 1; 
    elseif strcmp('lght',varargin{i})
        lght = varargin{i+1};
        if lght==1
            clear lght; 
            lght.Position = [1 -1 1];
            lght.Type = 'phong'; 
        end
        lght_exist = 1;
    else
        error(append("Unspecified argument: ", varargin{i}));
    end
end

if lmts_exist
    if ~isfield(lmt,'XLimit')
        lmt.XLimit = [-inf inf]; 
    end
    if ~isfield(lmt,'YLimit')
        lmt.YLimit = [-inf inf]; 
    end
    if ~isfield(lmt,'ZLimit')
        lmt.ZLimit = [-inf inf]; 
    end
end

if lgd_exist
    if ~isfield(lgd,'Location')
        lgd.Location = {'northeast'};  
    end
    if ~isfield(lgd,'Orientation')
        lgd.Orientation = {'vertical'};  
    end
    if ~isfield(lgd,'Interpreter')
        lgd.Interpreter = 'latex';  
    end
    if ~isfield(lgd,'FontSize')
        lgd.FontSize = 12;  
    end
    if ~isfield(lgd,'FontName')
        lgd.FontName = 'Times';  
    end
end

if ~bg_exist
    bg = 'white';
end

if ttl_exist
    if ~isfield(ttl,'Interpreter')
        ttl.Interpreter = 'latex';  
    end
    if ~isfield(ttl,'FontSize')
        ttl.FontSize = 18;  
    end
    if ~isfield(ttl,'FontName')
        ttl.FontName = 'Times';  
    end
end

if lbl_exist
    if ~isfield(lbl,'Interpreter')
        lbl.Interpreter = 'latex';  
    end
    if ~isfield(lbl,'FontSize')
        lbl.FontSize = 16;  
    end
    if ~isfield(lbl,'FontName')
        lbl.FontName = 'Times';  
    end
end

if strcmp(spacing, 'even')
    spacing = {};
    
    c = ceil(sqrt(numel(n)));
    r = ceil(numel(n) / c);
    
    % Calculate figure size
    fig_width = floor(w / c);
    fig_height = floor(h / r);
    
    % Loop to create figures
    for i = 1:numel(n)
        % Calculate row and column index
        ri = ceil(i / c);
        ci = mod(i - 1, c) + 1;
        
        % Calculate figure position
        fig_left = (ci - 1) * fig_width + (margin(1)/2);
        fig_bottom = (r - ri) * fig_height + (margin(2)/2);

        spacing{i} = [fig_left fig_bottom fig_width-margin(1) fig_height*.8-margin(2)]; 
    end
elseif strcmp (spacing, 'equal')
    spacing = {};
    for i = 1:numel(n)
        spacing{i} = [margin(1)/2 margin(2)/2 w-margin(1) h-margin(2)]; 
    end
else
    if numel(spacing)==1
        spacing_temp = {};
        for i = 1:numel(n)
            spacing_temp{i} = spacing{1};
        end
        spacing = spacing_temp; 
    end
end

if lght_exist
    if ~isfield(lght,'Position')
        lght.Position = [1 -1 1];  
    end
    if ~isfield(lght,'Type')
        lght.Type = 'phong';  
    end
end

for i = 1:numel(n)
    fi = figure(n(i)); clf; hold all;
    fi.Position = spacing{i}; 

    if lgd_exist
        if numel(lgd.Location)==1
            legend('Location', lgd.Location,'Orientation', lgd.Orientation,'Interpreter', lgd.Interpreter,'FontSize', lgd.FontSize,'FontName',lgd.FontName); 
        else
            legend('Location', lgd.Location{i},'Orientation', lgd.Orientation{i},'Interpreter', lgd.Interpreter,'FontSize', lgd.FontSize,'FontName',lgd.FontName);
        end
    end

    if ttl_exist
        if numel(ttl.String)==1
            title(ttl.String,'Interpreter', ttl.Interpreter,'FontSize', ttl.FontSize,'FontName',ttl.FontName);
        else
            title(ttl.String{i},'Interpreter', ttl.Interpreter,'FontSize', ttl.FontSize,'FontName',ttl.FontName);
        end
    end

    if lbl_exist
        if isfield(lbl,'XString')
            xlabel(lbl.XString,...
                   'Interpreter', lbl.Interpreter,...
                   'FontSize', lbl.FontSize,...
                   'FontName',lbl.FontName); 
        end
        if isfield(lbl,'YString')
            ylabel(lbl.YString,...
                   'Interpreter', lbl.Interpreter,...
                   'FontSize', lbl.FontSize,...
                   'FontName',lbl.FontName); 
        end
        if isfield(lbl,'ZString')
            zlabel(lbl.ZString,...
                   'Interpreter', lbl.Interpreter,...
                   'FontSize', lbl.FontSize,...
                   'FontName',lbl.FontName);     
        end
    end

    if axs_exist
        if numel(axs)==1
            if strcmp(axs, 'equal'), axis equal; 
            elseif strcmp(axs, 'square'), axis square; 
            end
        else
            if strcmp(axs{i}, 'equal'), axis equal; 
            elseif strcmp(axs{i}, 'square'), axis square; 
            end
        end
    end

    if vw_exist
        if numel(vw)==1
            vw1 = vw{1};
            view(vw1(1), vw1(2)); 
        else
            vwi = vw{i}; 
            view(vwi(1), vwi(2)); 
        end
    end

    if lght_exist
        light('Position', lght.Position);
        if strcmp(lght.Type, 'phong')
            lighting phong;
        elseif strcmp(lght.Type, 'gouraud')
            lighting gouraud;
        elseif strcmp(lght.Type, 'flat')
            lighting flat;
        end
    end

    if strcmp(bg, 'white')
        set(gca, 'Color', 'white', 'XColor', 'black', 'YColor', 'black', 'ZColor', 'black');
    elseif strcmp(bg, 'black')
        set(gca, 'Color', 'black', 'XColor', 'white', 'YColor', 'white', 'ZColor', 'white');
    elseif strcmp(bg, 'clear')
        set(gca, 'Color', 'none', 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none');
    end

    xlim(lmt.XLimit); ylim(lmt.YLimit); zlim(lmt.ZLimit); 
end