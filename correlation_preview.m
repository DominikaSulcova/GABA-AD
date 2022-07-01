function [fig, cor_coef, cor_p] = correlation_preview(data, varnames, varargin) 
% ---------------------------------------------------------
% Author: Dominika
% Fcn:  1) calculates correlation coefficients and p values across all
%       variables
%       2) identifies significantly correlated pairs of variables

% Input:    data = 2D matrix
%           --> rows = observations, columns = variables
%           varargin:
%               - 'alpha' --> required value (default 0.05)
%               - 'method' --> 'pearson' (default) or ''
% ---------------------------------------------------------
% set varargins
if isempty(varargin)
    alpha = 0.05;
    method = 'Pearson';
else
    a = find(strcmpi(varargin, 'alpha'));
    b = find(strcmpi(varargin, 'method'));
    % alpha
    if ~isempty(a)
        alpha = varargin{a + 1};
    else
        alpha = 0.05;
    end
    % method
    if ~isempty(b)
        method = varargin{b + 1};
    else
        alpha = 0.05;
    end
end

% calculate correlation coefficients and p values
[cor_coef, cor_p] = corr(data, 'Type', method);

% in case of Spearman correlation, rank the data
if strcmp(method, 'Spearman')
    for c = 1:size(data, 2)
        [temp, data(:, c)]  = ismember(data(:, c), unique(data(:, c)));
    end
end

% launch the common figure
fig = figure;
hold on
maximize(fig)

% plot relationship of each variable pair
t = tiledlayout(size(data, 2), size(data, 2));
for b = 1:size(data, 2)
    for c = 1:size(data, 2)
        % plot the data + least squares line
        if b ~= c
            nexttile
            scatter(data(:, b),data(:, c), 'MarkerFaceColor', [0.02 0.4 0.7], 'MarkerEdgeColor', 'none');
            h1 = lsline; h1.Color = [0 0 0]; 
            h2 = text(0.05, 1.15, sprintf('r = %1.3f', cor_coef(b, c)), 'Units', 'Normalized');  
            set(h2, 'fontsize', 12);
            
            % choose significant correlation
            p_value = cor_p(b, c);
            if p_value < alpha
                h3 = text(0.5, 1.15, sprintf('p = %1.3f', p_value), 'Units', 'Normalized');
                set(h3, 'fontweight', 'bold', 'fontsize', 12); 
                set(gca,'Color', [1 0.67 0.67]);
            end
        else
            nexttile 
        end       
        
        % add variable names
        if b == size(data, 2)
            xl = xlabel(varnames(c));
            xl.FontSize = 16; xl.FontWeight = 'bold';
        end        
        if mod(c, size(data, 2)) == 1
            yl = ylabel(varnames(b));
            yl.FontSize = 16; yl.FontWeight = 'bold';
            yl.Rotation = 0; yl.VerticalAlignment = 'middle'; yl.HorizontalAlignment = 'right';
        end
    end
end
t.TileSpacing = 'compact';
t.Padding = 'compact';

end

function maximize(fig)
menubar_px = 120;
taskbar_px = 40;
scrsz = get(groot,'ScreenSize');
set(fig,'Position',[1 scrsz(2)+taskbar_px scrsz(3) scrsz(4)-menubar_px]);
end