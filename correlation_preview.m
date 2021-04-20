function fig = correlation_preview(data, varnames, varargin) 
% ---------------------------------------------------------
% Author: Dominika
% Fcn:  1) calculates correlation coefficients and p values across all
%       variables
%       2) identifies significantly correlated pairs of variables

% Input:    data = 2D matrix
%           --> rows = observations, columns = variables
%           varargin = 'alpha' + required value (default 0.05)
% ---------------------------------------------------------
% set alpha
if isempty(varargin)
    alpha = 0.05;
else
    a = find(strcmpi(varargin, 'alpha'));
    if ~isempty(a)
        alpha = varargin{a + 1};
    end
end

% calculate correlation coefficients and p values
[cor_coef, cor_p] = corrcoef(data);

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
            
            % choose significant correlation
            p_value = cor_p(b, c);
            if p_value < alpha
                h3 = text(0.5, 1.15, sprintf('p = %1.3f', p_value), 'Units', 'Normalized');
                set(h3, 'fontweight', 'bold'); 
                set(gca,'Color', [1 0.67 0.67]);
            end
        else
            nexttile 
        end       
        
        % add variable names
        if b == size(data, 2)
            xl = xlabel(varnames(c));
            xl.FontSize = 20; xl.FontWeight = 'bold';
        end        
        if mod(c, size(data, 2)) == 1
            yl = ylabel(varnames(b));
            yl.FontSize = 20; yl.FontWeight = 'bold';
        end
    end
end

end

function maximize(fig)
menubar_px = 120;
taskbar_px = 40;
scrsz = get(groot,'ScreenSize');
set(fig,'Position',[1 scrsz(2)+taskbar_px scrsz(3) scrsz(4)-menubar_px]);
end