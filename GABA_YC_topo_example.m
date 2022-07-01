%% GABA-AD: TEPs - TOPOGRAPHICAL ANALYSIS - EXAMPLE AMPLITUDE EXTRACTION
% Written by Dominika for GABA-AD project (2022)
% 
% Colection of scripts to visualize the outcome of the topographical analysis performed in Ragu 
%   --> figures are saved in a folder 'GABA_YC_figures'
% 
% 1) load GFP data
% 2) compare GFP
% 3) plot GFP
% 4) load peak data
% 5) plot peak data
% 6) plot TEP - example peak/subject

clear all; clc

%% parameters
% ----- adjustable parameters -----
% dataset
prefix = 'GABA';
group = 'YC';
participant = 1:20;

% visualization (in ms)
time_window = [-50, 300];
analysis_window = [10, 300];
x_delta = 0.5;
shade = 0.2;
col = [0.0824, 0.3373, 0.5098; 0.8000, 0.1961, 0.1529];

% statistics
z = 1.96;
alpha = 0.05;
% --------------------------------

% navigate to the Git folder --> source of saved default files
folder_git = uigetdir(pwd, 'Choose the Git folder');

% visualization 
figure_counter = 1;
x = time_window(1):x_delta:time_window(2);

% input/output folders
folder_results = uigetdir(pwd, 'Choose the Results folder');
cd(folder_results)
folder_figures = [folder_results '\GABA_' group '_figures'];

%% 1) LOAD GFP DATA
% mean GFP 
load([folder_results '\GABA_' group '_variables\GABA_' group '_M1_TEPs.mat'], 'GABA_GFP_mean')
GFP_M1 = GABA_GFP_mean;
% load([folder_results '\GABA_' group '_variables\GABA_' group '_AG_TEPs.mat'], 'GABA_GFP_mean')
% GFP_AG = GABA_GFP_mean;

% choose GFP from AG and sub-threshold M1
data_mean(1, :, :) = squeeze(GFP_M1(:, 1, 1, :));
data_mean(2, :, :) = squeeze(GFP_M1(:, 1, 2, :));
disp(['data_mean - datasize: ' num2str(size(data_mean))])
clear GABA_GFP_mean GFP_AG GFP_M1

% TEP per subject
load([folder_results '\GABA_' group '_variables\GABA_' group '_M1_TEPs.mat'], 'GABA_data')
TEP_M1 = GABA_data;
% load([folder_results '\GABA_' group '_variables\GABA_' group '_AG_TEPs.mat'], 'GABA_data')
% TEP_AG = GABA_data;

% calculate individual GFP from AG and sub-threshold M1
data_subj(1, :, :, :) = squeeze(std(TEP_M1(:, 1, 1, :, 1:30, :), 1, 5));
data_subj(2, :, :, :) = squeeze(std(TEP_M1(:, 1, 2, :, 1:30, :), 1, 5));
% data_subj(2, :, :, :) = squeeze(std(TEP_AG(:, 1, :, 1:30, :), 1, 4));
disp(['data_subj - datasize: ' num2str(size(data_subj))])
clear TEP_M1 TEP_AG GABA_data

%% 2) COMPARE GFP 
% ----- adjustable parameters -----
dataset = [1, 2];
% ---------------------------------
n_perms = 1500;
threshold = norminv(1 - alpha/2);
n_subj = size(data_subj, 3);

% ----- identify time intervals of significant change in GFP -----
% create data matrix
data_perm = cat(2, squeeze(data_subj(1, dataset(1), :, :))', squeeze(data_subj(2, dataset(2), :, :))');

% generate true condition labels
labels = (1:n_subj*2) > n_subj;

% run the permutation
diffs_permuted = zeros(length(x), n_perms);
for p = 1:n_perms
    % shuffle condition label vector
    labels_mixed = labels(randperm(2*n_subj));

    % compute and store difference time series
    mean_1 = mean(data_perm(:, labels_mixed == 1), 2);
    mean_2 = mean(data_perm(:, labels_mixed == 0), 2);
    diffs_permuted(:, p) = mean_1 - mean_2;
end
clear labels_mixed mean_1 mean_2

% compute z-score difference
diff_observed = mean(data_perm(:, labels == 1), 2) - mean(data_perm(:, labels == 0), 2);
z_diff = (diff_observed - mean(diffs_permuted, 2)) ./ std(diffs_permuted, [], 2);

% statistically threshold the final result
z_thresholded = z_diff;
z_thresholded(abs(z_thresholded) < threshold) = 0;

% ----- correct for too small clusters -----
% find cluster sizes under the null hypothesis
cluster_size = zeros(n_perms, 1);
for p = 1:n_perms
    % compute z-score difference
    z_diff_perm = (diffs_permuted(:, p) - mean(diffs_permuted, 2)) ./ std(diffs_permuted,[],2);

    % threshold
    z_diff_perm(abs(z_diff_perm) < threshold) = 0;

    % identify clusters
    islands = bwconncomp(logical(z_diff_perm));

    % find cluster sizes
    if length(islands.PixelIdxList) == 0
        cluster_size(p) = 0;
    else
        clustNs = cellfun(@length, islands.PixelIdxList);
        cluster_size(p) = max(clustNs);
    end
end
clear z_diff_perm islands clustNs

% compute cluster threshold
cluster_threshold = prctile(cluster_size, 100 - alpha*100);

% show distribution of cluster sizes
fig = plot_hist(cluster_size, cluster_threshold, figure_counter);
figure_counter = figure_counter + 1;

% find islands
islands = bwconncomp(logical(z_thresholded));

% find and remove any subthreshold islands
z_thresholded_cc = z_thresholded;
if islands.NumObjects > 0
    for i = 1:islands.NumObjects
        if numel(islands.PixelIdxList{i}) < cluster_threshold
            z_thresholded_cc(islands.PixelIdxList{i}) = 0;
        end
    end
end
clear p i n_subj data_perm labels diffs_permuted diff_observed z_diff threshold n_perms threshold cluster_size islands fig

%% 3) PLOT GFP
% choose the data
for a = 1:2
    data_visual(a, :) = squeeze(data_mean(a, dataset(a), :, :));
end

% plot GFP timecourse
fig = plot_GFP(x, data_visual, figure_counter, analysis_window, z_thresholded, col);
figure_counter = figure_counter + 1 ;    

% name and save figure
figure_name = 'example_GFP';
savefig([folder_figures '\' figure_name '.fig'])
saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')   
clear a data_visual fig figure_name 

%% 4) LOAD PEAK DATA
% ----- adjustable parameters -----
dataset = [1, 2];
peak = [2, 2];
sub2include = logical(ones(1, 20));
% sub2include(2) = false;
% comparison = {'stimulus', 'sub-threshold', 'supra-threshold'};
comparison = {'session', '1', '2'};
% comparison = {'component', 'N17', 'P30'};
% sub2include = (data.intensity.IndFeature == 1)';
% ---------------------------------
clear data_amplitude data_latency
load([folder_results '\GABA_' group '_statistics\GABA_' group '_results.mat'])

% peak amplitude 
data_amplitude(1, :) = GABA_YC_results.TEP_M1(dataset(1)).amplitude.pre(sub2include, 2, peak(1))';
data_amplitude(2, :) = GABA_YC_results.TEP_M1(dataset(2)).amplitude.pre(sub2include, 2, peak(2))';
% data_amplitude(2, :) = GABA_YC_results.TEP_AG(dataset(2)).amplitude.pre(:, peak)';

% peak latency 
data_latency(1, :) = GABA_YC_results.TEP_M1(dataset(1)).latency.pre(sub2include, 2, peak(1))';
data_latency(2, :) = GABA_YC_results.TEP_M1(dataset(2)).latency.pre(sub2include, 2, peak(2))';
% data_latency(2, :) = GABA_YC_results.TEP_AG(dataset(2)).latency.pre(:, peak)';

% extract mean values
data_amplitude = data_amplitude(:, 1:end-1);
average = struct;
for a = 1:2
    average.amplitude(a, 1) = mean(data_amplitude(a, :));
    average.amplitude(a, 2) = std(data_amplitude(a, :));
    average.latency(a, 1) = mean(data_latency(a, :));
    average.latency(a, 2) = std(data_latency(a, :));
end
clear a

%% 5) PLOT PEAK DATA
% plot amplitude
data_visual = data_amplitude';
fig = plot_scatter(data_visual, figure_counter, col, 'datatype', 'amplitude', 'comparison', comparison);   
hold off

% name and save figure
figure_name = 'example_amplitude';
savefig([folder_figures '\' figure_name '.fig'])
saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')   

% update figure counter
figure_counter = figure_counter + 1;

% plot latency
data_visual = data_latency';
fig = plot_scatter(data_visual, figure_counter, col, 'datatype', 'latency', 'comparison', comparison);   
hold on

% name and save figure
figure_name = 'example_latency';
savefig([folder_figures '\' figure_name '.fig'])
saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')   

% update figure counter
figure_counter = figure_counter + 1;

clear data_visual fig figure_name sub2include comparison peak dataset

%% 6) PLOT TEP - EXAMPLE SUBJECT
% load TEP data
% ----- adjustable parameters -----
dataset = [1, 2];
peak = 2;
sub2include = 19;
window_visual = [-10, 60];
% ---------------------------------
% load data
load([folder_results '\GABA_' group '_variables\GABA_' group '_M1_TEPs.mat'], 'GABA_data', 'GABA_TEP_default')

% identify EOIs
load([folder_git '\GABA_header_default.mat'])
for e = 1:length(GABA_TEP_default.eoi.TS{peak})
    eois(e) = find(strcmp(GABA_TEP_default.eoi.TS{peak}{e}, {header.chanlocs.labels}));  
end

% select dataset
for d = 1:length(dataset)
    data_TEP(d, :) = mean(GABA_data(dataset(d), 1, 2, sub2include, eois, :), 5);
end

%plot TEPs
fig = plot_TEP(x, data_TEP, figure_counter, window_visual);

% name and save figure
figure_name = 'example_TEP';
savefig([folder_figures '\' figure_name '.fig'])
saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')  

clear e d header dataset peak sub2include fig

%% functions
function fig = plot_hist(cluster_size, cluster_threshold, figure_counter)
    fig = figure(figure_counter);
    hold on

    % plot histogram
    h = histogram(cluster_size);
    h.EdgeColor = 'white'; 
    
    % add threshold line
    xl = get(gca, 'xlim');
    l = line([xl(1), xl(2)], [cluster_threshold, cluster_threshold], 'color', [0.8000 0.1961 0.1529], 'linewidth', 2);

    % other parameters
    xlabel('Cluster size (time points)');
    ylabel('Count')
    set(gca, 'FontSize', 14) 
    set(gca, 'layer', 'top');
    hold off
end
function fig = plot_GFP(x, data_visual, figure_counter, analysis_window, z_thresholded, col)
    % launch the figure
    fig = figure(figure_counter);
    hold on

    % set limits of y
    yl = [0 3];
    ylim(yl)
           
    % plot intervals of significant difference
    for i = 1:2
        I(i) = area(x, logical(z_thresholded) * yl(i));
        I(i).FaceColor = [0.9882    0.7608    0.7608];
        I(i).EdgeColor = 'none';
        I(i).Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    
%     % shade area of GFP difference
%     diff = data_visual(2, :) - data_visual(1, :);
%     F = fill([x fliplr(x)],[diff zeros(1, length(x))], [0.5 0.5 0.5], 'linestyle', 'none');
%     F.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    % plot the GFP course
    for t = 1:size(data_visual, 1)
        plot(x, data_visual(t, :), 'Color', col(t, :), 'LineWidth', 4);
    end

    % other parameters
    xlabel('time (ms)')
    ylabel('amplitude(\muV)')
    set(gca, 'FontSize', 18)
    xlim([analysis_window(1) + 2, analysis_window(2)])   
    set(gca, 'layer', 'top');

    % change figure size
    fig.Position = [500 500 550 400];
    
    hold off
end
function fig = plot_scatter(data_visual, figure_counter, col, varargin)
    % determine data type
    a = find(strcmpi(varargin, 'datatype'));
    if ~isempty(a)
        datatype = varargin{a + 1};
    end

    % determine compared variables
    b = find(strcmpi(varargin, 'comparison'));
    if ~isempty(b)
        comparison = varargin{b + 1};
    end
    
    % launch the figure
    fig = figure(figure_counter);
    hold on
    
    % plot the markers
    data_group = cat(1, ones(size(data_visual, 1), 1), 2 * ones(size(data_visual, 1), 1));
    data_y = cat(1, data_visual(:, 1), data_visual(:, 2));
    beeswarm(data_group, data_y)
    
    % add boxplot
    for a = 1:2
        boxchart(a * ones(size(data_visual, 1), 1), data_visual(:, a), 'BoxFaceColor', col(a, :))
    end
    
    % reorder layers
    child = get(gca, 'Children');
    set(gca,'Children',[child(3) child(4) child(1) child(2)])
    
    % font
    set(gca, 'Fontsize', 18)
    
    % x label
    xlabel(comparison{1})
    set(gca, 'xtick', 1:size(data_visual, 2), 'xticklabel', comparison(2:end))
    
    % y label
    if strcmp(datatype, 'amplitude')
        ylabel('amplitude (\muV)')
    elseif strcmp(datatype, 'latency')
        ylabel('latency (ms)')
    end
end
function fig = plot_TEP(x, data_visual, figure_counter, window_visual)
    % launch the figure
    fig = figure(figure_counter);
    hold on
    
    % plot TEP 
    for p = 1:size(data_visual, 1)
        P(p) = plot(x, data_visual(p, :), 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    end
    
    % shade interpolated interval 
    yl = get(gca, 'ylim');
    rectangle('Position', [-5, yl(1), 15, yl(2) - yl(1)], 'FaceColor', [0.99 0.73 0.73], 'EdgeColor', 'none')
    
    % plot the zero line
    xl = get(gca, 'xlim');
    line(xl, [0, 0], 'color', [0, 0, 0], 'LineStyle', ':', 'LineWidth', 1.5)
    
    % TMS stimulus
    line([0, 0], yl, 'Color', [0.88 0.08 0.08], 'LineWidth', 3)
    
    % reorder layers
    child = get(gca, 'Children');
    set(gca,'Children',[child(1) child(2) child(4) child(5) child(3)])
    
    % crop
    xlim(window_visual)
    ylim([-5, yl(2)])
    
    % other parameters
    xlabel('time (ms)')
    ylabel('amplitude (\muV)')
    set(gca, 'FontSize', 18) 
    set(gca, 'layer', 'top');
    hold off
    
%     % change figure size
%     fig.Position = [500 500 500 300];
end
function beeswarm(x,y,varargin)
% function xbee = beeswarm(x,y)
%
% Input arguments:
%   x               column vector of groups (only tested for integer)
%   y               column vector of data
%
% Optional input arguments:
%   sort_style      ('nosort' - default | 'up' | 'down' | 'fan' | 'rand' | 'square' | 'hex')
%   corral_style    ('none' default | 'gutter' | 'omit' | 'rand')
%   dot_size        relative. default=1
%   overlay_style   (false default | 'box' | 'sd' | 'ci')
%   use_current_axes (false default | true)
%   colormap        (lines default | 'jet' | 'parula' | 'r' | Nx3 matrix of RGB values]
%
% Output arguments:
%   xbee            optimized layout positions
%
% Known Issues:
%       x locations depend on figure aspect ratio. resizing the figure window and rerunning may give different results
%       setting corral to 'none' still has a gutter when the width is large
%
% Usage example:
% 	x = round(rand(150,1)*5);
%   y = randn(150,1);
%   beeswarm(x,y,3,'sort_style','up','overlay_style','ci')
%
% % Ian Stevenson, CC-BY 2019

p = inputParser;
addRequired(p,'x')
addRequired(p,'y')
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addOptional(p,'sort_style','nosort')
addOptional(p,'corral_style','none')
addOptional(p,'dot_size',20/sqrt(length(x)),validScalarPosNum)
addOptional(p,'overlay_style',false)
addOptional(p,'use_current_axes',false)
addOptional(p,'colormap','lines')
addOptional(p,'MarkerFaceColor','')
addOptional(p,'MarkerFaceAlpha',1)
addOptional(p,'MarkerEdgeColor', 'black')
parse(p,x,y,varargin{:});

% extra parameters
rwid = .05; % width of overlay box/dash

dcut=8; % spacing factor
nxloc=512; % resolution for optimization
chanwid = .9; % percent width of channel to use
yl = [min(y) max(y)]; % default y-limits
asp_rat = 1;
keep_hold = false;

% get aspect ratio for a figure window
if isfinite(p.Results.dot_size)
    if ~p.Results.use_current_axes
        % make new axes
        s=scatter(x,y);
        xl=[min(x)-.5 max(x)+.5];
    else
        xl=xlim();
    end
    yl=ylim();
    pasp_rat = get(gca,'PlotBoxAspectRatio');
    dasp_rat = get(gca,'DataAspectRatio');
    asp_rat = pasp_rat(1)/pasp_rat(2);
    
    % pix-scale
    pf = get(gcf,'Position');
    pa = get(gca,'Position');
    as = pf(3:4).*pa(3:4); % width and height of panel in pixels
    dcut = dcut*sqrt(p.Results.dot_size)/as(1)*(range(unique(x))+1);
    cla
end

% sort/round y for different plot styles
yorig=y;
switch lower(p.Results.sort_style)
    case 'up'
        [y,sid]=sort(y);
    case 'fan'
        [~,sid]=sort(abs(y-mean(y)));
        sid=[sid(1:2:end); sid(2:2:end)];
        y=y(sid);
    case 'down'
        [y,sid]=sort(y,'descend');
    case 'rand'
        sid=randperm(length(y));
        y=y(sid);
    case 'square'
        nxloc=.9/dcut;
%         [~,e,b]=histcounts(y,ceil((range(x)+1)*chanwid*nxloc/2/asp_rat));
        edges = linspace(min(yl),max(yl),ceil((range(x)+1)*chanwid*nxloc/asp_rat));
        [~,e,b]=histcounts(y,edges);
        y=e(b)'+mean(diff(e))/2;
        [y,sid]=sort(y);
    case 'hex'
        nxloc=.9/dcut;
%         [~,e,b]=histcounts(y,ceil((range(x)+1)*chanwid*nxloc/2/sqrt(1-.5.^2)/asp_rat));
        edges = linspace(min(yl),max(yl),ceil((range(x)+1)*chanwid*nxloc/sqrt(1-.5.^2)/asp_rat));
        [n,e,b]=histcounts(y,edges);
        oddmaj=0;
        if sum(mod(n(1:2:end),2)==1)>sum(mod(n(2:2:end),2)==1),
            oddmaj=1;
        end
        y=e(b)'+mean(diff(e))/2;
        [y,sid]=sort(y);
        b=b(sid);
    otherwise
        sid=1:length(y);
end
x=x(sid);
yorig=yorig(sid);
[ux,~,ic] = unique(x);
% rmult=(range(ux)+1)*2;
rmult=5;

% for each group...
for i=1:length(ux)
    fid = find(ic==i);   
    
    % set of possible x locations
    xi = linspace(-chanwid/2*rmult,chanwid/2*rmult,nxloc*rmult+(mod(nxloc*rmult,2)==0))'+ux(i);

    % rescale y to that things are square visually
    zy=(y(fid)-min(yl))/(max(yl)-min(yl))/asp_rat*(range(ux)+1)*chanwid;
    
    % precalculate y distances so that we only worry about nearby points
    D0=squareform(pdist(zy))<dcut*2;    
    
    if length(fid)>1
        % for each data point in the group sequentially...
        for j=1:length(fid)
            if strcmp(lower(p.Results.sort_style),'hex')
                xi = linspace(-chanwid/2*rmult,chanwid/2*rmult,nxloc*rmult+(mod(nxloc*rmult,2)==0))'+ux(i);
                if mod(b(fid(j)),2)==oddmaj
                    xi = linspace(-chanwid/2*rmult,chanwid/2*rmult,nxloc*rmult+(mod(nxloc*rmult,2)==0))'+ux(i)+mean(diff(xi))/2;
                end
            end
            zid = D0(j,1:j-1);
            e = (xi-ux(i)).^2; % cost function
            if ~strcmp(lower(p.Results.sort_style),'hex') && ~strcmp(lower(p.Results.sort_style),'square')
                if sum(zid)>0
                    D = pdist2([xi ones(length(xi),1)*zy(j)], [x(fid(zid)) zy(zid)]);
                    D(D<=dcut)=Inf;
                    D(D>dcut & isfinite(D))=0;
                    e = e + sum(D,2) + randn(1)*10e-6; % noise to tie-break
                end
            else
                if sum(zid)>0
                    D = pdist2([xi ones(length(xi),1)*zy(j)], [x(fid(zid)) zy(zid)]);
                    D(D==0)=Inf;
                    D(D>dcut & isfinite(D))=0;
                    e = e + sum(D,2) + randn(1)*10e-6; % noise to tie-break
                end
            end

            if strcmp(lower(p.Results.sort_style),'one')
                e(xi<ux(i))=Inf;
            end
            [~,mini] = min(e);
            if mini==1 && rand(1)>.5, mini=length(xi); end
            x(fid(j)) = xi(mini);
        end
    end
%     x(fid)=x(fid)-median(x(fid))+ux(i); % center x locations by median
end

if strcmp(lower(p.Results.sort_style),'randn')
    x=ux(ic)+randn(size(ic))/4;
end

% corral any points outside of the channel
out_of_range = abs(x-ux(ic))>chanwid/2;
switch lower(p.Results.corral_style)
    case 'gutter'
        id = (x-ux(ic))>chanwid/2;
        x(id)=chanwid/2+ux(ic(id));
        id = (x-ux(ic))<-chanwid/2;
        x(id)=-chanwid/2+ux(ic(id));
    case 'omit'
        x(out_of_range)=NaN;
    case 'random'
        x(out_of_range)=ux(ic(out_of_range))+rand(sum(out_of_range),1)*chanwid-chanwid/2;
end

% plot groups and add overlay
if isfinite(p.Results.dot_size)
    if isnumeric(p.Results.colormap)
        cmap=p.Results.colormap;
    else
        cmap = feval(p.Results.colormap,length(ux));
    end
    for i=1:length(ux)
        if isempty(p.Results.MarkerFaceColor')
            scatter(x(ic==i),y(ic==i),p.Results.dot_size*36,'filled','MarkerFaceAlpha',p.Results.MarkerFaceAlpha,'MarkerEdgeColor',p.Results.MarkerEdgeColor,'MarkerFaceColor',cmap(i,:))
        else
            scatter(x(ic==i),y(ic==i),p.Results.dot_size*36,'filled','MarkerFaceAlpha',p.Results.MarkerFaceAlpha,'MarkerEdgeColor',p.Results.MarkerEdgeColor,'MarkerFaceColor',p.Results.MarkerFaceColor)
        end
        hold on
        iqr = prctile(yorig(ic==i),[25 75]);
        switch lower(p.Results.overlay_style)
            case 'box'
                rectangle('Position',[ux(i)-rwid iqr(1) 2*rwid iqr(2)-iqr(1)],'EdgeColor','k','LineWidth',2)
                line([ux(i)-rwid ux(i)+rwid],[1 1]*median(yorig(ic==i)),'LineWidth',3,'Color',cmap(i,:))
            case 'sd'
                line([1 1]*ux(i),mean(yorig(ic==i))+[-1 1]*std(yorig(ic==i)),'Color',cmap(i,:),'LineWidth',2)
                line([ux(i)-2*rwid ux(i)+2*rwid],[1 1]*mean(yorig(ic==i)),'LineWidth',3,'Color',cmap(i,:))
            case 'ci'
                line([1 1]*ux(i),mean(yorig(ic==i))+[-1 1]*std(yorig(ic==i))/sqrt(sum(ic==i))*tinv(0.975,sum(ic==i)-1),'Color',cmap(i,:),'LineWidth',2)
                line([ux(i)-2*rwid ux(i)+2*rwid],[1 1]*mean(yorig(ic==i)),'LineWidth',3,'Color',cmap(i,:))
        end
        
    end
    hold on
    xlim(xl)
    ylim(yl)
end

% unsort so that output matches the original y data
x(sid)=x;
end
