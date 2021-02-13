function TEP_topoplot(header,data,x_pos, map_lims)
% ------------------------------------------------------------------------
% Fnc: Draws the topographical signal distribution in the current figure
% 
% Input: 
%   header - letswave header with dataset metadata
%   data - letswave 6D data matrix
%   x_pos - index of target timepoint
%   map_lims - colormap limits [-y, +y]
% 
% Author: Dominika (2020) 
% ------------------------------------------------------------------------

varargin = {'maplimits' map_lims 'shading' 'interp' 'whitebk' 'on'};

%fetch data to display
x_visual = ceil((x_pos - header.xstart)/header.xstep);
% x_visual = ceil((peak_widths.latency(peak) - header_visual.xstart)/header_visual.xstep);
vector = squeeze(data(1,:,1,1,1,x_visual));

%fetch chanlocs
chanlocs = header.chanlocs;

%parse data and chanlocs 
k=1;
for chanpos=1:size(chanlocs,2);
    vector2(k)=double(vector(chanpos));
    chanlocs2(k)=chanlocs(chanpos);
    k=k+1;
end;

topoplot(vector2,chanlocs2,varargin{:});
set(gcf,'color',[1 1 1]);