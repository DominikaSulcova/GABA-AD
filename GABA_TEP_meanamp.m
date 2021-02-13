function [amplitude, averaged_x, averaged_data] = mean_TEP(data, polarity, center, span, percent, step, xstart)
% ------------------------------------------------------------------------
% Fnc: Calculates mean amplitude of the most prominent <percent> of the TOI 
%       --> TOI defined by the center value and the span of the time window
% 
% Input: 
%   data - vector of values (double), signal from target electrode
%   polarity -  string, values: 'positive', 'negative' - according to TEP peak  
%   center, span - num values that define the TOI
%   percent - how many top percent of datapoints will be included in the average
%   step, xstart - defines properties of time axes 
% 
% Author: Dominika (2020) 
% ------------------------------------------------------------------------

% prepare the interval vector 
interval = false(1, length(data));

% define the boundaries of the TOI - both limits are included in the interval!
start = round(((center-span/2) - xstart)/step);
stop = round(((center+span/2) - xstart)/step);
data_crop = data(start : stop);

% calculate number of points to average
points_number = ceil((percent/100) * length(data_crop));

% sort the data according to peak polarity
switch polarity
    case 'positive'
        data_sorted = sort(data_crop, 'descend');        
    case 'negative'
        data_sorted = sort(data_crop); 
end

% calculate the mean value
points_included = data_sorted(1:points_number);
amplitude = mean(points_included);

% index averaged points
for i = 1:points_number
    point_pos = find(data_crop == points_included(i)); 
    % in case there are two same values, take the earlier one
    index(i) = point_pos(1);
end
index = sort(index);

% fill ones in the interval vector at indexed positions
for i = 1:points_number
    interval(start + index(i) -1) = true;
end

% calculate the outcome vectors for visualisation
averaged_data = data .* interval;                                   % keep only values ov averaged datapoints, the rest is set to 0
averaged_index = find(averaged_data);                               % index the averaged datapoints
averaged_x = averaged_index * step + xstart;                        % set the time interval
averaged_data = averaged_data(find(averaged_data));                 % get rid of the zeros

end