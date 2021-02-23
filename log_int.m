function [XXi, YYi] = log_int(XX, YY)
% ---------------------------------------------------------
% Author: Michele, Dominika
% Fcn:  1) log-transforms data 
%       2) upsamples 4 times
% ---------------------------------------------------------

% avoid zeros 
ix= YY==0; YY(ix)=[]; XX(ix)=[]; 
ix= XX==0; YY(ix)=[]; XX(ix)=[]; 

% log-transform
 X= log10(XX);
 Y= log10(YY);

% interpolate in log-log: X, Y --> Xi, Yi
 stepsFrex= (diff( XX)); stepsFrex= round(stepsFrex.*1000)./1000;
 if length( unique(stepsFrex))==1 %IF frex is a vector of linearly increasing frex, try
     XXi= logspace( X(1), X(end),  length(X)*4); % 2^10*2;--> auto-chosen
      Xi= log10(XXi); 
     Yi= interpn( X, Y,  (Xi) );% interpn ,'spline', 'nearest' 
     YYi= 10.^(Yi); 
 else % catch
     XXi= XX;  Xi= log10(XXi); 
     YYi= YY;  Yi= log10(YYi); disp('could not upsample PSD')
 end
end