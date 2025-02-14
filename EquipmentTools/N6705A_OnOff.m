function [ ] = N6705A_OnOff( N6705A, Channel, ChangeTo )
%N6705A_OnOff Allows N6705A supply channels to be turned off and on
%   N/A

if strcmp(ChangeTo, 'On') | strcmp(ChangeTo, 'on')
    message = sprintf('OUTP ON,(@%.0f)', Channel);
elseif strcmp(ChangeTo, 'Off') | strcmp(ChangeTo, 'off')
   message = sprintf('OUTP OFF,(@%.0f)', Channel);
else
    disp('Unknown ChangeTo for N6705A_GetVal')
end
fprintf(N6705A, message);
end

