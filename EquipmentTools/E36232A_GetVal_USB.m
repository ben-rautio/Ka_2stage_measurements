function [ Value ] = E36232A_GetVal_USB( E36232A, ValueType )
%E36232A_GetVal_USB Returns the current or voltage on the specified Channel
%   N/A
Channel = 1;
if ValueType == 'voltage' | ValueType == 'Voltage'
    message = sprintf('MEAS:VOLT? (@%.0f)', Channel);
elseif ValueType == 'current' | ValueType == 'Current'
    message = sprintf('MEAS:CURR? (@%.0f)', Channel);
else
    disp('Unknown Value Type for E36232A_GetVal')
end
Value = str2num(query(E36232A, message));
end

