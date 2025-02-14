function [ Value ] = N6705A_GetVal( N6705A, Channel, ValueType )
%N6705A_ReadVal Returns the current or voltage on the specified Channel
%   N/A

if ValueType == 'voltage' | ValueType == 'Voltage'
    message = sprintf('MEAS:VOLT? (@%.0f)', Channel);
elseif ValueType == 'current' | ValueType == 'Current'
    message = sprintf('MEAS:CURR? (@%.0f)', Channel);
else
    disp('Unknown Value Type for N6705A_GetVal')
end
Value = str2num(query(N6705A, message));
end

