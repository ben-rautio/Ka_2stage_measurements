function [ ] = N6705A_SetVal( N6705A, Channel, ValueType, Value )
%N6705A_SetVal Sets the current or voltage to Value on the specified
%Channel
%   N/A
if ValueType == 'voltage' | ValueType == 'Voltage'
    message = sprintf('VOLT %.2f, (@%.0f)', Value, Channel);
elseif ValueType == 'current' | ValueType == 'Current'
    message = sprintf('CURR %.2f, (@%.0f)', Value, Channel);
else
    disp('Unkown Value Type for N6705A_SetVal')
end
fprintf(N6705A, message)
end

