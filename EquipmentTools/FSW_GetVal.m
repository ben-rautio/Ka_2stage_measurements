function [ Value ] = FSW_GetVal(FSW, ValueType)
%N6705A_ReadVal Returns the current or voltage on the specified Channel
%   N/A

if ValueType == 'channel' | ValueType == 'Channel'
    message = sprintf('CALC:MARK:FUNC:POW:RES? CPOW'); % Channel Power
elseif ValueType == 'Sidebar' 
    message = sprintf('CALC:MARK:FUNC:POW:RES? ACP'); % Adj Power
end


Value = str2num(query(FSW, message));
end

