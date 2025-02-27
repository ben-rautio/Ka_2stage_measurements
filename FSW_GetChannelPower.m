function power = FSW_GetChannelPower(FSW)
message = sprintf('CALC:MARK:FUNC:POW:RES? CPOW'); %Get Tx channel pwoer
power = str2double(query(FSW,message));
end