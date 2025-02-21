function [Power] = NRP67T_ReadPower(NRP67T, freq)

% fwrite(NRP67T,"*RST");
pause(1)
fwrite(NRP67T,['SENS:FREQ ' num2str(freq)]);
fwrite(NRP67T,"INIT:ALL");
% fwrite(NRP67T,"TRIG:LEV -25 DBM");
% pause(5)
% fwrite(NRP67T,"TRIG:COUN 20");

message =  sprintf('FETCH%d?',1);
Power = str2double(query(NRP67T, message));

if isempty(Power)
    Power = 0;
end

end
