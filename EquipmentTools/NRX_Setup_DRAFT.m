% Written by Grant G. - Updated by Adam Der
% Connects to NRX - Power Meter
% Reads from NRP67T

NRP67T = visa('ni', "TCPIP::192.168.0.10::INSTR");
fopen(NRP67T);

dev_name = query(NRP67T, "*IDN?")

if ~contains(dev_name, "Rohde&Schwarz")
   disp("Ah crap it connected to the wrong thing. Ooops"); 
   return;
else
    disp("Ay it connected!");
end

%{
ChannelOut = 2;
message =  sprintf('FETC%d?',ChannelOut);
Power = str2num(query(dev, message));
%}

fwrite(NRP67T,"*RST");
%testing = query(dev,"TEST:SENSor?")
%
%fwrite(NRP67T,"SENS:POW:BURS"); % This is not working for new power sensor
%fwrite(NRP67T,"SENS:POW:BURS:DTOL 1e-6");
fwrite(NRP67T,"SENS:FREQ 5.8e9"); % Need to append a different frequency for each sweep
%fwrite(NRP67T,"CALC:POW:BURS:AVERAGE")
fwrite(NRP67T,"INIT:ALL");
fwrite(NRP67T,"TRIG:LEV -25 DBM");
pause(5);
fwrite(NRP67T,"TRIG:COUN 20");

message =  sprintf('FETCH%d?',1);
Power = str2num(query(NRP67T, message))
%}
disp("Aye, it ran!")

fclose(NRP67T);