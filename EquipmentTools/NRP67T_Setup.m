function [ NRP67T ] = NRP67T_Setup()
%NRP67T_Setup Connects to the NRP67T Power Sensor (NRX Power Meter) and returns an
%object that can access the instrument.

% NRP67T = visa('ni', "TCPIP::192.168.0.10::INSTR"); %
NRP67T = visa('keysight', "TCPIP::192.168.0.10::INSTR"); %
% NRP67T = visa('keysight', "USB0::0x0AAD::0x0158::101117::0::INSTR"); %
fopen(NRP67T);

% Says whether or not it connected to the correct power sensor
dev_name = query(NRP67T, "*IDN?") % 
    if ~contains(dev_name, "Rohde&Schwarz")
        disp("Ah crap it connected to the wrong thing. Ooops"); 
        return;
    else
        disp("Ay it connected!");
    end
end