function [ FSQ26 ] = FSQ26_Setup( centerFrequency)%, bandwidth, measureTime,...
    %mechAttenuation )
%FSQ26_Setup Connects to the FSQ26 Vector Signal Generator and returns an
%object that can access the instrument.
%   centerFrequency: Initial Center Frequency of Operation(Hz)
%   bandwidth: Resolution bandwidth of signal(Hz)
%   measureTime: Sets the integration time(sec)
%   mechAttenuation: Mechanical attenuation Setting(dB)

% Bandwidth of the signal (Hz)
bandwidth = 25e6;
% Measurement time (s)
measureTime = 8e-3;
% Mechanical attenuation in the signal analyzer(dB)
mechAttenuation = 0;


FSQ26 = gpib('AGILENT', 7, 20);
FSQ26.InputBufferSize = 30e6;
FSQ26.Timeout = 20;
fopen(FSQ26)

% Set up signal analyzer mode to Basic/IQ mode
%fprintf(FSQ26,':INSTrument:SELect BASIC');

% Set the center frequency
fprintf(FSQ26,[':SENSe:FREQuency:CENTer ' num2str(centerFrequency)]);

% Set the resolution bandwidth
%fprintf(FSQ26,[':SENSe:WAVEform:BANDwidth:RESolution ' num2str(bandwidth)]);

% Turn off averaging
fprintf(FSQ26,':SENSe:WAVeform:AVER OFF');

% set to take one single measurement once the trigger line goes high
fprintf(FSQ26,':INIT:CONT OFF');

% Set the trigger to external source 1 with positive slope triggering
fprintf(FSQ26,':TRIGger:WAVeform:SOURce IMMediate');
fprintf(FSQ26,':TRIGger:LINE:SLOPe POSitive');

% Set the time for which measurement needs to be made
%fprintf(FSQ26,[':WAVeform:SWE:TIME '  num2str(measureTime)]);

% Turn off electrical attenuation.
fprintf(FSQ26,':SENSe:POWer:RF:EATTenuation:STATe OFF');

% Set mechanical attenuation level
fprintf(FSQ26,[':SENSe:POWer:RF:ATTenuation ' num2str(mechAttenuation)]);

% Turn IQ signal ranging to auto
%fprintf(FSQ26,':SENSe:VOLTage:IQ:RANGe:AUTO ON');

% Set the endianness of returned data
fprintf(FSQ26,':FORMat:BORDer NORMal');

% Set the format of the returned data
fprintf(FSQ26,':FORMat:DATA REAL,32');

end

