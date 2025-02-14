function [ SMJ100A ] = SMJ100A_Setup( Board, Instrument )
%SMJ100A_Setup Configures the SMJ100A Signal Generator to be used in MatLab
%   Detailed explanation goes here

%SMJ100A = visa('ni', "GPIB2::28::INSTR");
SMJ100A = gpib('AGILENT', Board, Instrument);
fopen(SMJ100A)

KeepLocalFile    = 0;                       % waveform only temporarily saved
LocalFileName    = 'awgn.wv';               % The local and remote file name
InstrTargetPath = 'D:\';                   % MS Windows based, e.g. SMU, SMJ, SMATE
StartARB         = 1;                       % start in path A
end

