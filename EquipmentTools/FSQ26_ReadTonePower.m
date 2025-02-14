function [ Power ] = FSQ26_ReadTonePower( FSQ26, Freq )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Clear markers
% Trigger the instrument and initiate measurement
 fprintf(FSQ26,'*TRG');
 fprintf(FSQ26,':INITiate:WAVeform');

%Set Marker to Freq
fprintf(FSQ26,['CALC:MARK1:X ' num2str(Freq)]);
%Read Marker
Power = str2num(query(FSQ26, 'CALC:MARK1:Y?'));

end
