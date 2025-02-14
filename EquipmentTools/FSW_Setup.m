function [ FSW ] = FSW_Setup( Board, Instrument)
%N6705A_Setup Connects to the N6705A DC Power Supply and returns an
%object that can access the instrument.
FSW = gpib('AGILENT', Board, Instrument);
fopen(FSW)

end

