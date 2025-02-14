function [ N6705A ] = N6705A_Setup( Board, Instrument)
%N6705A_Setup Connects to the N6705A DC Power Supply and returns an
%object that can access the instrument.
N6705A = gpib('AGILENT', Board, Instrument);
fopen(N6705A)

end

