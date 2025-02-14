function [ E4419B ] = E4419B_Setup(Board, Instrument)
%FSQ26_Setup Connects to the E4419B Power Meter and returns an
%object that can access the instrument.

E4419B = gpib('AGILENT', Board, Instrument);
fopen(E4419B)
end

