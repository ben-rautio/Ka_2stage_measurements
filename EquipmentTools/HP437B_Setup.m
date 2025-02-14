function [ HP437B ] = HP437B_Setup(Board, Instrument)
%HP437B_SETUP Summary of this function goes here
%   Detailed explanation goes here
HP437B = gpib('KEYSIGHT', Board, Instrument);
fopen(HP437B)
end

