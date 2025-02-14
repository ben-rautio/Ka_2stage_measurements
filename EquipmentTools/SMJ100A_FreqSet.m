function [ ] = SMJ100A_FreqSet( SMJ100A, Fo )
%SMJ100A_FreqSet Sets the Frequency of the SMJ100A
%   Detailed explanation goes here

message =  sprintf('FREQ %.6f GHz; *OPC?',Fo/1e9);
[status, Result] = rs_send_query( SMJ100A, message );
if( ~status ); clear;  return; end

end

