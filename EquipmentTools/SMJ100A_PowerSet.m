function [  ] = SMJ100A_PowerSet( SMJ100A,  Power)
%SMJ100A_PowerSet Sets output power of SMJ100A in dBm
%   N/A
message =  sprintf('POW %f dBm; *OPC?',Power);
[status, Result] = rs_send_query( SMJ100A, message );
    if( ~status ); clear;  return; end

end

