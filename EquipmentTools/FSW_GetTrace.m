function [  ] = FSW_GetTrace(FSW, TraceNumber)
%N6705A_ReadVal Returns the current or voltage on the specified Channel
%   N/A

% message = sprintf('TRAC1:DATA? TRACE1');

message = sprintf('FORM:DEXP:DSEP POIN');
fprintf(FSW,message)

message = sprintf('FORM:DEXP:XDIS STAR');
fprintf(FSW,message)

message = sprintf('FORM:DEXP:FORM CSV');
fprintf(FSW,message)

message = sprintf('FORM:DEXP:HEAD ON');
fprintf(FSW,message)

message = sprintf('FORM:DEXP:TRAC SING');
fprintf(FSW,message)

message = sprintf('MMEM:STOR1:TRAC 1,''C:\\TESTTRACE%u.CSV''', TraceNumber);
fprintf(FSW,message)

% message = sprintf("MMEM:STOR1:TRAC 1, %s",'C:\TESTTRACE3.CSV');
% Value = str2num(query(FSW, message));


%
% -> FORM:DEXP:DSEP POIN
% -> FORM:DEXP:XDIS STAR
% -> FORM:DEXP:FORM CSV
% -> FORM:DEXP:HEAD ON
% -> FORM:DEXP:TRAC SING
% -> MMEM:STOR1:TRACE1,'TRACETEST.CSV'
% MMEM:STOR1:TRAC 1,'TESTTRACE2.CSV'
% MMEM:STOR1:TRAC 1,'C:\TESTTRACE2.CSV'
end

