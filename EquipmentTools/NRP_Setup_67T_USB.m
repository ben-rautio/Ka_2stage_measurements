function [nrp] = NRP_Setup_67T_USB()
% nrp = visa('keysight',"USB::0x0aad::0x0145::102483::INSTR") % original
% nrp = visa('keysight',"USB0::0x0AAD::0x0138::102270::0::INSTR")
% nrp = visa('keysight',"USB0::0x0AAD::0x0145::102483::0::INSTR")
nrp = visa('keysight',"USB0::0x0AAD::0x0158::101117::0::INSTR") %for the 67T (currently on input)

fopen(nrp)
fprintf(nrp, 'INIT:CONT OFF');
pause(2)
fprintf(nrp, 'SENS:FUNC "POW:AVG"');
pause(2)
fprintf(nrp, 'SENS:AVER:COUN:AUTO OFF');
pause(2)
fprintf(nrp, 'SENS:AVER:COUN 16');
pause(2)
fprintf(nrp, 'SENS:AVER:STAT ON');
pause(2)
fprintf(nrp, 'SENS:AVER:TCON REP');
pause(2)
fprintf(nrp, 'FORMAT ASCII');
end