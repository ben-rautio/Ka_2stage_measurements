function [power] = NRP_ReadPower(nrp, frequency)
fprintf(nrp, 'SENS:FREQ %f', frequency);
fprintf(nrp, 'INIT:IMM')
tic
measuring  = 1;
while measuring >0
    measuring = str2num(query(nrp, 'STAT:OPER:COND?'));
    if toc >5.0
        break
    end
end
power_watts = str2num( query(nrp, 'FETCH?'));
power = 10.*log10(power_watts.*1000);

% power = str2num( query(nrp, 'FETCH?'));

end