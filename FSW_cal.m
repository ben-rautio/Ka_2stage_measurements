%% CW Measurement over Frequency
close all; clear all; clc
%% External Tools
dirPath = pwd;
addpath(strcat(dirPath,'\EquipmentTools'))
addpath(strcat(dirPath,'\RSTools'))

%% load cal data
avg_inp_off_s = load("MIWAVE_input_cal.mat","avg_inp_off");
avg_inp_off = avg_inp_off_s.avg_inp_off;
avg_output_atten_s = load("outputCal_data1_raw.mat","avg_output_atten");
avg_output_atten = avg_output_atten_s.avg_output_atten;
avg_GSG_atten_s = load("GSG_cal.mat","avg_gsg_loss");
avg_gsg_atten = avg_GSG_atten_s.avg_gsg_loss;

%calculate input offset to get input power
inp_offset = avg_inp_off - avg_gsg_atten; %add this to inp power meter to get inp power
%calculate output attenuation to get output power
output_atten = avg_output_atten + avg_gsg_atten; %add this to output power meter to get output power
%avg the last 5-10 samples to get a single value vs power
inp_off_avg = mean(inp_offset(size(inp_offset,2)-5:end,:),1);
output_off_avg = mean(output_atten(size(output_atten,2)-5:end,:),1);

%% sweep params
Fo = 38e9:0.2e9:40e9; %Center Frequency (Hz)
Pdes = -35:1:-18; %Power Applied to PA (dBm)
Pavs = repmat(Pdes, length(Fo), 1)';

inp_off_avg = inp_off_avg(inp_off_avg~=0);
output_off_avg = output_off_avg(output_off_avg~=0);

%% add in estimated adapter losses
inp_off_avg = inp_off_avg + 0.1;
output_off_avg = output_off_avg - 0.17;

% Scale these as necessary
Attenuation_scaled = repmat(output_off_avg, length(Pdes), 1); % Output Attenuation, adapter loss estimate
InputOffsets_scaled = repmat(inp_off_avg, length(Pdes), 1); % Input Offsets, adapter loss estimate

%% Prepare Instruments
NRP67T = NRP_Setup_67T_USB(); %coupler meter
% power67 = NRP_ReadPower(NRP67T,40e9)
% NRP_Close( NRP67T )
NRX = NRP67T_Setup(); %output meter for output atten, wait ~5s for read after changing power
% powerNRX = NRP67T_ReadPower(NRX,40e9)
% NRP67T_Close( NRX )

% N6705A Power Supply Settings
Gate1 = 1; %channel
Drain1 = 2; %channel
Gate2 = 3; %channel
Drain2 = 4; %channel
N6705A = N6705A_Setup(7,5); %DC Power Supply

FSW = FSW_Setup(7,22);

%% Prepare Engineering Datasets

size_eng = [length(Pavs), length(Fo)];

FundTone_In = zeros(size_eng);
FundTone_Out = zeros(size_eng);
FundTone_InFreq = zeros(size_eng);
FundTone_OutFreq = zeros(size_eng);
FSW_PM_diff1 = zeros(size_eng);
FSW_PM_diff2 = zeros(size_eng);

%% Measurement Time
%NOTE: needed to add a C7238 Adapter on the FSW to attach to the 2.4mm
%cable
idx = 0; %temp value for recording progress on waitbar
h = waitbar(0, 'Testing Sweep in Progress...');

message = sprintf('INST:SEL AMPL'); % select the amplifier channel
fprintf(FSW,message);
% message = sprintf(['DISP:WIND:TRAC:Y:SCAL:RLEV ', num2str(RLev_amp(i))]); %
% fprintf(FSW,message)
message = sprintf('CONF:DDPD:APPL:STAT OFF'); % Turn off DPD
fprintf(FSW,message);
pause(2) % wait a second
message = sprintf('INIT:CONT ON'); % set amplifier channel to continuous mode
fprintf(FSW,message);
pause(2) % wait a second

%plotting stuff
tiledlayout(1,4)
nexttile
hold on
xlabel("SMW Power")
ylabel("Input Coupled Raw, Output FSW Raw")
hold off
nexttile
hold on
xlabel("Freq")
ylabel("Input Coupled Raw, Output FSW Raw")
hold off
nexttile
hold on
xlabel("SMW Power")
ylabel("Offset between expected PM and FSW reading")
hold off
nexttile
hold on
xlabel("Freq")
ylabel("Offset between expected PM and FSW reading")
hold off

for j = 1:numel(Fo)
    for i = 1:1:length(Pavs)
            
        waitbar(idx./prod(size_eng),h,'Testing Sweep in Progress...')

        message = sprintf('INST:SEL AMPL'); % select the amplifier channel
        fprintf(FSW,message);
        message = sprintf(['SENS:FREQ:CENT ',num2str(Fo(j))]); %set freq
        fprintf(FSW,message);
        pause(2)
        message = sprintf('CONF:GEN:POW:LEV %f;*WAI',Pavs(i)); % set the power
        fprintf(FSW,message);
        pause(2) % wait a second
        message = sprintf('CONF:GEN:RFO:STAT ON;*WAI'); % turn on RF
        fprintf(FSW,message);
        pause(2) % wait a second
        
        FundTone_In(i,j) = NRP_ReadPower(NRP67T, Fo(j)); %Get input power
        FundTone_InFreq(i,j) = NRP_ReadPower(NRP67T, Fo(j)); %Get input power 2nd sample

        message = sprintf('INST:SEL SAN'); % hop back to the spectrum
        fprintf(FSW,message);
        pause(2) % wait a second
        message = sprintf('CALC1:MARK1:MAX'); %move marker to peak power of CW tone
        fprintf(FSW,message);
        message = sprintf('CALC1:MARK1:Y?'); %Grab the power level at the peak power
        FundTone_Out(i,j) = str2double(query(FSW,message));

        %get 2nd sample
        message = sprintf('CALC1:MARK1:MAX'); %move marker to peak power of CW tone
        fprintf(FSW,message);
        message = sprintf('CALC1:MARK1:Y?'); %Grab the power level at the peak power
        FundTone_OutFreq(i,j) = str2double(query(FSW,message));

        %plotting
        nexttile(1)
        hold on
        scatter(Pavs(i,j), FundTone_In(i,j),'red');
        scatter(Pavs(i,j), FundTone_InFreq(i,j),'red');
        scatter(Pavs(i,j), FundTone_Out(i,j),'blue');
        scatter(Pavs(i,j), FundTone_OutFreq(i,j),'blue');
        grid on
        hold off
        nexttile(2)
        hold on
        scatter(Fo(j), FundTone_In(i,j),'red');
        scatter(Fo(j), FundTone_InFreq(i,j),'red');
        scatter(Fo(j), FundTone_Out(i,j),'blue');
        scatter(Fo(j), FundTone_OutFreq(i,j),'blue');
        grid on
        hold off
        
        %calculate difference between expected output power reading and the
        %one measured by the FSW
        expected1 = FundTone_In(i,j) + inp_off_avg(j) - output_off_avg(j); 
        expected2 = FundTone_InFreq(i,j) + inp_off_avg(j) - output_off_avg(j);
        FSW_PM_diff1(i,j) = FundTone_Out(i,j) - expected1;
        FSW_PM_diff2(i,j) = FundTone_OutFreq(i,j) - expected2;

        nexttile(3)
        hold on
        scatter(Pavs(i,j), FSW_PM_diff1(i,j),'red');
        scatter(Pavs(i,j), FSW_PM_diff2(i,j),'blue');
        grid on
        hold off
        nexttile(4)
        hold on
        scatter(Fo(j), FSW_PM_diff1(i,j),'red');
        scatter(Fo(j), FSW_PM_diff2(i,j),'blue');
        grid on
        hold off

        waitbar(idx./prod(size_eng),h,'Testing Sweep in Progress...')
        idx = idx + 1;
        
    end
end

message = sprintf('INST:SEL AMPL'); % select the amplifier channel
fprintf(FSW,message);
message = sprintf('CONF:GEN:POW:LEV %f;*WAI',-60); % set the power to something low
fprintf(FSW,message);
message = sprintf('CONF:GEN:RFO:STAT OFF;*WAI'); %turn off RF
fprintf(FSW,message);

close(h) %close waitbar
            
%% Safe Bench Equipment
NRP_Close( NRP67T ) % coupler power meter
NRP67T_Close( NRX ) %output power meter
% SMJ100A_Close( SMJ100A ) % Vector Signal Generator
FSW_Close(FSW)
% SMJ100A_Close(SMJ100A) %Vector Signal Generator

%% Plotting
close all
tiledlayout(1,6)
%predicted PM output plot vs. what is read out by the FSW
nexttile
hold on
xlabel("SMW Power")
ylabel("Expected Power Meter reading and FSW reading")
inp_off_avg_scaled = repmat(inp_off_avg,numel(Pdes),1);
output_off_avg_scaled = repmat(output_off_avg,numel(Pdes),1);
PM1 = FundTone_In + inp_off_avg_scaled - output_off_avg_scaled;
PM2 = FundTone_InFreq + inp_off_avg_scaled - output_off_avg_scaled;
%predicted output power reading
plot(Pavs,PM1,'r','LineWidth',2)
plot(Pavs,PM2,'b','LineWidth',2)
%FSW reading
plot(Pavs,FundTone_Out,'--r','LineWidth',2)
plot(Pavs,FundTone_OutFreq,'--b','LineWidth',2)
grid on
hold off
nexttile
hold on
xlabel("Freq")
ylabel("Expected Power Meter reading and FSW reading")
for i = 1:numel(Pdes)
    plot(Fo,PM1(i,:),'r','LineWidth',2)
    plot(Fo,PM1(i,:),'b','LineWidth',2)
    plot(Fo, FundTone_Out(i,:),'--r','LineWidth',2)
    plot(Fo, FundTone_OutFreq(i,:),'--b','LineWidth',2)
end
grid on
hold off
%plot offsets again
nexttile
hold on
xlabel("SMW Power")
ylabel("Offset BTW PM and FSW")
plot(Pavs,FSW_PM_diff1,'r','LineWidth',2)
plot(Pavs,FSW_PM_diff2,'b','LineWidth',2)
grid on
hold off
%vs freq
nexttile
hold on
xlabel("Freq")
ylabel("Offset between PM and FSW")
for i = 1:numel(Pdes)
    plot(Fo,FSW_PM_diff1(i,:),'r','LineWidth',2)
    plot(Fo,FSW_PM_diff2(i,:),'b','LineWidth',2)
end
grid on
hold off
% something weird is happening at 38GHz, PAVS = -19dBm, don't take avg of
% this sample!
fsw_off_avg = (FSW_PM_diff2 + FSW_PM_diff2) ./ 2;
fsw_off_avg = mean(fsw_off_avg);
special_case = (FSW_PM_diff2(:,1) + FSW_PM_diff2(:,1)) ./ 2;
special_case(end-1) = [];
special_case = mean(special_case);
fsw_off_avg(1) = special_case;
fsw_off_avg = abs(fsw_off_avg);
%compare predicted output power of PM to FSW power with offset applied
fsw_off_avg_rep = repmat(fsw_off_avg,numel(Pdes),1);
fsw_pred_pout1 = FundTone_Out + fsw_off_avg_rep;
fsw_pred_pout2 = FundTone_OutFreq + fsw_off_avg_rep;
error1 = PM1 - fsw_pred_pout1;
error2 = PM2 - fsw_pred_pout2;
%vs Pdes
nexttile
hold on
xlabel("SMW Power")
ylabel("Error")
plot(Pavs,error1,'r','LineWidth',2)
plot(Pavs,error2,'b','LineWidth',2)
grid on
hold off
%vs freq
nexttile
hold on
xlabel("Freq")
ylabel("Error")
for i = 1:numel(Pdes)
    plot(Fo,error1(i,:),'r','LineWidth',2)
    plot(Fo,error2(i,:),'b','LineWidth',2)
end
grid on
hold off

