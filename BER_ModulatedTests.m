%% CW Measurement over Frequency

close all; clear all; clc
%% External Tools
dirPath = pwd;
addpath(strcat(dirPath,'\EquipmentTools'))
addpath(strcat(dirPath,'\RSTools'))

%% Load data
avg_inp_off_s = load("MIWAVE_input_cal.mat","avg_inp_off");
avg_inp_off = avg_inp_off_s.avg_inp_off;
avg_output_atten_s = load("outputCal_data1_raw.mat","avg_output_atten");
avg_output_atten = avg_output_atten_s.avg_output_atten;
avg_GSG_atten_s = load("GSG_cal.mat","avg_gsg_loss");
avg_gsg_atten = avg_GSG_atten_s.avg_gsg_loss;
fsw_avg_offset_s = load("FSW_cal_offset.mat","fsw_off_avg");
fsw_avg_offset = fsw_avg_offset_s.fsw_off_avg; %add this to FSW reading with output atten to get output power

%calculate input offset to get input power
inp_offset = avg_inp_off - avg_gsg_atten; %add this to inp power meter to get inp power
%calculate output attenuation to get output power
output_atten = avg_output_atten + avg_gsg_atten; %add this to output power meter to get output power
%avg the last 5-10 samples to get a single value vs power
inp_off_avg = mean(inp_offset(size(inp_offset,2)-5:end,:),1);
inp_off_avg = inp_off_avg(inp_off_avg~=0);
output_off_avg = mean(output_atten(size(output_atten,2)-5:end,:),1);
output_off_avg = output_off_avg(output_off_avg~=0);

%% Prepare Sweeps
Fo = 38e9:0.2e9:40e9; %Center Frequency (Hz)

Pdes = -33:1:-33; %Power set on signal generator (dBm)
Pavs = repmat(Pdes, length(Fo), 1)'; % shaping

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

Vd1_eng = zeros(size_eng);
Vg1_eng = zeros(size_eng);
Id1_eng = zeros(size_eng);
Ig1_eng = zeros(size_eng);
Vd2_eng = zeros(size_eng);
Vg2_eng = zeros(size_eng);
Id2_eng = zeros(size_eng);
Ig2_eng = zeros(size_eng);

Vd1_eng_dpd = zeros(size_eng);
Vg1_eng_dpd = zeros(size_eng);
Id1_eng_dpd = zeros(size_eng);
Ig1_eng_dpd = zeros(size_eng);
Vd2_eng_dpd = zeros(size_eng);
Vg2_eng_dpd = zeros(size_eng);
Id2_eng_dpd = zeros(size_eng);
Ig2_eng_dpd = zeros(size_eng);

FundTone_In = zeros(size_eng);
FundTone_Out = zeros(size_eng);
FundTone_InFreq = zeros(size_eng);
FundTone_OutFreq = zeros(size_eng);

FundTone_In_dpd = zeros(size_eng);
FundTone_Out_dpd = zeros(size_eng);
FundTone_InFreq_dpd = zeros(size_eng);
FundTone_OutFreq_dpd = zeros(size_eng);

adjpow_upper = zeros(size_eng);
adjpow_lower = zeros(size_eng);
adjpow_upper_dpd = zeros(size_eng);
adjpow_lower_dpd = zeros(size_eng);

EVM = zeros(size_eng);
EVM_dpd = zeros(size_eng);

CFAC_out = zeros(size_eng);
CFAC_out_dpd = zeros(size_eng);

AMAM = zeros(size_eng);
AMAM_dpd = zeros(size_eng);

AMPM = zeros(size_eng);
AMPM_dpd = zeros(size_eng);

%% Measurement Time
close all
% Grab quiescent values before starting sweep
Vgsq1 =  N6705A_GetVal(N6705A, Gate1, 'voltage');
Igsq1 = N6705A_GetVal(N6705A, Gate1, 'current');
Idsq1 = N6705A_GetVal(N6705A, Drain1, 'current');
Vdsq1 = N6705A_GetVal(N6705A, Drain1, 'voltage');
Vgsq2 =  N6705A_GetVal(N6705A, Gate2, 'voltage');
Igsq2 = N6705A_GetVal(N6705A, Gate2, 'current');
Idsq2 = N6705A_GetVal(N6705A, Drain2, 'current');
Vdsq2 = N6705A_GetVal(N6705A, Drain2, 'voltage');

idx = 0; %temp value for recording progress on waitbar
h = waitbar(0, 'Testing Sweep in Progress...');

%plot powers in case...
tiledlayout(2,2)
nexttile
hold on
xlabel("SMW Power")
ylabel("Input Power")
hold off
nexttile
hold on
xlabel("SMW Power")
ylabel("Channel Power")
hold off
nexttile
hold on
xlabel("Center Freq")
ylabel("Input Power")
hold off
nexttile
hold on
xlabel("Center Freq")
ylabel("Output Power")
hold off

for j = 1:1:numel(Fo)
    for i = 1:1:length(Pdes)
        %setup ref levels, generate signal, set freq and power, setup
        %channels
        FSW_ModulatedSetup(FSW,Fo(j),Pavs(i),200e6)
        message = sprintf('INST:SEL AMPL'); % back to the amplifier channel
        fprintf(FSW,message);
        pause(2)
        message = sprintf('INIT:CONT ON'); % set amplifier channel to continuous mode
        fprintf(FSW,message);
        pause(2) % wait a second
        message = sprintf('FETC:MACC:REVM:CURR:RES?'); % get EVM with no DPD
        EVM(i,j) = str2double(query(FSW, message)); % get EVM with no DPD
        pause(2) % wait a second
        message = sprintf('FETC:POW:CFAC:OUT:CURR:RES?'); % GRAB THE OUTPUT CREST FACTOR
        CFAC_out(i,j) = str2double(query(FSW, message));
        pause(2) % wait a second
        message = sprintf('FETC:AMAM:CWID:CURR:RES?'); % GRAB AM/AM
        AMAM(i,j) = str2double(query(FSW, message));
        pause(2) % wait a second
        message = sprintf('FETC:AMPM:CWID:CURR:RES?'); % GRAB AM/PM
        AMPM(i,j) = str2double(query(FSW, message));
        pause(2) % wait a second
        message = sprintf('INST:SEL SAN'); % hop back to the spectrum
        fprintf(FSW,message);
        pause(2) % wait a second
        
        % Grab NDPD values
        FundTone_In(i,j) = NRP_ReadPower(NRP67T, Fo(j));
        FundTone_Out(i,j) = FSW_GetVal(FSW,'Channel');
        FundTone_InFreq(i,j) = NRP_ReadPower(NRP67T, Fo(j));
        FundTone_OutFreq(i,j) = FSW_GetVal(FSW,'Channel');

        Vd1_eng(i,j) = N6705A_GetVal(N6705A, Drain1, 'voltage');
        Vg1_eng(i,j) = N6705A_GetVal(N6705A, Gate1, 'voltage');
        Id1_eng(i,j) = N6705A_GetVal(N6705A, Drain1, 'current');
        Ig1_eng(i,j) = N6705A_GetVal(N6705A, Gate1, 'current');
        Vd2_eng(i,j) = N6705A_GetVal(N6705A, Drain2, 'voltage');
        Vg2_eng(i,j) = N6705A_GetVal(N6705A, Gate2, 'voltage');
        Id2_eng(i,j) = N6705A_GetVal(N6705A, Drain2, 'current');
        Ig2_eng(i,j) = N6705A_GetVal(N6705A, Gate2, 'current');
         
        adjpow_data = FSW_GetVal(FSW, 'Sidebar'); % Ask the FSW now
        adjpow_lower(i,j) = adjpow_data(2);
        adjpow_upper(i,j) = adjpow_data(3);
        
        % DPD Time
        message = sprintf('INST:SEL AMPL'); % back to the amplifier channel
        fprintf(FSW,message);
        %ask ricky, is there a reason why this needs to be reset once DPD
        %was turned on?
        %using default direct dpd settings, may want to adjust this
        %later...
        message = sprintf('CONF:DDPD:STAR;*WAI'); % start DPD
        fprintf(FSW,message);
        pause(15)   % give it time to breathe
        message = sprintf('INIT:CONT ON'); % set amplifier channel to continuous mode
        fprintf(FSW,message);
        pause(2) % wait a second
        message = sprintf('FETC:MACC:REVM:CURR:RES?'); % get EVM with DPD
        EVM_dpd(i,j) = str2double(query(FSW, message)); % get EVM with DPD
        pause(2) % wait a second
        message = sprintf('FETC:POW:CFAC:OUT:CURR:RES?'); % GRAB THE OUTPUT CREST FACTOR
        CFAC_out_dpd(i,j) = str2double(query(FSW, message));
        pause(2) % wait a second
        message = sprintf('FETC:AMAM:CWID:CURR:RES?'); % GRAB AM/AM
        AMAM_dpd(i,j) = str2double(query(FSW, message));
        pause(2) % wait a second
        message = sprintf('FETC:AMPM:CWID:CURR:RES?'); % GRAB AM/PM
        AMPM_dpd(i,j) = str2double(query(FSW, message));
        pause(2) % wait a second
        message = sprintf('INST:SEL SAN'); % back to the Spectrum
        fprintf(FSW,message);
	    pause(2)
        
        % Grab DPD Values
        FundTone_In_dpd(i,j) = NRP_ReadPower(NRP67T, Fo(j));
        FundTone_Out_dpd(i,j) = FSW_GetVal(FSW,'Channel');
        FundTone_InFreq_dpd(i,j) = NRP_ReadPower(NRP67T, Fo(j));
        FundTone_OutFreq_dpd(i,j) = FSW_GetVal(FSW,'Channel');
        
        Vd1_eng_dpd(i,j) = N6705A_GetVal(N6705A, Drain1, 'voltage');
        Vg1_eng_dpd(i,j) = N6705A_GetVal(N6705A, Gate1, 'voltage');
        Id1_eng_dpd(i,j) = N6705A_GetVal(N6705A, Drain1, 'current');
        Ig1_eng_dpd(i,j) = N6705A_GetVal(N6705A, Gate1, 'current');
        Vd2_eng_dpd(i,j) = N6705A_GetVal(N6705A, Drain2, 'voltage');
        Vg2_eng_dpd(i,j) = N6705A_GetVal(N6705A, Gate2, 'voltage');
        Id2_eng_dpd(i,j) = N6705A_GetVal(N6705A, Drain2, 'current');
        Ig2_eng_dpd(i,j) = N6705A_GetVal(N6705A, Gate2, 'current');
        
        adjpow_data = FSW_GetVal(FSW, 'Sidebar'); % Ask the FSW now
        adjpow_lower_dpd(i,j) = adjpow_data(2); %ACPR lower
        adjpow_upper_dpd(i,j) = adjpow_data(3); %ACPR upper
        
        %plot stuff during sweeps
        nexttile(1)
        hold on
        scatter(Pavs(i),FundTone_In(i,j)+inp_off_avg(j),'r')
        scatter(Pavs(i),FundTone_InFreq(i,j)+inp_off_avg(j),'+r')
        scatter(Pavs(i),FundTone_In_dpd(i,j)+inp_off_avg(j),'b')
        scatter(Pavs(i),FundTone_InFreq_dpd(i,j)+inp_off_avg(j),'+b')
        hold off
        nexttile(2)
        hold on
        scatter(Pavs(i),FundTone_Out(i,j)+fsw_avg_offset(j),'r')
        scatter(Pavs(i),FundTone_OutFreq(i,j)+fsw_avg_offset(j),'+r')
        scatter(Pavs(i),FundTone_Out_dpd(i,j)+fsw_avg_offset(j),'b')
        scatter(Pavs(i),FundTone_OutFreq_dpd(i,j)+fsw_avg_offset(j),'+b')
        hold off
        nexttile(3)
        hold on
        scatter(Fo(j),FundTone_In(i,j)+inp_off_avg(j),'r')
        scatter(Fo(j),FundTone_InFreq(i,j)+inp_off_avg(j),'+r')
        scatter(Fo(j),FundTone_In_dpd(i,j)+inp_off_avg(j),'b')
        scatter(Fo(j),FundTone_InFreq_dpd(i,j)+inp_off_avg(j),'+b')
        hold off
        nexttile(4)
        hold on
        scatter(Fo(j),FundTone_Out(i,j)+fsw_avg_offset(j),'r')
        scatter(Fo(j),FundTone_OutFreq(i,j)+fsw_avg_offset(j),'+r')
        scatter(Fo(j),FundTone_Out_dpd(i,j)+fsw_avg_offset(j),'b')
        scatter(Fo(j),FundTone_OutFreq_dpd(i,j)+fsw_avg_offset(j),'+b')
        hold off

        waitbar(idx./prod(size_eng),h,'Testing Sweep in Progress...')
        idx = idx + 1;
    end
end

message = sprintf('INST:SEL AMPL'); % select the amplifier channel
fprintf(FSW,message);
message = sprintf('CONF:GEN:RFO:STAT OFF;*WAI'); %turn off RF
fprintf(FSW,message);
message = sprintf('CONF:GEN:POW:LEV %f;*WAI',-60); % set the power to something low
fprintf(FSW,message);

close(h) %close waitbar
            
%% Safe Bench Equipment
NRP_Close( NRP67T ) % coupler power meter
NRP67T_Close( NRX ) %output power meter
N6705A_Close( N6705A ) % DC Supply (this does not turn it off - it just terminates the remote connection)
FSW_Close(FSW)

%% Calculate Pout, Gain, PAE no DPD
% Scale these as necessary
Attenuation_scaled = repmat(output_off_avg, length(Pdes), 1); % Output Attenuation vs freq
InputOffsets_scaled = repmat(inp_off_avg, length(Pdes), 1); % Input Offsets vs freq
FSWOffsets_scaled = repmat(fsw_avg_offset,length(Pdes), 1); %FSW offsets vs freq

Pout1 = FundTone_Out + Attenuation_scaled + FSWOffsets_scaled;
Pout1_watts = 10.^(Pout1./10)./1000;
Pout2 = FundTone_OutFreq + Attenuation_scaled + FSWOffsets_scaled;
Pout2_watts = 10.^(Pout2./10)./1000;

Pin1 = FundTone_In + InputOffsets_scaled;
Pin1_watts = 10.^(Pin1./10)./1000;
Pin2 = FundTone_InFreq + InputOffsets_scaled;
Pin2_watts = 10.^(Pin2./10)./1000;

Pdc_Drain1 = Vd1_eng.*Id1_eng;
Pdc_Gate1 = Vg1_eng.*Ig1_eng;
Pdc_tot1 = Pdc_Drain1 + Pdc_Gate1;
Pdc_Drain2 = Vd2_eng.*Id2_eng;
Pdc_Gate2 = Vg2_eng.*Ig2_eng;
Pdc_tot2 = Pdc_Drain2 + Pdc_Gate2;
Pdc_tot = Pdc_tot1 + Pdc_tot2;

DE1 = (Pout_watts1./Pdc_tot).*100;
PAE1 = ((Pout_watts1-Pin_watts1)./Pdc_tot).*100;
Gain1 = Pout1 - Pin1;

DE2 = (Pout_watts2./Pdc_tot).*100;
PAE2 = ((Pout_watts2-Pin_watts2)./Pdc_tot).*100;
Gain2 = Pout2 - Pin2;

%% Calculate Pout, Gain, PAE with DPD
Pout1_dpd = FundTone_Out_dpd + Attenuation_scaled + FSWOffsets_scaled;
Pout1_watts_dpd = 10.^(Pout1_dpd./10)./1000;
Pout2_dpd = FundTone_OutFreq_dpd + Attenuation_scaled + FSWOffsets_scaled;
Pout2_watts_dpd = 10.^(Pout2_dpd./10)./1000;

Pin1_dpd = FundTone_In_dpd + InputOffsets_scaled;
Pin1_watts_dpd = 10.^(Pin1_dpd./10)./1000;
Pin2_dpd = FundTone_InFreq_dpd + InputOffsets_scaled;
Pin2_watts_dpd = 10.^(Pin2_dpd./10)./1000;

DE1_dpd = (Pout1_watts_dpd./Pdc_tot).*100;
PAE1_dpd = ((Pout1_watts_dpd-Pin1_watts_dpd)./Pdc_tot).*100;
Gain1_dpd = Pout1_dpd - Pin1_dpd;

DE2_dpd = (Pout2_watts_dpd./Pdc_tot).*100;
PAE2_dpd = ((Pout2_watts_dpd-Pin2_watts_dpd)./Pdc_tot).*100;
Gain2_dpd = Pout2_dpd - Pin2_dpd;
%% Create Data Summary
%quiescent points
data.Vgsq1 = Vgsq1;
data.Igsq1 = Igsq1;
data.Vdsq1 = Vdsq1;
data.Idsq1 = Idsq1;
data.Vgsq2 = Vgsq2;
data.Igsq2 = Igsq2;
data.Vdsq2 = Vdsq2;
data.Idsq2 = Idsq2;
%powers
data.Pavs = Pavs;
data.Pin1 = Pin1;
data.Pin2 = Pin2;
data.Pin1_dpd = Pin1_dpd;
data.Pin2_dpd = Pin2_dpd;
data.Gain1 = Gain1;
data.Gain2 = Gain2;
data.Gain1_dpd = Gain1_dpd;
data.Gain2_dpd = Gain2_dpd;
data.Pout1 = Pout1;
data.Pout2 = Pout2;
data.Pout1_dpd = Pout1_dpd;
data.Pout2_dpd = Pout2_dpd;
data.PAE1 = PAE1;
data.PAE2 = PAE2;
data.PAE1_dpd = PAE1_dpd;
data.PAE2_dpd = PAE2_dpd;
data.Fund_in = FundTone_In;
data.Fund_in_dpd = FundTone_In_dpd;
data.Fund_out = FundTone_Out;
data.Fund_out_dpd = FundTone_Out_dpd;
%gate/drain voltages/currents
data.Vgate1 = Vg1_eng;
data.Vgate1_dpd = Vg1_eng_dpd;
data.Igate1 = Ig1_eng;
data.Igate1_dpd = Ig1_eng_dpd;
data.Vdrain1 = Vd1_eng;
data.Vdrain1_dpd = Vd1_eng_dpd;
data.Idrain1 = Id1_eng;
data.Idrain1_dpd = Id1_eng_dpd;
data.Vgate2 = Vg2_eng;
data.Vgate2_dpd = Vg2_eng_dpd;
data.Igate2 = Ig2_eng;
data.Igate2_dpd = Ig2_eng_dpd;
data.Vdrain2 = Vd2_eng;
data.Vdrain2_dpd = Vd2_eng_dpd;
data.Idrain2 = Id2_eng;
data.Idrain2_dpd = Id2_eng_dpd;
%modulated performance
data.EVM = EVM;
data.EVM_dpd = EVM_dpd;
data.CFAC_out = CFAC_out;
data.CFAC_out_dpd = CFAC_out_dpd;
data.AMAM = AMAM;
data.AMAM_dpd = AMAM_dpd;
data.AMPM = AMPM;
data.AMPM_dpd = AMPM_dpd;
data.ACPR_Lower = adjpow_lower;
data.ACPR_Upper = adjpow_upper;
data.ACPR_Lower_dpd = adjpow_lower_dpd;
data.ACPR_Upper_dpd = adjpow_upper_dpd;
%% Plotting

set(0,'DefaultLineLineWidth', 2)

figure(11)
set(gcf, 'color', 'w');
hold on 

yyaxis left
plot(Pout, Gain)
plot(Pout_dpd, Gain_dpd)
plot(Pout, CFAC_out)
plot(Pout_dpd, CFAC_out_dpd)
xlabel('Pout (dBm)')
ylabel('Gain (dB)')
title('20 MHz, 10 dB PAPR Modulated Driveup')

yyaxis right
plot(Pout, PAE)
plot(Pout_dpd, PAE_dpd)
ylabel('Efficiency (%)')
legend('Gain no DPD','Gain DPD','Crest Factor Out no DPD','Crest Factor Out DPD','PAE no DPD','PAE DPD')

grid on
ax = gca;
ax.GridAlphaMode = 'manual';
ax.GridAlpha = 0.35;
ax.YColor = 'k';

hold off

figure(22)
set(gcf, 'color', 'w');
hold on 


plot(Pout,adjpow(:,2))
plot(Pout_dpd,adjpow_dpd(:,2))


legend('ACPR_L No DPD (dBc)','ACPR_L DPD (dBc)')
xlabel('Pout (dBm)')
ylabel('ACPR (dBc)')
title('20 MHz, 10 dB PAPR Modulated Driveup')
grid on
hold off

figure(23)
set(gcf, 'color', 'w');
hold on 
grid on

plot(Pout,EVM)
plot(Pout_dpd,EVM_dpd)

legend('EVM No DPD (%)','EVM DPD (%)')
xlabel('Pout (dBm)')
ylabel('EVM (%)')
title('20 MHz, 10 dB PAPR Modulated Driveup')
hold off

%% test