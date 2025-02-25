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

%% Load data
% load('10-3-24_full_cal.mat');

% Cal = 10;

% Cal_Fo = Fo;

% Cal Freq Index:
% 1 = 2.50 GHz
% 2 = 2.55 GHz
% 3 = 2.60 GHz
% 4 = 2.65 GHz
% 5 = 2.70 GHz
% 6 = 2.75 GHz
% 7 = 2.80 GHz
% 8 = 2.85 GHz
% 9 = 2.90 GHz
% 10 = 2.95 GHz
% 11 = 3.00 GHz
% 12 = 3.05 GHz
% 13 = 3.10 GHz

% Fo = Cal_Fo(Cal); %Center Frequency (Hz)

Fo = 38e9:0.2e9:40e9; %Center Frequency (Hz)


%% Prepare Sweeps

% Pdes = -40:1.0:-20; %Power Applied to UPA (dBm)
% RLev_amp = linspace(-10,10,length(Pdes)); % Spectrum analyzer amp channel reference level range UPA
% RLev_san = linspace(-30,10,length(Pdes)); % Spectrum analyzer san channel reference level range UPA

Pdes = -35:1:-18; %Power Applied to PA (dBm)
RLev_amp = linspace(-15,10,length(Pdes)); % Spectrum analyzer amp channel reference level range PA
RLev_san = linspace(-40,10,length(Pdes)); % Spectrum analyzer san channel reference level range PA


Pavs = repmat(Pdes, length(Fo), 1)'; 


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

FSW = FSW_Setup(9,22);

%% Prepare Engineering Datasets

size_eng = [length(Pavs), length(Fo)];

Vd1_eng = zeros(size_eng);
Vg1_eng = zeros(size_eng);
Id1_eng = zeros(size_eng);
Ig1_eng = zeros(size_eng);

Vd1_eng_dpd = zeros(size_eng);
Vg1_eng_dpd = zeros(size_eng);
Id1_eng_dpd = zeros(size_eng);
Ig1_eng_dpd = zeros(size_eng);

FundTone_In = zeros(size_eng);
FundTone_Out = zeros(size_eng);
FundTone_InFreq = zeros(size_eng);
FundTone_OutFreq = zeros(size_eng);

FundTone_In_dpd = zeros(size_eng);
FundTone_Out_dpd = zeros(size_eng);
FundTone_InFreq_dpd = zeros(size_eng);
FundTone_OutFreq_dpd = zeros(size_eng);

adjpow = zeros(size_eng(1) , 3);
adjpow_dpd = zeros(size_eng(1) , 3);

EVM = zeros(size_eng);
EVM_dpd = zeros(size_eng);

CFAC_out = zeros(size_eng);
CFAC_out_dpd = zeros(size_eng);

AMAM = zeros(size_eng);
AMAM_dpd = zeros(size_eng);

AMPM = zeros(size_eng);
AMPM_dpd = zeros(size_eng);

%% Measurement Time
%NOTE: needed to add a C7238 Adapter on the FSW to attach to the 2.4mm
%cable
idx = 0; %temp value for recording progress on waitbar
h = waitbar(0, 'Testing Sweep in Progress...')

for j = 1:numel(Fo)
    for i = 1:1:length(Pavs)
            
        waitbar(idx./prod(size_eng),h,'Testing Sweep in Progress...')
        
        message = sprintf('INST:SEL AMPL'); % select the amplifier channel
        fprintf(FSW,message)
        message = sprintf(['DISP:WIND:TRAC:Y:SCAL:RLEV ', num2str(RLev_amp(i))]); %
        fprintf(FSW,message)
        message = sprintf('CONF:DDPD:APPL:STAT OFF'); % Turn off DPD
        fprintf(FSW,message)
        pause(2) % wait a second
    
        message = sprintf(['SENS:FREQ:CENT ',num2str(Fo(j))]); %set freq
        fprintf(FSW,message)
        pause(2)
    
        message = sprintf('CONF:GEN:POW:LEV %f;*WAI',Pavs(i)); % set the power
        fprintf(FSW,message) 
        pause(2) % wait a second
        message = sprintf('CONF:GEN:RFO:STAT ON;*WAI'); % turn on RF
        fprintf(FSW,message)
        pause(5) % wait a second
        message = sprintf('INST:SEL SAN'); % back to the Spectrum
        fprintf(FSW,message)
        message = sprintf(['DISP:WIND:TRAC:Y:SCAL:RLEV ', num2str(RLev_san(i))]); %
        fprintf(FSW,message)
	    pause(5)
        message = sprintf('INST:SEL AMPL'); % back to the amplifier channel
        fprintf(FSW,message)
        pause(5)
        message = sprintf('INIT:CONT ON'); % set amplifier channel to continuous mode
        fprintf(FSW,message)
        pause(2) % wait a second

        message = sprintf('CALC1:MARK1:MAX')

        % message = sprintf('FETC:MACC:REVM:CURR:RES?'); % get EVM with no DPD
        % EVM(i) = str2num(query(FSW, message)); % get EVM with no DPD
        % pause(2) % wait a second
        % message = sprintf('FETC:POW:CFAC:OUT:CURR:RES?'); % GRAB THE OUTPUT CREST FACTOR
        % CFAC_out(i) = str2num(query(FSW, message));
        % pause(2) % wait a second
        % message = sprintf('FETC:AMAM:CWID:CURR:RES?'); % GRAB AM/AM
        % AMAM(i) = str2num(query(FSW, message));
        % pause(2) % wait a second
        % message = sprintf('FETC:AMPM:CWID:CURR:RES?'); % GRAB AM/PM
        % AMPM(i) = str2num(query(FSW, message));
        % pause(2) % wait a second
        % message = sprintf('INST:SEL SAN'); % hope back to the spectrum
        % fprintf(FSW,message)
        % pause(2) % wait a second
        
        % Grab NDPD values
        % FundTone_In(i) = NRP_ReadPower(NRP18S, Fo);
        % FundTone_InFreq(i) = NRP_ReadPower(NRP18S, Fo);
        % FundTone_Out(i) = NRP_ReadPower(NRP18SN, Fo);
        % 
        % Vg1_eng(i) = E36232A_GetVal(gate_supply, 'voltage');
        % Vd1_eng(i) = E36232A_GetVal(drain_supply, 'voltage');
        % Ig1_eng(i) = E36232A_GetVal(gate_supply, 'current');
        % Id1_eng(i) = E36232A_GetVal(drain_supply, 'current');
        % 
        % adjpow(i,:) = FSW_GetVal(FSW, 'Sidebar'); % Ask the FSW now
        % 
        % % DPD Time
        % message = sprintf('INST:SEL AMPL'); % back to the amplifier channel
        % fprintf(FSW,message)
        % message = sprintf(['DISP:WIND:TRAC:Y:SCAL:RLEV ', num2str(RLev_amp(i))]); %
        % fprintf(FSW,message)
        % message = sprintf('CONF:DDPD:STAR;*WAI'); % start DPD
        % fprintf(FSW,message)
        % pause(15)   % give it time to breathe
        % message = sprintf('INST:SEL SAN'); % back to the Spectrum
        % fprintf(FSW,message)
        % message = sprintf(['DISP:WIND:TRAC:Y:SCAL:RLEV ', num2str(RLev_san(i))]); %
        % fprintf(FSW,message)
	    % pause(5)
        % message = sprintf('INST:SEL AMPL'); % back to the amplifier channel
        % fprintf(FSW,message)
        % pause(5)
        % message = sprintf('INIT:CONT ON'); % set amplifier channel to continuous mode
        % fprintf(FSW,message)
        % pause(2) % wait a second
        % message = sprintf('FETC:MACC:REVM:CURR:RES?'); % get EVM with DPD
        % EVM_dpd(i) = str2num(query(FSW, message)); % get EVM with DPD
        % pause(2) % wait a second
        % message = sprintf('FETC:POW:CFAC:OUT:CURR:RES?'); % GRAB THE OUTPUT CREST FACTOR
        % CFAC_out_dpd(i) = str2num(query(FSW, message));
        % pause(2) % wait a second
        % message = sprintf('FETC:AMAM:CWID:CURR:RES?'); % GRAB AM/AM
        % AMAM_dpd(i) = str2num(query(FSW, message));
        % pause(2) % wait a second
        % message = sprintf('FETC:AMPM:CWID:CURR:RES?'); % GRAB AM/PM
        % AMPM_dpd(i) = str2num(query(FSW, message));
        % pause(2) % wait a second
        % message = sprintf('INST:SEL SAN'); % back to the Spectrum
        % fprintf(FSW,message)
	    % pause(5)
        % 
        % % Grab DPD Values
        % FundTone_In_dpd(i) = NRP_ReadPower(NRP18S, Fo);
        % FundTone_InFreq_dpd(i) = NRP_ReadPower(NRP18S, Fo);
        % FundTone_Out_dpd(i) = NRP_ReadPower(NRP18SN, Fo);
        % 
        % Vg1_eng_dpd(i) = E36232A_GetVal(gate_supply, 'voltage');
        % Vd1_eng_dpd(i) = E36232A_GetVal(drain_supply, 'voltage');
        % Ig1_eng_dpd(i) = E36232A_GetVal(gate_supply, 'current');
        % Id1_eng_dpd(i) = E36232A_GetVal(drain_supply, 'current');
        
        % adjpow_dpd(i,:) = FSW_GetVal(FSW, 'Sidebar');
        
        waitbar(idx./prod(size_eng),h,'Testing Sweep in Progress...')
        idx = idx + 1;
        
    end
end

message = sprintf('INST:SEL AMPL'); % select the amplifier channel
fprintf(FSW,message)
message = sprintf('CONF:GEN:RFO:STAT OFF;*WAI');
fprintf(FSW,message)

close(h) %close waitbar
            
%% Safe Bench Equipment
NRP_Close(NRP18S) %Power Meter
NRP_Close(NRP18SN) %Power Meter
E36232A_Close(gate_supply) %Power supply
E36232A_Close(drain_supply) %Power supply
HP437B_Close(HP437B) %Power Meter
FSW_Close(FSW)
% SMJ100A_Close(SMJ100A) %Vector Signal Generator

%% Calculate Pout, Gain, PAE no DPD

FundTone_OutFreq = adjpow(:,1);

Pout = FundTone_Out + output_delta_2_avg(10) + output_attenuation_avg(10);
Pout_watts = 10.^(Pout./10)./1000;

Pin = FundTone_InFreq + input_delta_avg(10);
Pin_watts = 10.^(Pin./10)./1000;

Pdc_Drain1 = abs(Vd1_eng.*Id1_eng);
Pdc_Gate1 = abs(Vg1_eng.*Ig1_eng);
Pdc_tot = Pdc_Drain1 + Pdc_Gate1;

DE = (Pout_watts./Pdc_tot).*100;
PAE = ((Pout_watts-Pin_watts)./Pdc_tot).*100;

Gain = Pout - Pin;

%% Calculate Pout, Gain, PAE with DPD

FundTone_OutFreq_dpd = adjpow_dpd(:,1);

Pout_dpd = FundTone_Out_dpd + output_delta_2_avg(10) + output_attenuation_avg(10);
Pout_watts_dpd = 10.^(Pout_dpd./10)./1000;

Pin_dpd = FundTone_InFreq_dpd + input_delta_avg(10);
Pin_watts_dpd = 10.^(Pin_dpd./10)./1000;

Pdc_Drain1_dpd = abs(Vd1_eng_dpd.*Id1_eng_dpd);
Pdc_Gate1_dpd = abs(Vg1_eng_dpd.*Ig1_eng_dpd);
Pdc_tot_dpd = Pdc_Drain1_dpd + Pdc_Gate1_dpd;

DE_dpd = (Pout_watts_dpd./Pdc_tot_dpd).*100;
PAE_dpd = ((Pout_watts_dpd-Pin_watts_dpd)./Pdc_tot_dpd).*100;

Gain_dpd = Pout_dpd - Pin_dpd;

%% Create Data Summary

data.Pavs = Pavs;
data.Pin = Pin;
data.Pin_dpd = Pin_dpd;
data.Gain = Gain;
data.Gain_dpd = Gain_dpd;
data.Pout = Pout;
data.Pout_dpd = Pout_dpd;
data.PAE = PAE;
data.PAE_dpd = PAE_dpd;
data.Fund_in = FundTone_In;
data.Fund_in_dpd = FundTone_In_dpd;
data.Fund_out = FundTone_Out;
data.Fund_out_dpd = FundTone_Out_dpd;
data.Vgate = Vg1_eng;
data.Vgate_dpd = Vg1_eng_dpd;
data.Igate = Ig1_eng;
data.Igate_dpd = Ig1_eng_dpd;
data.Vdrain = Vd1_eng;
data.Vdrain_dpd = Vd1_eng_dpd;
data.Idrain = Id1_eng;
data.Idrain_dpd = Id1_eng_dpd;
data.EVM = EVM;
data.EVM_dpd = EVM_dpd;
data.CFAC_out = CFAC_out;
data.CFAC_out_dpd = CFAC_out_dpd;
data.AMAM = AMAM;
data.AMAM_dpd = AMAM_dpd;
data.AMPM = AMPM;
data.AMPM_dpd = AMPM_dpd;

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
% 
% close all; clear all; clc
% 
% FSW = FSW_Setup(9,22);
% pause(5) % wait a second
% 
% 
% message = sprintf('FETC:MACC:REVM:CURR:RES?'); % get EVM with no DPD
% EVM = str2num(query(FSW, message)) % get EVM with no DPD
% 
% message = sprintf('FETC:POW:CFAC:OUT:CURR:RES?'); % GRAB THE OUTPUT CREST FACTOR
% CFAC_out = str2num(query(FSW, message))
% 
% message = sprintf('FETC:AMAM:CWID:CURR:RES?'); % GRAB AM/AM
% AMAM = str2num(query(FSW, message))
% 
% message = sprintf('FETC:AMPM:CWID:CURR:RES?'); % GRAB AM/PM
% AMPM = str2num(query(FSW, message))
% 
% pause(5) % wait a second
% FSW_Close(FSW);





