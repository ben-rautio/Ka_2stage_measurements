%% CW Measurement over Frequency

close all; clear all; clc
%% External Tools
dirPath = pwd;
addpath(strcat(dirPath,'\EquipmentTools'))
addpath(strcat(dirPath,'\RSTools'))
% addpath(strcat(dirPath,'\Scripts'))

% load('CalibrationValues.mat')


%% Prepare Sweeps
%max of SMW200A is 40GHz!!
Fo = 38e9:0.2e9:40e9; %Center Frequency (Hz)

% THIS RANGE IS SOMETHING YOU NEED TO FIGURE OUT YOURSELF
Pdes = -45:1:-25; %Power set on signal generator (dBm) (for input cal)
Pavs = repmat(Pdes, length(Fo), 1)'; % shaping

%% Using NRX for output
NRP67T = NRP_Setup_67T_USB(); %coupler meter
% power67 = NRP_ReadPower(NRP67T,40e9)
% NRP_Close( NRP67T )

NRX = NRP67T_Setup(); %output meter for output atten, wait ~5s for read after changing power
% powerNRX = NRP67T_ReadPower(NRX,40e9)
% NRP67T_Close( NRX )

%% Using SMW200A
SMJ100A = SMJ100A_Setup(7,28); %Vector Signal Generator - control


%% Prepare Engineering Datasets
% size_eng = [length(Pavs(:,1)), length(Pctrl_phaseOffset)];
size_eng = [length(Pavs(:,1))];
wbar_len = numel(Pavs);

% Three different measured power values at each sweep point
FundTone_In = zeros(size_eng);
FundTone_Out = zeros(size_eng);

FundTone_InFreq = zeros(size_eng);
FundTone_OutFreq = zeros(size_eng);

FundTone_InFreq2 = zeros(size_eng);
FundTone_OutFreq2 = zeros(size_eng);

%% Conduct Sweeps
idx = 0; %temp value for recording progress on waitbar
h = waitbar(0, 'Testing Sweep in Progress...');
SMJ100A_PowerSet(SMJ100A,-60)
SMJ100A_OnOff(SMJ100A, 'on')

%have plots running while sweep running
%vs freq and Pavs
tiledlayout(1,2)
nexttile
hold on
xlabel("PAVS")
ylabel("Coupled, Output")
hold off
nexttile
hold on
xlabel("Freq")
ylabel("Coupled, Output")
hold off
for j = 1:1:length(Fo)
    SMJ100A_FreqSet(SMJ100A, Fo(j))
    SMJ100A_PowerSet(SMJ100A,-50)
    pause(5) % wait five seconds
    for i = 1:1:length(Pavs(:,1))
        SMJ100A_PowerSet(SMJ100A, Pavs(i,j))
        % pause(5e-1) % wait half a second
        %RECORD RESULTS
        pause(5) %wait 5s, saw this settling time by changing SMJ with knob
        %read coupler power
        FundTone_In(i,j) = NRP_ReadPower(NRP67T, Fo(j));
        %read power to output meter
        FundTone_Out(i,j) = NRP67T_ReadPower(NRX,Fo(j));

        FundTone_InFreq(i,j) = NRP_ReadPower(NRP67T, Fo(j));
        FundTone_OutFreq(i,j) = NRP67T_ReadPower(NRX,Fo(j));

        FundTone_InFreq2(i,j) = NRP_ReadPower(NRP67T, Fo(j));
        FundTone_OutFreq2(i,j) = NRP67T_ReadPower(NRX,Fo(j));

        %plot
        nexttile(1)
        hold on
        scatter(Pavs(i,j), FundTone_In(i,j));
        scatter(Pavs(i,j), FundTone_InFreq(i,j));
        scatter(Pavs(i,j), FundTone_InFreq2(i,j));
        scatter(Pavs(i,j), FundTone_Out(i,j));
        scatter(Pavs(i,j), FundTone_OutFreq(i,j));
        scatter(Pavs(i,j), FundTone_OutFreq2(i,j));
        hold off
        nexttile(2)
        hold on
        scatter(Fo(j), FundTone_In(i,j));
        scatter(Fo(j), FundTone_InFreq(i,j));
        scatter(Fo(j), FundTone_InFreq2(i,j));
        scatter(Fo(j), FundTone_Out(i,j));
        scatter(Fo(j), FundTone_OutFreq(i,j));
        scatter(Fo(j), FundTone_OutFreq2(i,j));
        hold off
        
        % UPDATE WAITBAR
        waitbar(idx./wbar_len,h,'Testing Sweep in Progress...')
        idx = idx + 1;
    end %Pavs
    freqChangeFlag = 1;
end %Fo

close(h) %close waitbar

%% turn off everything

SMJ100A_OnOff(SMJ100A, 'off') % Turn off sig gen
SMJ100A_PowerSet(SMJ100A,-60) % set to a low power value

% Safe Bench Equipment
NRP_Close( NRP67T ) % coupler power meter
NRP67T_Close( NRX ) %output power meter
SMJ100A_Close( SMJ100A ) % Vector Signal Generator

%% Plotting - change as necessary

% for this cal
%%
%calculate output attenuation for each sample
close all

avg_inp_off_s = load("MIWAVE_input_cal.mat","avg_inp_off");
avg_inp_off = avg_inp_off_s.avg_inp_off;

avg_output_atten_s = load("outputCal_data1_raw.mat","avg_output_atten");
avg_output_atten = avg_output_atten_s.avg_output_atten;

%calculate GSG loss for each, plot vs Pavs and Fo

Pred_Pin1 = FundTone_In + avg_inp_off;
GSG_loss1 = Pred_Pin1 - avg_output_atten - FundTone_Out;
Pred_Pin2 = FundTone_InFreq + avg_inp_off;
GSG_loss2 = Pred_Pin2 - avg_output_atten - FundTone_OutFreq;
Pred_Pin3 = FundTone_InFreq2 + avg_inp_off;
GSG_loss3 = Pred_Pin3 - avg_output_atten - FundTone_OutFreq2;

%plot GSG loss vs. freq
tiledlayout("flow")
nexttile
hold on
xlabel("Pavs")
ylabel("GSG loss for each sample")
plot(Pavs,GSG_loss1(:,1:numel(Fo)))
plot(Pavs,GSG_loss2(:,1:numel(Fo)))
plot(Pavs,GSG_loss3(:,1:numel(Fo)))
grid on
hold off
%plot vs freq
nexttile
hold on
xlabel("Freq")
ylabel("GSG loss for each sample")
for i=1:size(Pavs,1)
    plot(Fo,GSG_loss1(i,1:numel(Fo)));
    plot(Fo,GSG_loss2(i,1:numel(Fo)));
    plot(Fo,GSG_loss3(i,1:numel(Fo)));
end 
grid on
hold off

%predict Pout using GSG loss
avg_gsg_loss = (GSG_loss3+GSG_loss2+GSG_loss1)./3;
err1 = (FundTone_In + avg_inp_off - avg_gsg_loss - avg_output_atten) - FundTone_Out;
err2 = (FundTone_InFreq + avg_inp_off - avg_gsg_loss - avg_output_atten) - FundTone_OutFreq;
err3 = (FundTone_InFreq2 + avg_inp_off - avg_gsg_loss - avg_output_atten) - FundTone_OutFreq2;
nexttile
hold on
xlabel("Pavs")
ylabel("Error")
plot(Pavs,err1(:,1:numel(Fo)))
plot(Pavs,err2(:,1:numel(Fo)))
plot(Pavs,err3(:,1:numel(Fo)))
grid on

% GSG loss is half since this is a thru measurement
avg_gsg_loss = avg_gsg_loss ./ 2;

