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
%plot all freq vs Pavs

close all
tiledlayout('flow')
nexttile
hold on
for i = 1:1:numel(Fo)
    %all input coupler values, 3 samples
    plot(Pavs(:,i),FundTone_In(:,i),'DisplayName',"IN s1 " + num2str(Fo(i)))
    plot(Pavs(:,i),FundTone_InFreq(:,i),'DisplayName',"IN s2 " + num2str(Fo(i)))
    plot(Pavs(:,i),FundTone_InFreq2(:,i),'DisplayName',"IN s2 " + num2str(Fo(i)))
    %all output coupler values, 3 samples
    plot(Pavs(:,i), FundTone_Out(:,i),'DisplayName',"OUT s1 " + num2str(Fo(i)))
    plot(Pavs(:,i), FundTone_OutFreq(:,i), 'DisplayName',"OUT s2 " + num2str(Fo(i)))
    plot(Pavs(:,i), FundTone_OutFreq2(:,i), 'DisplayName',"OUT s3 " + num2str(Fo(i)))
end
grid on
% legend
xlabel("PAVS")
ylabel("Coupled Power, Output Power")
hold off

% for this cal
%%
%calculate output attenuation for each sample
close all

avg_inp_off_s = load("MIWAVE_input_cal.mat","avg_inp_off");
avg_inp_off = avg_inp_off_s.avg_inp_off;

out_atten1 = FundTone_In + avg_inp_off - FundTone_Out;
out_atten2 = FundTone_InFreq + avg_inp_off - FundTone_OutFreq;
out_atten3 = FundTone_InFreq2 + avg_inp_off - FundTone_OutFreq2;

tiledlayout("flow")
nexttile
hold on
xlabel("Pavs")
ylabel("Output Atten")
plot(Pavs,out_atten1(:,1:numel(Fo)))
plot(Pavs,out_atten2(:,1:numel(Fo)))
plot(Pavs,out_atten3(:,1:numel(Fo)))

grid on

hold off

nexttile
hold on
xlabel("Freq")
ylabel("Output Atten")
%plotting each sample
for i=15:size(Pavs,1)
    plot(Fo,out_atten1(i,1:numel(Fo)));
    plot(Fo,out_atten2(i,1:numel(Fo)));
    plot(Fo,out_atten3(i,1:numel(Fo)));
end
grid on
hold off

%average the output atten for each sample
avg_output_atten = (out_atten1+out_atten2+out_atten3) ./ 3;
%predict output coupler reading based on this for each sample
pred_output1 = FundTone_In + avg_inp_off - avg_output_atten;
pred_output2 = FundTone_InFreq + avg_inp_off - avg_output_atten;
pred_output3 = FundTone_InFreq2 + avg_inp_off - avg_output_atten;
%calculate error
err1 = FundTone_Out - pred_output1;
err2 = FundTone_OutFreq - pred_output2;
err3 = FundTone_OutFreq2 - pred_output3;
%plot error vs pavs and freq
nexttile
hold on
xlabel("Pavs");
ylabel("Error");
plot(Pavs,err1(:,1:numel(Fo)));
plot(Pavs,err2(:,1:numel(Fo)));
plot(Pavs,err3(:,1:numel(Fo)))
grid on
hold off
%vs freq
nexttile
hold on
for i=15:size(Pavs,1)
    plot(Fo,err1(i,1:numel(Fo)));
    plot(Fo,err2(i,1:numel(Fo)));
    plot(Fo,err3(i,1:numel(Fo)));
end
xlabel("Freq")
ylabel("Error")
grid on
hold off
%think this looks acceptable
%next, do the thru measurement to get probe loss
