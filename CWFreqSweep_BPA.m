%% CW Measurement over Frequency

close all; clear all; clc
%% External Tools
dirPath = pwd;
addpath(strcat(dirPath,'\EquipmentTools'))
addpath(strcat(dirPath,'\RSTools'))
addpath(strcat(dirPath,'\Scripts'))

% load('CalibrationValues.mat')


%% Prepare Sweeps
Fo = 5.6e9:0.2e9:12.4e9; %Center Frequency (Hz)

% THIS RANGE IS SOMETHING YOU NEED TO FIGURE OUT YOURSELF
Pdes = -30:1:-10; %Power set on signal generator (dBm)
Pavs = repmat(Pdes, length(Fo), 1)'; % shaping

%% Prepare Instruments

% N6705A Power Supply Settings
Gate1 = 1; %channel
Drain1 = 2; %channel
Gate2 = 3; %channel
Drain2 = 4; %channel

NRP = NRP_Setup() %Output Power Meter
N6705A = N6705A_Setup(8,5) %DC Power Supply
NRP_control = NRP_50SN_setup() % Control Input Power Meter

SMJ100A = SMJ100A_Setup(8,28) %Vector Signal Generator - control

%% Prepare Engineering Datasets
size_eng = [length(Pavs(:,1)), length(Pctrl_phaseOffset)];

% assuming two devices
Vd1_eng = zeros(size_eng);
Vg1_eng = zeros(size_eng);
Id1_eng = zeros(size_eng);
Ig1_eng = zeros(size_eng);

Vd2_eng = zeros(size_eng);
Vg2_eng = zeros(size_eng);
Id2_eng = zeros(size_eng);
Ig2_eng = zeros(size_eng);

% Three different measured power values at each sweep point
FundTone_In = zeros(size_eng);
FundTone_Out = zeros(size_eng);

FundTone_InFreq = zeros(size_eng);
FundTone_OutFreq = zeros(size_eng);

FundTone_InFreq2 = zeros(size_eng);
FundTone_OutFreq2 = zeros(size_eng);

%% Conduct Sweeps
idx = 0; %temp value for recording progress on waitbar
h = waitbar(0, 'Testing Sweep in Progress...')
SMJ100A_PowerSet(SMJ100A,-60)
SMJ100A_OnOff(SMJ100A, 'on')

for j = 1:1:length(Fo)
    SMJ100A_FreqSet(SMJ100A, Fo(j))
    SMJ100A_PowerSet(SMJ100A,-30)
    pause(5) % wait five seconds
    for i = 1:1:length(Pavs(:,1))
        SMJ100A_PowerSet(SMJ100A, Pavs(i,j))
        pause(5e-1) % wait half a second
        %RECORD RESULTS

        FundTone_In(i,j) = NRP_ReadPower(NRP_control, Fo(j));
        FundTone_Out(i,j) = NRP_ReadPower(NRP, Fo(j));

        FundTone_InFreq(i,j) = NRP_ReadPower(NRP_control, Fo(j));
        FundTone_OutFreq(i,j) = NRP_ReadPower(NRP, Fo(j));

        FundTone_InFreq2(i,j) = NRP_ReadPower(NRP_control, Fo(j));
        FundTone_OutFreq2(i,j) = NRP_ReadPower(NRP, Fo(j));
        
        Vd1_eng(i,j) = N6705A_GetVal(N6705A, Drain1, 'voltage');
        Vg1_eng(i,j) = N6705A_GetVal(N6705A, Gate1, 'voltage');
        Id1_eng(i,j) = N6705A_GetVal(N6705A, Drain1, 'current');
        Ig1_eng(i,j) = N6705A_GetVal(N6705A, Gate1, 'current');
        
        Vd2_eng(i,j) = N6705A_GetVal(N6705A, Drain2, 'voltage');
        Vg2_eng(i,j) = N6705A_GetVal(N6705A, Gate2, 'voltage');
        Id2_eng(i,j) = N6705A_GetVal(N6705A, Drain2, 'current');
        Ig2_eng(i,j) = N6705A_GetVal(N6705A, Gate2, 'current');
            

        % UPDATE WAITBAR
        waitbar(idx./prod(size_eng),h,'Testing Sweep in Progress...')
        idx = idx + 1;
    end %Pavs
    freqChangeFlag = 1;
end %Fo

close(h) %close waitbar
SMJ100A_OnOff(SMJ100A, 'off') % Turn off sig gen
SMJ100A_PowerSet(SMJ100A,-60) % set to a low power value

%% Safe Bench Equipment
N6705A_Close( N6705A ) % DC Supply (this does not turn it off - it just terminates the remote connection)
NRP_Close( NRP ) % Power Meter
NRP_Close(NRP_control)
SMJ100A_Close( SMJ100A ) % Vector Signal Generator

%% Calculate Pout, Gain, PAE etc.

load('CalibrationValues.mat') % Contains Output Attenuation and Input Offsets, or whatever you want

% Scale these as necessary
Attenuation_scaled = repmat(OutputAttenuation, length(Pdes), 1); % Output Attenuation
InputOffsets_scaled = repmat(InputOffsets, length(Pdes), 1); % Input Offsets

Pout = FundTone_OutFreq2 + Attenuation_scaled; % Use whichever measured power you want. Or all of them.
Pout_watts = 10.^(Pout./10)./1000;
Pin = FundTone_InFreq2 + InputOffsets_scaled;
Pin_watts = 10.^(Pin./10)./1000;

Pdc_Drain1 = Vd1_eng.*Id1_eng;
Pdc_Gate1 = Vg1_eng.*Ig1_eng;
Pdc_tot1 = Pdc_Drain1 + Pdc_Gate1;

Pdc_Drain2 = Vd2_eng.*Id2_eng;
Pdc_Gate2 = Vg2_eng.*Ig2_eng;
Pdc_tot2 = Pdc_Drain2 + Pdc_Gate2;

Pdc_tot = Pdc_tot1 + Pdc_tot2;

DE = (Pout_watts./Pdc_tot).*100;
PAE = ((Pout_watts-Pin_watts)./Pdc_tot).*100;

Gain = Pout - Pin;

%% Plotting - change as necessary

% Input Power Driveups
figure(1)
hold on; box on; grid on; grid minor; set(gcf,'color','w')
yyaxis left
plot(Pin,Pout)
plot(Pin,Gain)
ylabel('Pout (dBm) and Gain (dB)')
yyaxis right
plot(Pin,PAE)
ylabel('PAE (%)')
xlabel('Pin (dBm)')
% legend('6 GHz','8 GHz','10 GHz','12 GHz') % or however you're doing it
hold off

% Output Power Driveups
figure(2)
hold on; box on; grid on; grid minor; set(gcf,'color','w')
yyaxis left
plot(Pout,Gain)
ylabel('Gain (dB)')
% ylim([4 13])
yyaxis right
plot(Pout,PAE)
ylabel('PAE (%)')
xlabel('Pout (dBm)')
% legend('6 GHz','8 GHz','10 GHz','12 GHz')
% ylim([0 40])
hold off

% CW
figure(3)
hold on; box on; grid on; grid minor; set(gcf,'color','w')
plot(Fo, Pout)
plot(Fo, Gain)
plot(Fo, PAE)
hold off



