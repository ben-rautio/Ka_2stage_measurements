function [  ] = SMJ100A_TwoToneGen( SMJ100A, Fo, DFreq1, DFreq2, OS )
%SMJ100A_TwoTone Generates a Two Tone Signal of frequencies DFreq1 and
%DFreq2

KeepLocalFile    = 0;                       % waveform only temporarily saved
LocalFileName    = 'awgn.wv';               % The local and remote file name
InstrTargetPath = 'D:\';                   % MS Windows based, e.g. SMU, SMJ, SMATE
StartARB         = 1;                       % start in path A


FMax = DFreq1;
if( DFreq2>FMax )
    FMax = DFreq2;
end

% calculate clock rate for ARB
Fsample = FMax * OS;
disp( ['  Clock Rate      : ', num2str(Fsample/1e6,5), ' MHz'] );
if( Fsample > 300e6)
  error ('*** Error: Max. clock rate of 300 MHz exceeded !');
  return;
end

%ggT of the two frequencies
a = int32(DFreq1);
b = int32(DFreq2);
if( a>0 && b>0 )
    while( a ~= b )
        if( a>b )
            a = a-b;
        else
            b = b-a;
        end
    end
end
iggt = a;

disp( ['  ggT             : ', num2str(iggt), ' Hz' ] );

% useful period width
Tggt = 1.0 / double(iggt);
disp( ['  Signal Duration : ', num2str(Tggt,5), ' s'] );

% calculate number of samples
Nsamples = Tggt * Fsample;
disp( ['  Samples         : ', num2str(Nsamples,5)] );
% Make sure to use less than 16 Msamples
if( Nsamples > 16*1024*1024 )
  error ('*** Error: Maximum number of samples exceeded !');
  return;
end

disp (['  ']);

% generate array of sample points
k = 1:1:Nsamples;
% and convert to time scale
Time = (k-1) * 1.0/Fsample;

% signal 1 and 2
Sig1 = cos( 2.0 * pi * DFreq1 * Time );
Sig2 = cos( 2.0 * pi * DFreq2 * Time );

Sig3 = sin( 2.0 * pi * DFreq1 * Time );
Sig4 = sin( 2.0 * pi * DFreq2 * Time );

% sum up all carriers along the time dimension in our
% result matrix
I_data = (Sig1 + Sig2);
Q_data = (Sig3 + Sig4);

% scale to 1.0 for the max envelope 
MaxEnvelope = max( abs(complex(I_data,Q_data)) );
I_data = (I_data / MaxEnvelope) / sqrt(2);
Q_data = (Q_data / MaxEnvelope) / sqrt(2);

% *************************************************************************
% Setup Waveform Structure
% *************************************************************************

IQInfo.I_data         = I_data;         % I signal
IQInfo.Q_data         = Q_data;         % Q signal
IQInfo.comment        = 'Two Tone';     % optional comment
IQInfo.copyright      = 'Rohde&Schwarz';% optional copyright
IQInfo.clock          = Fsample;        % sample rate
IQInfo.no_scaling     = 0;              % scale automatically
IQInfo.path           = InstrTargetPath;% location on instrument     
IQInfo.filename       = LocalFileName;  % local and remore file name

% *************************************************************************
% Plot Data
% *************************************************************************

rs_visualize( Fsample, IQInfo.I_data, IQInfo.Q_data );

% *************************************************************************
% Instrument Setup
% *************************************************************************

% check for R&S device, we also need the *IDN? result later...
disp( 'Checking instrument...' );
[status, InstrIDN] = rs_send_query( SMJ100A, '*IDN?' );
if( ~status ); clear;  return; end
if( isempty( strfind( InstrIDN, 'Rohde&Schwarz' ) ) )
    disp('This is not a Rohde&Schwarz device.');
    clear; return;
end

% reset the instrument
[status, OPCResponse] = rs_send_query( SMJ100A, '*RST; *OPC?' );
if( ~status ); clear;  return; end
[status] = rs_send_command( SMJ100A, '*CLS' );
if( ~status ); clear;  return; end

% generate and send waveform
[status] = rs_generate_wave( SMJ100A, IQInfo, StartARB, KeepLocalFile );
if( ~status ); clear;  return; end

% Apply RF settings if this is not a pure baseband source
if( isempty( strfind( InstrIDN, 'AFQ' ) ) && isempty( strfind( InstrIDN, 'AMU' ) ) )

    disp( 'Setting up RF...' );
        
    % apply some RF settings
    message =  sprintf('FREQ %.6f GHz; *OPC?',Fo/1e9)
    [status, Result] = rs_send_query( SMJ100A, message );
    if( ~status ); clear;  return; end

    [status, Result] = rs_send_query( SMJ100A, 'POW -30.0 dBm; *OPC?' );
    if( ~status ); clear;  return; end

    % switch output ON
    [status, Result] = rs_send_query( SMJ100A, 'OUTP:STAT ON; *OPC?' );
    if( ~status ); clear;  return; end

end

% Read the instruments error queue
status = rs_check_instrument_errors( SMJ100A );
if( ~status ); clear;  return; end

end

