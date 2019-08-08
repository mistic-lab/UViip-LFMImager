%%  Header
% Tx combined into one file since I'm having variable sharing issues

%% Initializations
clear all; %normally 'all' is bad but I'm being paranoid about this toolbox
radio_direction = 'Tx';

%% discover_USRP
connectedRadios = findsdru;
if connectedRadios(1).Status == 'Success' %these (1) need attention
    radioFound = true;
    platform = connectedRadios(1).Platform;
    address = connectedRadios(1).IPAddress;
else
    radioFound = false;
    disp('No radio found');
    address = '192.168.10.2';
    platform = 'N200/N210/USRP2';
end

%% SDRU_params
%% Sets MasterClockRate - leaving the switch for when we have to adapt to limeSDR or whatever
switch platform
  case {'B200','B210'}
    RF_params.RadioMasterClockRate = 20e6; %Hz
  case {'X300','X310'}
    RF_params.RadioMasterClockRate = 200e6; %Hz
  case {'N200/N210/USRP2'}
    RF_params.RadioMasterClockRate = 100e6; %Hz
  otherwise
    error(message('sdru:examples:UnsupportedPlatform', platform))
end

%% Sets radio params shared by Tx and Rx
RF_params.CenterFrequency = 3.55e6; %Putting this here so it can be shared by Rx and Tx
RF_params.RadioSampleRate = 200e3; % Hz
RF_params.StopTime = 5; % Simulation stop time in seconds

%% Sets radio params specific to Tx or Rx
RF_params.RadioGain = 0; % apparently via radio.info() this is both min and max
% if radio_direction == 'Rx'
%     RF_params.RadioGain = 0;
% elseif radio_direction == 'Tx'
%     RF_params.RadioGain = 100;
% end
% There is some more commented out stuff here

% RF_params.ChirpSweeptime = 0.5; % seconds I think
% RF_params.ChirpInitialFreq = 3.5e6; % Hz I think
% RF_params.ChirpTargetFreq = 4e6;
% RF_params.SourceFrameTime = 0.06;



%% configure_radio
%% Pull object identifiers from SDRuReceiver
if radio_direction == 'Rx'
    switch platform
        case {'X300','X310'} %again just leaving this here for when we have to adapt to limeSDR or whatever
            radio = comm.SDRuReceiver(...
                'Platform', platform, ...
                'IPAddress', address, ...
                'MasterClockRate', RF_params.RadioMasterClockRate);
        case {'N200/N210/USRP2'}
            radio = comm.SDRuReceiver(...
                'Platform', platform, ...
                'IPAddress', address);
    end
elseif radio_direction == 'Tx'
    switch platform
        case {'X300','X310'} %again just leaving this here for when we have to adapt to limeSDR or whatever
            radio = comm.SDRuTransmitter(...
                'Platform', platform, ...
                'IPAddress', address, ...
                'MasterClockRate', RF_params.RadioMasterClockRate);
        case {'N200/N210/USRP2'}
            radio = comm.SDRuTransmitter(...
                'Platform', platform, ...
                'IPAddress', address);
    end
%     radio.InterpolationFactor =

end

%% Set some experiment specific parameters for the radio
radio.CenterFrequency  = RF_params.CenterFrequency;
radio.Gain = RF_params.RadioGain;
% radio.DecimationFactor = RF_params.RadioDecimationFactor;
% radio.SamplesPerFrame = RF_params.RadioFrameLength;
radio.TransportDataType = 'int16'; % I think this implies 16bits/sample
%radio.OutputDataType = 'Same as transport data type';

%% These should be in Tx file since Rx won't use them
% radio.ChirpSweeptime = RF_params.ChirpSweeptime; % seconds I think
% radio.ChirpSweeptime = RF_params.ChirpInitialFreq; % Hz I think
% radio.ChirpSweeptime = RF_params.ChirpTargetFreq;


%% Display the hwInfo
hwInfo = info(radio)






%% From main

%%Create or load
load = 1;
if load == 1
    %% Get file for transmission
    default_Tx_file = '../Transmission-files/100kHz_1.5s_fs200k_2vec.wav';
    
    %Check whether to use the normal 0.5s sweep
    prompt = ['Use default Tx file? [y/n] (default = ',default_Tx_file,') \n~>'];
    usedefaultTx = input(prompt,'s');
    
    if usedefaultTx == 'y'
        Tx_file = default_Tx_file;
    elseif usedefaultTx == 'n'
        [transmitted_filename, transmitted_path] = uigetfile('*.wav','Select transmitted WAV file');
        Tx_file = [transmitted_path transmitted_filename];
    end
    
    [transmission_sweep, fs] = audioread(Tx_file);
    transmission_sweep = complex(transmission_sweep(:,1),transmission_sweep(:,2));
else
    %% Create sweep
    fs = 200e3;     % sampling frequency
    T = 1.5;        % length of sweep
    t = 0:1/fs:T;
    f1 = -50e3;    % initial frequency
    f2 = 50e3;     % final frequency
    sweep = exp(1j*2*pi*(f1.*t + ((f2-f1)/(2*T)).*t.^2));
    
%     transmission_sweep = transpose(sweep);
    transmission_sweep = transpose([real(sweep); imag(sweep)]);
    audiowrite('../Transmission-files/100kHz_1.5s_fs200k_2vec.wav', transmission_sweep, fs);
    
end

% This is absolutely killing the computer. Crashes school comp one and takes ~1min
% on my macbook
% spectrogram(transmission_sweep, 2024, 2023, 512, fs,'yaxis');

gain = 0.75;
transmission_sweep = transmission_sweep*gain;

%% Transmit file to USRP
if radioFound
    for counter = 1:3
        radio(transmission_sweep)
    end
else
  warning(message('sdru:sysobjdemos:MainLoop'))
end



release(radio)
