function RF_params = SDRU_params(platform, radio_direction)
%% Sets all of the RF Options for the USRP link to matlab

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
RF_params.CenterFrequency = 3.55e6; % Putting this here so it can be shared by Rx and Tx
RF_params.RadioSampleRate = 200e3;  % Hz
RF_params.StopTime = 5;             % Simulation stop time in seconds


%% Sets radio params specific to Tx or Rx
RF_params.RadioGain = 0; % apparently via radio.info() this is both min and max
% if radio_direction == 'Rx'
%     RF_params.RadioGain = 0;
% elseif radio_direction == 'Tx'
%     RF_params.RadioGain = 100;
% end