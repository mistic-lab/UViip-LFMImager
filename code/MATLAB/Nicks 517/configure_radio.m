 function radio = configure_radio(platform, address, RF_params, radio_direction)
%% Header
% This sets up the radio object for both transmission and reception
% Set up radio object to use the found radio

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
end

%% Set some experiment specific parameters for the radio
radio.CenterFrequency  = RF_params.CenterFrequency;
radio.Gain = RF_params.RadioGain;
radio.TransportDataType = 'int16'; % I think this implies 16bits/sample

%% Display the hwInfo
hwInfo = info(radio)