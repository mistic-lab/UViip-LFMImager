%% Header
% Rx parent file - this file calls all of the scripts and functions
% necessary to receive whatever the chosen spectrum area is to the chosen
% SDR is

%% Initializations
clear all; %normally 'all' is bad but I'm being paranoid about this toolbox
radio_direction = 'Rx';

%% Discover the radio
disp('Looking for available USRP...');
discover_USRP

%% Set radio parameters
disp('Setting USRP parameters ...');
RxParams = SDRU_params(platform, radio_direction)

%% Configure radio object
disp('Creating radio object for transmission ...');
rx_radio = configure_radio(platform, address, RxParams, radio_direction);


%% Receive spectrum from USRP
if radioFound
    timeCounter = 0;
    while timeCounter < RxParams.StopTime
        [received_spectrum, timeCounter] = step(rx_radio);
        timeCounter = timeCounter + 1/RxParams.RadioSampleRate;
        % Check to see the size of received_spectrum
    end
else
    warning(message('sdru:sysobjdemos:MainLoop'))
end

%% Release(rx_radio)
%  I think it may just release the radio so it's available to a different
%  service
release(rx_radio)