%% Header
% Tx parent file - this file calls all of the scripts and functions
% necessary to send whatever the chosen transmission is to the chosen SDR
% is

%% Initializations
clear all; %normally 'all' is bad but I'm being paranoid about this toolbox
radio_direction = 'Tx';

%% Discover the radio
disp('Looking for available USRP...');
discover_USRP

%% Set radio parameters
disp('Setting USRP parameters ...');
TxParams = SDRU_params(platform, radio_direction)

%% Configure radio object
disp('Creating radio object for transmission ...');
tx_radio = configure_radio(platform, address, TxParams, radio_direction);

%% Create transmission or load transmission file
disp('Creating data to be transmitted ...');
[transmission_sweep, sweep_time] = createtxsweep(TxParams);

%% Create finite transmission time
prompt = 'For how many seconds should the transmission repeat: \n~>';
transmission_time = input(prompt);
number_of_sweeps = round(transmission_time/sweep_time);
 
 
%% Transmit file to USRP
if radioFound
    for counter = 1:number_of_sweeps
        tx_radio(transmission_sweep)
    end
else
  warning(message('sdru:sysobjdemos:MainLoop'))
end


%% Release(radio)
%  I think it may just release the radio so it's available to a different
%  service
release(tx_radio)