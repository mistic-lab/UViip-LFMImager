function [transmission_sweep, T] = createtxsweep(TxParams)

prompt = 'Transmission options \n 1: Create signal now (there will be prompts) \n 2: Load existing .wav file. \n~>';
loadfilebool = input(prompt);

if loadfilebool == 2
    %% Get file for transmission
    default_Tx_file = '../Transmission-files/100kHz_1.5s_fs200k_2vec.wav';
    
    %Check whether to use the normal 0.5s sweep
    prompt = ['Use default Tx file? [y/n] (default = ',default_Tx_file,') \n~>'];
    defaultTxbool = input(prompt,'s');
    
    if defaultTxbool == 'y'
        Tx_file = default_Tx_file;
    elseif defaultTxbool == 'n'
        [transmitted_filename, transmitted_path] = uigetfile('*.wav','Select transmitted WAV file');
        Tx_file = [transmitted_path transmitted_filename];
    end
    
    [transmission_sweep, fs] = audioread(Tx_file);
    transmission_sweep = complex(transmission_sweep(:,1),transmission_sweep(:,2));
else
    %% Create sweep
    fs = TxParams.RadioSampleRate; % sampling frequency taken from SDRU_params file
    disp(['Sampling frequency is taken from radio settings (fs = ',fs,')']);
    
    prompt = 'Sweep time (s): \n~>';
    T = input(prompt);	% length of sweep  
           
    t = 0:1/fs:T;
    
    prompt = 'Initial frequency (Hz): \n~>';
    f1 = input(prompt);	% initial frequency
    
    prompt = 'Bandwidth of sweep (Hz): \n~>';
    sweepBW = input(prompt);
    
    f2 = f1+sweepBW;     % final frequency
    disp(['Creating ',T,'s sweep from ',f1,'Hz to ',f2,'Hz']);
    
    sweep = exp(1j*2*pi*(f1.*t + ((f2-f1)/(2*T)).*t.^2));
    
    transmission_sweep = transpose(sweep);   
end


% Potentially another way to create the chirp
% tx_radio.ChirpSweeptime = RF_params.ChirpSweeptime; % seconds I think
% tx_radio.ChirpSweeptime = RF_params.ChirpInitialFreq; % Hz I think
% tx_radio.ChirpSweeptime = RF_params.ChirpTargetFreq;

%% Add gain to transmission sweep
%  without this there is noise across the spectrum
gain = 0.75;
prompt = ['Use default gain value? [y/n] (default = ',gain,') \n~>'];
defaultgainbool = input(prompt,'s');
if defaultgainbool == 'n'
    prompt = 'Gain: \n~>';
    gain = input(prompt);
end

transmission_sweep = transmission_sweep*gain;