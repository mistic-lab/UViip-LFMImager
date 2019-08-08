function [sweep, t] = generate_chirp(fs,f1,f2,T,initial_phase_shift,shift_entire_spectrum)
%GENERATE_CHIRP Generates a frequency sweep (outputs a complex column)
%   fs: sample frequency (Hz)
%   f1: start of frequency sweep (Hz)
%   f2: end of frequency sweep (Hz)
%   T: length of sweep (s)

%% Function management
% Check number of inputs.
if nargin > 6
    error('generate_chirp:TooManyInputs', ...
        'requires at most 4 mandatory & 2 optional inputs');
end
% Fill in unset optional arguments.
switch nargin
    case 4
        initial_phase_shift = 1;
        shift_entire_spectrum = 1;
    case 5
        shift_entire_spectrum=1;
end

%% Functional code
% Make time vector
t = 0:1/fs:T;

% Make main sweep
sweeping_spectrum = exp(1j*(pi*((f2-f1)/T).*t.^2));

% Combine main sweep with initial phase and spectrum shift
sweep  = initial_phase_shift.*sweeping_spectrum.*shift_entire_spectrum;

end
