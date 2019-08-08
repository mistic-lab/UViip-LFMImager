close all
clear 
clc
%% Parameters
delta=10e-4; % Pulse width
fs=5e6;       % Sampling Frequenacy ( it has to be less than the attenuator)
noise_fac=1;  % Noise factor 1~ means noise like signal
B=1e6;
t_taw=0:1/fs:delta-1/fs;
u=B/delta;

%% Signal

tx=exp(1j*pi*u*t_taw.^2);
RefSig=tx.';
orgSigSz=size(RefSig);
M=10;R=2;
% x_tx=[zeros(2*length(RefSig),1); RefSig ;zeros(2*length(RefSig),1)];
x_tx=(RefSig);
figure(1)
subplot(211)
plot(real(x_tx));
 title('Ref signal without phase noise')
code=ones(M,1);
TcSigOrg=TC_OLA_tx(x_tx,M,R,1,code);
rng(200)
code_noise=2*pi*rand(1,length(TcSigOrg));
PhaseSignal1=exp(1j*noise_fac*code_noise);
TcSig=TcSigOrg.*PhaseSignal1;
TcSigSz=size(TcSig);
pulsed=TcSig;
pulsed=[pulsed ];%pulsed pulsed pulsed  pulsed];
pulse_sz=size(pulsed) ;
subplot(212)
plot(real(pulsed))
 title('TCsignal')
write_complex_binary(pulsed.','chirpgnu');