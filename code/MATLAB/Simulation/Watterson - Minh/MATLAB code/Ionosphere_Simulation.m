% =======================================================================
% =======================================================================
% Initialize ============================================================
clear
close all
clc
warning off 
% basic inputs ==========================================================

fc=5;         % MHz  Carrier frequency
F=4000;             % sampling rate: fraction of wave length
V=1e-9;            % m/s MS1 speed 
NFFT=64;         % Number of points in FFT
Nsamples=400;    % Number of route samples 
avPower=-20;     % sigma^2  Raverage power
delaystep=1e-7   % delay discretization setep in s
step_f=0.01;     % Freq axis step kHz

% geometry inputs ========================================================

dBS=5000;     
angleBS=180;
BSx=dBS*cosd(angleBS) % location of transmitter (BS) x-coordinate
BSy=dBS*sind(angleBS)  % location of transmitter (BS) y-coordinate

% locations of point scatterers =========================================

fig=figure;
plot(BSx,BSy,'k^'), hold on
plot([0 0],'ko');

% indirect parameters ===================================================

lambdac=300/fc;    % m wavelength
Dx=lambdac/F;      % m sampling spacing 
ts=Dx/V;           % s time sampling interval
fs=1/ts;           % Hz sampling frequency
kc=2*pi/lambdac;   % propagation constant
fm=V/lambdac       % max Doppler shift
cc=3e8;            % speed of light

%========================================================================
% axes
% =======================================================================
timeaxis=ts.*[0:Nsamples-1];
Doppleraxis=([0:NFFT-1]-NFFT/2)*(fs/(NFFT-1));
faxis=[3.5:step_f:4];    % Freq axis in  kHz
% DELAY AXIS DEPENDS ON MAX DELAY, SET LATER 

% ========================================================================
MS0=-V*timeaxis(end)/2;   % initial location of receiver (MS) x-coordinate

MSx=MS0+V.*timeaxis;    % MS route along x-axis
MSy=zeros(Nsamples,1)';  % MS route along x-axis (y=0)
plot(MSx,MSy,'k','LineWidth',5)

MINx=min(min(BSx,MSx))-1000;
MAXx=max(max(BSx,MSx))+1000;
MINy=min(min(BSy,MSy))-1000;
MAXy=max(max(BSy,MSy))+15000;
% axis([MINx MAXx MINy MAXy])
plot([0 0],[MINy MAXy], 'k:')
plot([MINx MAXx],[0 0], 'k:')

%=========================================================================
% SCENARIO EDITOR 
% ========================================================================
% placing point-scatterers in propagation scenario.


DLayerHeight=90;
ELayerHeight=120;
F1LayerHeight=220;
F2LayerHeight=350;
ElectronDensityNormalized=1000;

DElec=1;
EElec=ceil(ElectronDensityNormalized/100*5);      
F1Elec=ceil(ElectronDensityNormalized/100*40);
F2Elec=ElectronDensityNormalized-DElec-EElec-F1Elec;

DThick=0.01;
EThick=2;
F1Thick=2;
F2Thick=4;
SCx=randn(ElectronDensityNormalized,1)*1e3-2e3;

SCy(1:DElec,1)=DLayerHeight+rand(DElec,1)*DThick*10;    %D Layer
SCy(DElec+1:DElec+EElec,1)=ELayerHeight+rand(EElec,1)*EThick*10;    %E Layer
SCy(DElec+EElec+1:DElec+EElec+F1Elec,1)=F1LayerHeight+rand(F1Elec,1)*F1Thick*10;    %F1 Layer
SCy(DElec+EElec+F1Elec:ElectronDensityNormalized,1)=F2LayerHeight+rand(F2Elec+1,1)*F2Thick*10;    %F2 Layer

SCy=SCy*1e3;
NSC=length(SCx);
scatter(SCx,SCy,'k+');
xlabel('Distance (m)')
ylabel('Distance (m)')

% =======================================================================
% calculate distance matrix 
% =======================================================================

distBSSC=sqrt((BSx-SCx).^2+(BSy-SCy).^2);

distBSSCext=repmat(distBSSC,1,Nsamples);

distSCMS=zeros(NSC,Nsamples);
for ii=1:Nsamples
    distSCMS(:,ii)=sqrt((SCx-MSx(ii)).^2+SCy.^2);
end

distBSSCMS=distBSSCext+distSCMS;

% ======================================================================

distBSMS1aux=sqrt((BSx-MSx).^2+(BSy-MSy).^2);   
distBSMS1=min(min(distBSMS1aux));               % Ref distance is min BSMS dist 

% a=(distBSMS1./distBSSC(:)).*(distBSMS1./distSCMS(:,1));
a=(distBSMS1./sqrt(distBSSC(:))).*(distBSMS1./sqrt(distSCMS(:,1)));  % <-----

DeltaPower=avPower-10*log10(sum(a.^2));
deltaa=10.^(DeltaPower/20);             % to achieve reference power
a=deltaa*a;

% =====================================================================
% Define time-varying complex magnitudes of point scatterer contributions 
% amplitudes remain constant while phases change

aa=zeros(NSC,Nsamples);     % create variable 

for k1=1:Nsamples           % scan route points
    for k2=1:NSC            % scan scatterers
        aa(k2,k1)=a(k2)*exp(-j*kc*distBSSCMS(k2,k1));  % time-varying phase
    end
end

% ======================================================================

distBSSCMS1=distBSSCMS-distBSMS1;     % set a new refernece for delays wrt to 
DelaysNormalized=distBSSCMS1/cc;      % arrival of direct ray, here assumed 
                                      % to be totally blocked 
% DelaysNormalized=distBSSCMS/cc; 
                                      
DelaysNormalized=round(DelaysNormalized/delaystep);  % quantify delays delaystep (s)

auxx=size(DelaysNormalized);

auxx2=max(max(DelaysNormalized))+1;    % to include 0 ns delay
ImpulseResponse=zeros(auxx(2),auxx2);     % Create delay profile with step delaystep (s)

for jj=1:auxx(2)           % scan route locations
    for ii=1:auxx(1)       % scan scatterers
        indexx=DelaysNormalized(ii,jj)+1;               
        ImpulseResponse(jj,indexx)=ImpulseResponse(jj,indexx)+aa(ii,jj);     
                                    % put in corresponding delay bin 
                                    % complex amplitude of delta
    end
end

axisdelayprofile=[0:auxx2-1];   % axis in delaystep units
figure, hold
for ii=1:Nsamples
 stem(axisdelayprofile*delaystep*1e6,abs(ImpulseResponse(ii,:)))  
 % delays in us
 % accumulate deltas with time and delay on same plot
end
xlabel('Delay (\mus)')
ylabel('Relative signal level (lin.units)')
title('Absolute value of time varying impulse response, h(\tau;t)')

FreqResp=zeros(Nsamples,length(faxis));

for k1=1:Nsamples               % scan route points
    for k2=1:length(faxis)      % scan frequencies
        for k3=1:NSC
            wl2=0.3/faxis(k2);
            FreqResp(k1,k2)=[FreqResp(k1,k2) + a(k3)*exp(-j*(2*pi/wl2)*distBSSCMS(k3,k1))];
        end
    end
end


figure;mesh(faxis,timeaxis,20*log10(abs(FreqResp)))
ylabel('Time (s)')
xlabel('Frequency (MHz)')
zlabel('Level (dB)')
title('Time-varying frequency response')

figure;plot(faxis,20*log10(abs(FreqResp(1,:))),'k')
xlabel('frequency (MHz)')
ylabel('level (dB)')
title('Frequency response for first route point')

%=======================================================================
% Impulse response through IFFT

ImpResp=zeros(Nsamples,length(faxis));
for k4=1:Nsamples
    ImpResp(k4,:)=ifft(FreqResp(k4,:));
end
taumax=1/(step_f.*1e6);
step_tau=taumax/(length(faxis)-1);
step_tau=taumax/(length(faxis));               %<--------------?????
tauaxis=[0:length(faxis)-1].*step_tau;

figure;mesh(tauaxis,MSx,abs(ImpResp))
xlabel('delay (s)')
ylabel('route point (m)')
zlabel('level (l.u.)')
title('Time-varying impulse response. Magnitude')

figure;plot(tauaxis,abs(ImpResp(1,:)),'k')
xlabel('delay (s)')
zlabel('level (l.u.)')
title('Impulse response for the first route point. Magnitude')


% =======================================================================
% Tapped delay line
% =======================================================================

apdp=abs(ImpulseResponse(1,:)).^2;
APDP=10*log10(apdp)
minAPDP=-60;
maxAPDP=max(APDP);
maxDel=max(axisdelayprofile*delaystep*1e6);
figure, stem2D(axisdelayprofile*delaystep*1e6,APDP,minAPDP) 
axis([-0.5 maxDel+0.5 minAPDP maxAPDP+10])

xlabel('Delay (\mus)')
ylabel('Relative signal level (dB)')
title('Averaged power delay profile')
 
