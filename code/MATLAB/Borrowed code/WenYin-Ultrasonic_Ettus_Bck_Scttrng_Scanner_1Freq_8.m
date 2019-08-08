function varargout=Ultrasonic_Ettus_Bck_Scttrng_Scanner_1Freq_2(varargin)
% ULTRASONIC_ETTUS_BCK_SCTTRNG_SCANNER_1FREQ_2 MATLAB code for Ultrasonic_Ettus_Bck_Scttrng_Scanner_1Freq_2.fig
%      ULTRASONIC_ETTUS_BCK_SCTTRNG_SCANNER_1FREQ_2, by itself, creates a new ULTRASONIC_ETTUS_BCK_SCTTRNG_SCANNER_1FREQ_2 or raises the existing
%      singleton*.

%      H = ULTRASONIC_ETTUS_BCK_SCTTRNG_SCANNER_1FREQ_2 returns the handle to a new ULTRASONIC_ETTUS_BCK_SCTTRNG_SCANNER_1FREQ_2 or the handle to
%      the existing singleton*.

%      ULTRASONIC_ETTUS_BCK_SCTTRNG_SCANNER_1FREQ_2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ULTRASONIC_ETTUS_BCK_SCTTRNG_SCANNER_1FREQ_2.M with the given input arguments.

%      ULTRASONIC_ETTUS_BCK_SCTTRNG_SCANNER_1FREQ_2('Property','Value',...) creates a new ULTRASONIC_ETTUS_BCK_SCTTRNG_SCANNER_1FREQ_2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Ultrasonic_Ettus_Bck_Scttrng_Scanner_1Freq_2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Ultrasonic_Ettus_Bck_Scttrng_Scanner_1Freq_2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help Ultrasonic_Ettus_Bck_Scttrng_Scanner_1Freq_2
% Last Modified by GUIDE v2.5 20-Jul-2015 18:11:00
% Begin initialization code - DO NOT EDIT
gui_Singleton=1;
gui_State=struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Ultrasonic_Ettus_Bck_Scttrng_Scanner_1Freq_2_OpeningFcn, ...
    'gui_OutputFcn',  @Ultrasonic_Ettus_Bck_Scttrng_Scanner_1Freq_2_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback=str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}]=gui_mainfcn(gui_State,varargin{:});
else
    gui_mainfcn(gui_State,varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before Ultrasonic_Ettus_Bck_Scttrng_Scanner_1Freq_2 is made visible.
function Ultrasonic_Ettus_Bck_Scttrng_Scanner_1Freq_2_OpeningFcn(hObject,eventdata,handles,varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Ultrasonic_Ettus_Bck_Scttrng_Scanner_1Freq_2 (see VARARGIN)

% clear all;  % Workspace
% close all;  % Figures
% clc;  % Clear Command Window

% ***************************
% *** Setup/Configuration ***
% ***************************

% Data Points/Array Pointers/Step Counters/Samples/Number of Pixels
% num_pnts = num_smpls = num_pixels
%          = num_intrvls+1
%          = num_steps+1

% Translation
% x(i) , i=1:num_x_steps+1
% y(j) , j=1:num_y_steps+1
% x(k) , k=1:num_z_steps+1

% Angular Rotation
% theta(l) , l=1:num_theta_steps+1

% Newport Translation Stages xyz 0-50 mm travel each
% Combine/Couple x-stage and y-stage in tandem to provide 0-100 mm excursion
% Configure xy-Stages in tandem along the x-axis or
% Align the y-axis stge to the x-axis stage
% Let: x1 = x
%      x2 = y
%      zp = z
%      xp = x1 + x2

% x1(i) , i=1:num_x1_steps+1
% x2(i) , i=1:num_x2_steps+1
% xp(i) , i=1:num_x1_steps+num_x2_steps+1
% zp(k) , k=1:num_z_steps+1

% Newport ESP 300 Universal Motion Controller/Driver
% Translation Stage Control:
% The "?" sign Reports current/actual/instantaneos status which is sams as TP or TP?
% MD? Motion Done status
% PA? Position Absolute command initiates an absolute frame motion and reports current position with "?" symbol.
% PR? Position Relative command intitiates a relative frame motion and reports current position with "?" symbol.
% PA Position Absolute Move Command
% PR Position Relative Move Command
% TP Read Actual Instantaneous/Current Position
% SR Set Right travel limit
% MO Motor On
% WS Wait for motion Stop
% OH Set Home Search to High Speed
% VA Set Velocity
% OR Search for Home: Find negative limit signal with 'or4' or home signal with 'or2'

% Angular Rotation Stage
% theta(l) , l=1:num_theta_steps+1

% Stepper Motor Rotation Stage using Arduino Uno plus DF Robot Motor Shield 2A plus Stepper Motor
% theta = 0 to 360 degrees at 400 steps/revolution given Half-Step Increments
% Full-Step Resolution: 200 steps/revolution

% Scan Mode (options)
% 1) xyz,theta axis scanning (Mode 0 - scan_mode=0)
%    Where x1=x, x2=y, zp=z and xp=x1+x2
%    Set x1_min = x_min
%        x2_min = y_min
%    Set x1_max = x_max
%        x2_max = y_max
% 2) z,theta axis scanning (Mode 1 - scan_mode=1)
%    Where x is constant and translation stages along the x-axis are held stationary at a fixed position.
%          zp=z

% Acoustic Ultrasound Sensor Probe (prb) Signal

% Note
% ****
% num_xSteps or num_x1Steps
% num_ySteps or num_x2Steps
% num_zSteps
% num_thetaSteps

% xStepSize or x1StepSize = (xMax-xMin)/num_xSteps
% yStepSize or x2StepSize = (yMax-yMin)/num_ySteps
% zStepSize = (zMax-zMin)/num_zSteps
% thetaStepSize = (thetaMax-thetaMin)/num_thetaSteps

% Step-Size does not have to be a handle.
% Write to edit box
% set(handles.xStepSize_display,'String',num2str(xStepSize));

% Then read from the edit box in a function
% x1StepSize=str2double(get(handles.xStepSize_display,'String'));

% First Experiment:
% Scan Mode: z,theta axis (or zp,theta) axis Scanning
% zMin=0      zMax=32 mm         num_zSteps=64
% thetaMin=0  thetaMax=57.6 deg  num_thetaSteps=64
% zStepSize=delta_z=abs(zMax-zMin)/num_zSteps=0.5 mm/step
% thetaStepSize=delta_theta=abs(thetaMax-thetaMin)/num_thetaSteps=0.90 deg/step
% ***********
% *** end ***
% ***********

% Setup Serial Communication for Newport ESP 300 Translation Stage Controller
out=instrfind('Port','COM4');  % Check COM4 Port Status
if ~isempty(out)  % Check for valid serial connection on Port - COM4
    fclose(out);  % Close COM4 serial connection if serial object exists.
end

% Newport Universal Motion Controller/Driver ESP 300 - Translation Stage Controller
handles.s=serial('COM4','BaudRate',19200,'FlowControl','hardware','Terminator','CR');
fopen(handles.s);  % Open Serial Port then Home All Translation Stages
fprintf(handles.s,'1mo;2mo;3mo;1oh0.5;2oh0.5;3oh0.5;1or4;1ws;2or4;2ws;3or4;3ws');  % Turn On All Axes and Zero to Negative Hardware Limit
fprintf(handles.s,'1VA2;2VA2;3VA3');  % Set Axes Speed/Velocity to 2 mm/sec and 0.5 mm/sec
% fprintf(handles.s,'3AC2');  % Set Axes acceleration 2 mm/sec^2
fprintf(handles.s,'1SR50;2SR50;3SR50');  % Set Travel Limit to 50 mm Excursion

% Setup Serial Communication for Rotation Stage (Stepper Motor) Controller
% Arduino Rotational Stage Motion Stepper Control Simulation Driver Program
% find Port /dev/tty.usbmodem1421 Status
% out=instrfind('Port','/dev/tty.usbmodem1411');  % Find the Serial Port
out=instrfind('Port','COM3');  % Find the Serial Port
if isa(out,'serial') && isvalid(out)
    % Check whether the Serial Port is currently used by Matlab
    if ~isempty(out)
        % Port Activated and Reserved
        disp('WARNING: Port COM3 is open and in use. Now closing.');
        % Check if the Serial Port status is open
        if strcmpi(get(out,'status'),'open')
            fclose(out);  % If open then close
        end
        delete(out);
        clear out;  % clear('out')
    end
else
    disp('Looks like board is already connected to port COM3');
    disp('Delete comPort object to force disconnection before attempting to connect to a diffrent port.');
    delete(out);
    clear out;
end

% Create Serial Port communication object for Arduino USB interface
% handles.ser=serial('/dev/tty.usbmodem1411');
handles.ser=serial('COM3');
% Set serial port commnunication object property value
set(handles.ser, 'BytesAvailableFcnMode','Terminator');
% set(handles.ser, 'BytesAvailableFcn',@stepNum_update,handles);  original
% set(handles.ser, 'BytesAvailableFcn',{@stepNum_update,handles});  % Note: set(handles.ser, 'BytesAvailableFcn',@stepNum_update,handles) does not work!
set(handles.ser,'Terminator','LF')  % handles.ser.Terminator=LF (Line Feed or New Line char '\n' or ASCII Dec 10)
set(handles.ser,'DataBits',8);   % handles.ser.Databits=8 bits
set(handles.ser,'StopBits',1);   % handles.ser.Stopbits=1 bits
set(handles.ser,'BaudRate',9600);  % handles.ser.BaudRate=9600 bps
set(handles.ser,'Parity','none');  % handles.ser.parity=none
set(handles.ser,'TimerPeriod',4);  % handles.ser.TimerPeriod=4 sec
set(handles.ser,'InputBufferSize',256);  % handles.ser.InputBufferSize=256 bytes
set(handles.ser,'OutputBufferSize',256);  % handles.ser.OutputBufferSize=256 bytes
fopen(handles.ser);  % It takes several seconds before any operation could be attempted
disp(handles.ser);

% Get Translation Stage Minimum Position Values
x1Min=str2double(get(handles.xMin_display,'String'));
x2Min=str2double(get(handles.yMin_display,'String'));
zpMin=str2double(get(handles.zMin_display,'String'));
% thetaMin=str2double(get(handles.thetaMin_display,'String'));

% Move All Translation Stages to Minimum Position
x1_pos=move_abs_trans_stage(hObject,eventdata,handles,1,x1Min);  % Commanded Position: x1Min
x2_pos=move_abs_trans_stage(hObject,eventdata,handles,2,x2Min);  % Commanded Position: x2Min
zp_pos=move_abs_trans_stage(hObject,eventdata,handles,3,zpMin);  % Commanded Position: zpMin

% Rotation Stage Stepper Motor is preset to start position when configured in experiment.

% Display Home (Current) Position
set(handles.x_pos_display,'String',num2str(x1_pos));
set(handles.y_pos_display,'String',num2str(x2_pos));
set(handles.z_pos_display,'String',num2str(zp_pos));
% set(handles.theta_pos_display,'String',num2str(theta_pos));

% Set Array Size to Low Value (2) as GUI can update value by
% num_x1Steps=str2double(get(handles.num_xSteps_display,'String'))
% num_x2Steps=str2double(get(handles.num_ySteps_display,'String'))
% num_zpSteps=str2double(get(handles.num_zSteps_display,'String'))
% num_thetaSteps=str2double(get(handles.num_thetaSteps_display,'String'))
num_x1Steps=1;
num_x2Steps=1;
num_zpSteps=1;
num_thetaSteps=1;

% Initialize Index Pointers to Arrays for Step Counters
handles.indx_xp=1;
handles.indx_zp=1;
handles.indx_theta=1;

% Array Dimension Size and Initialization
% Scanning Variables
% x1(1:num_x1Steps+1,1:num_zSteps+1,1:num_thetaSteps+1)=0;
% x2(1:num_x2Steps+1,1:num_zSteps+1,1:num_thetaSteps+1)=0;
xp(1:num_x1Steps+num_x2Steps+1,1:num_zpSteps+1,1:num_thetaSteps+1)=0;
zp(1:num_x1Steps+num_x2Steps+1,1:num_zpSteps+1,1:num_thetaSteps+1)=0;
th(1:num_x1Steps+num_x2Steps+1,1:num_zpSteps+1,1:num_thetaSteps+1)=0;

% Output Data Variables for Sensor Probe (prb)  Note: No Reference (ref)
mag_ft(1:num_x1Steps+num_x2Steps+1,1:num_zpSteps+1,1:num_thetaSteps+1)=0;  % Fourier Transform (FT) Method
phz_ft(1:num_x1Steps+num_x2Steps+1,1:num_zpSteps+1,1:num_thetaSteps+1)=0;
mag_iq(1:num_x1Steps+num_x2Steps+1,1:num_zpSteps+1,1:num_thetaSteps+1)=0;  % In-Phase Quadrature (IQ) Method
phz_iq(1:num_x1Steps+num_x2Steps+1,1:num_zpSteps+1,1:num_thetaSteps+1)=0;

% Select Scan Mode:
% 1) xyz,theta axis (or xp,zp,theta axis) Scanning
% 2) z,theta axis (or zp,theta axis) Scanning
% handles.scan_mode=0  for xyz,theta axis Scanning
% handles.scan_mode=1  for z,theta axis Scanning
handles.scan_mode=1;

% Select Data Processing Type
% 1) Fourier Transform (FT) Method
% 2) In-Phase Quadrature (IQ) Method
% handles.data_processing_method=0  for Fourier Transform (FT) Data Processing
% handles.data_processing_method=1  for In-Phase Quadrature (IQ) Data Processing
handles.data_processing_method=0;  % Fourier Transform Method

% Scrolling Data Window - Strip Chart Plot
handles.numWinPnts=26;  % Window Size

% Time: time(1:handles.numWinPnts)
delta_time=200/1000;  % Arbitrary Time Step Interval or Sample Time: delta_time = smpl_tme (sec/sample)
time_min=0;  % Initialization
time_max=time_min+(handles.numWinPnts-1)*delta_time;
time=time_min:delta_time:time_max;  % Time: time(1:handles.numWinPnts)
% time(1:handles.numWinPnts)=0;  % time=zeros(1,handles.numWinPnts)
% time_sum=time_min;  % Initialization Start Condition: time(1)=time_min
% for i=1:handles.numWinPnts
%     time(i)=time_sum;
%     time_sum=time_sum+delta_time;  % Update
% end
% time(handles.numWinPnts) = time_max = tme_sum-delta_tme = 0+(26-1)*0.20 = 25*0.20 =  5 sec
% time=time_min:delta_time:time_max;  % time=time_max*linspace(0,1,handles.numWinPnts);
handles.time_plt=time;  % Time Array Handle: handles.time_plt(1:handles.numWinPnts)

% Strip-Chart Scrolling Data Plot Arrays
mag_ft_plt(1:handles.numWinPnts)=0;  % Fourier Transform (FT) Method
phz_ft_plt(1:handles.numWinPnts)=0;
mag_iq_plt(1:handles.numWinPnts)=0;  % In-Phase Quadrature (IQ) Method
phz_iq_plt(1:handles.numWinPnts)=0;

% Time and Frequency Arrays
% Why 128?
handles.tme(1:128)=0;  % handles.tme=zeros(1,128);
handles.freq(1:128)=0;  % handles.freq=zeros(1,128);

% Index Pointers at Peak in Power Spectrum
% default index of f_max=206?
handles.indx_at_pwr_spec_pk_prb=16;  % Probe (prb)

% Assign Variables/Parameters to Workspace - Initial Dimension Size
% Input Scanning Variables
% assignin('base','X or X1',x1);
% assignin('base','Y or X2',x2);
% assignin('base','Z',zp);
assignin('base','Xp',xp);
assignin('base','Zp',zp);
assignin('base','Theta',th);

% Output Variables
assignin('base','Magn_FT',mag_ft);
assignin('base','Phze_FT',phz_ft);
assignin('base','Magn_IQ',mag_iq);
assignin('base','Phze_IQ',phz_iq);

% Strip-Chart Scrolling Data Plot
assignin('base','Mag_FT_Plt',mag_ft_plt);
assignin('base','Phz_FT_Plt',phz_ft_plt);
assignin('base','Mag_IQ_Plt',mag_iq_plt);
assignin('base','Phz_IQ_Plt',phz_iq_plt);

pause on;  % Enable pause command but does not initiate a pause or delay

% Setup & Configure Ettus N210 USRP Data Acquisition Box
% Quadrature Output at 100 MSps = 100e6 samples/sec using LFTX and LFRX Daughterboards
% Center RF Frequency: Transmitter (Tx) and Receiver (Rx)
paramsTx.CenterFrequency=str2double(get(handles.freq_disp,'String'))*1e6;  % (Hz) (Typical: 1.5 MHz Center Frequency)
paramsRx.CenterFrequency=str2double(get(handles.freq_disp,'String'))*1e6;  % (Hz) (Typical: 1.5 MHz Center Frequency)

% Transmitter (Tx) Parameters
paramsTx.USRPClockRate            = 100e6;  % (samples/sec) Sample Rate or Clock Rate
paramsTx.USRPInterpolationFactor  = 400    % (steps or intervals per tick) Interpolation Scale-Factor
paramsTx.USRPGain                 = 0;      % (dB) Gain x1
paramsTx.USRPFrameLength          = 1000;   % (samples or data point) Number of Samples to Acquire (record length)
paramsTx.StopTime                 = 3;      % (sec) Stop Time Interval for Transmitter to Receiver Streaming
paramsTx.FrontEndSampleRate       = paramsTx.USRPClockRate ...
    /paramsTx.USRPInterpolationFactor;  % (samples/sec) Interpolation Factor Sample Rate

% Receiver (Rx) Parameters
paramsRx.USRPClockRate        = 100e6;  % (samples/sec) Sample Rate or Clock Rate
paramsRx.USRPDecimationFactor = 400;    % (steps or intervals per tick) Decimation Scale-Factor
paramsRx.USRPGain             = 0;      % (dB) Gain x1
paramsRx.USRPFrameLength      = 1000;   % (samples or data point) Number of Samples to Acquire (record length)
paramsRx.FrontEndSampleRate   = paramsRx.USRPClockRate...
    /paramsRx.USRPDecimationFactor;  % (samples/sec) Decimation Factor Sample Rate
paramsRx.USRPFrameTime        = paramsRx.USRPFrameLength ...
    /paramsRx.FrontEndSampleRate;  % (sec) Frame Time

handles.hSDRuTx=comm.SDRuTransmitter(...
    'IPAddress',              '192.168.10.2', ...
    'CenterFrequency',        paramsTx.CenterFrequency,...
    'Gain',                   paramsTx.USRPGain, ...
    'InterpolationFactor',    paramsTx.USRPInterpolationFactor);

handles.hSDRuRx=comm.SDRuReceiver(...
    'IPAddress',             '192.168.10.2', ...
    'CenterFrequency',       paramsRx.CenterFrequency,...
    'Gain',                  paramsRx.USRPGain, ...
    'DecimationFactor',      paramsRx.USRPDecimationFactor, ...
    'SampleRate',            paramsRx.FrontEndSampleRate, ...
    'FrameLength',           paramsRx.USRPFrameLength, ...
    'OutputDataType',        'single');

radio=findsdru(handles.hSDRuRx.IPAddress);  % Find USRP/SDRu Radio
hwInfo = info(handles.hSDRuRx);
disp(hwInfo);
% USRP/SDRu Radio Handles
handles.num_data_smpls=paramsRx.USRPFrameLength;  % handles.num_data_smpls=paramsTx.USRPFrameLength
handles.paramsTxStopTime=paramsTx.StopTime;
handles.paramsRxUSRPFrameTime=paramsRx.USRPFrameTime;
handles.smpl_rate=paramsRx.FrontEndSampleRate;
handles.radio=radio;

% USRP/SDRu Radio Communication to IP Address Status
if (strcmp(radio.Status,'Success'))
    disp('USRP/SDRu IP Address Connection Status: Successful');
else
    warning(message('USRP/SDRu IP Address Connection Status: Not Successful'));
end

% Release all SDRu Resources
% release(handles.hSDRuTx); clear handles.hSDRuTx;
% release(handles.hSDRuRx); clear handles.hSDRuRx;

% Other
handles.messiah=1;
fprintf(['Handels Messiah is #' num2str(handles.messiah,'%2.0f') '\n']);
fprintf('Building a GUI Mystery!\n');

% Number of Data Samples per Trigger or Acquisition Frame (samples or data points)
num_fft_smpls=2^(nextpow2(handles.num_data_smpls));  % Number of FFT Samples: num_fft_smpls = 2^n = 64, 128, 256, 512, 1024 etc.
% Number of Nyquist Samples at Nyquist Limit: num_nyq_smpls
handles.num_nyq_smpls=num_fft_smpls/2+1;  % Index Pointer to Nyquist Limit
% freq_nyq=smpl_rate/2;  % Nyquist Frequency: nyq_freq (cycles/sec or Hz)
% Note: Nyquist Frequency at 2 samples/cycle
% freq(1)=0;  % DC component
% delta_freq = (freq_nyq-freq(1))/(num_nyq_smpls-1)
%            = (smpl_rate/2-0)/((num_fft_smpls/2+1)-1)
%            = (smpl_rate/2)/(num_fft_smpls/2)
%            = smpl_rate/num_fft_smpls
% delta_freq=smpl_rate/num_fft_smpls;
% freq=freq(1):delta_freq:freq_nyq;  % freq=freq_nyq*linspace(0,1,num_nyq_smpls);
% freq_nyq = freq(num_nyq_smpls)
fprintf('Number of Samples per Trigger or Frame %6.0f\n',handles.num_data_smpls);
fprintf('Front-End Sample Rate (kHz or x1000 samples/sec) %6.1f\n',handles.smpl_rate/1000);
fprintf('FFT point number %6.0f\n',num_fft_smpls);

% Data Window
num_data_pnts=handles.num_data_smpls;  % num_smpls=num_pnts
num_win_pnts=50;
num_zero_pnts=10;
data_win(1:num_data_pnts)=1;  % data_win=ones(1,num_data_pnts)
data_win(1:num_zero_pnts+1)=0;  % Zero Asymptote LHS
data_win(num_data_pnts-num_zero_pnts:num_data_pnts)=0;  % Zero Asymptote RHS
slope=1/(num_win_pnts-num_zero_pnts);
dx=1;
dy=slope*dx;
y_sum=dy;  % Initialization
for i=1:num_win_pnts-num_zero_pnts-1
    data_win(i+num_zero_pnts+1)=y_sum;
    data_win((num_data_pnts-num_zero_pnts)-i)=y_sum;
    y_sum=y_sum+dy;  % Update
end
handles.win=data_win;  % Data Window Handle

% *** Diagnostics ***
% save('debug_info_1','handles');  % For debugging
% *** end ***

% Close-Out Commands
handles.output=hObject;  % Choose default command line output for Ultrasonic_Multi_Stage_Scanner_1
guidata(hObject,handles);  % Update GUI handles structure

% function stepNum_update(obj,event,handles)
% % Pointer (Ptr) to Stepper Motor Step Control Function
% % Maintains the pointer or index value in the modulus function for multiple stepping calls.
% % global Ptr;  % Modulus: Ptr = stepCount % 8  in Arduino Code
% % Ptr=0;  % Initial Condition
% % handles.rot_stage_speed=2.0;  % Rotation Stage Speed: rot_stage_speed (deg/sec)
% arduino=fscanf(handles.ser,'%s');
% global Num;
% [step,Num]=strread(arduino,'%s%u',1,'delimiter',',');
% assignin('base','stepNumber',Num);
% set(handles.readData_indicator,'String',Num,'Value',Num);
% drawnow;
% end

% UIWAIT makes Ultrasonic_Ettus_Bck_Scttrng_Scanner_1Freq_2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% --- Outputs from this function are returned to the command line.
function varargout=Ultrasonic_Ettus_Bck_Scttrng_Scanner_1Freq_2_OutputFcn(hObject,eventdata,handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1}=handles.output;

function xMin_display_Callback(hObject,eventdata,handles)
% hObject    handle to xMin_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of xMin_display as text
%        str2double(get(hObject,'String')) returns contents of xMin_display as a double
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function xMin_display_CreateFcn(hObject,eventdata,handles)
% hObject    handle to xMin_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function yMin_display_Callback(hObject,eventdata,handles)
% hObject    handle to yMin_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of yMin_display as text
%        str2double(get(hObject,'String')) returns contents of yMin_display as a double
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function yMin_display_CreateFcn(hObject,eventdata,handles)
% hObject    handle to yMin_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function zMin_display_Callback(hObject,eventdata,handles)
% hObject    handle to zMin_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of zMin_display as text
%        str2double(get(hObject,'String')) returns contents of zMin_display as a double
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function zMin_display_CreateFcn(hObject,eventdata,handles)
% hObject    handle to zMin_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function thetaMin_display_Callback(hObject,eventdata,handles)
% hObject    handle to thetaMin_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of thetaMin_display as text
%        str2double(get(hObject,'String')) returns contents of thetaMin_display as a double
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function thetaMin_display_CreateFcn(hObject,eventdata,handles)
% hObject    handle to thetaMin_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function xMax_display_Callback(hObject,eventdata,handles)
% hObject    handle to xMax_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of xMax_display as text
%        str2double(get(hObject,'String')) returns contents of xMax_display as a double
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function xMax_display_CreateFcn(hObject,eventdata,handles)
% hObject    handle to xMax_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function yMax_display_Callback(hObject,eventdata,handles)
% hObject    handle to yMax_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of yMax_display as text
%        str2double(get(hObject,'String')) returns contents of yMax_display as a double
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function yMax_display_CreateFcn(hObject,eventdata,handles)
% hObject    handle to yMax_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function zMax_display_Callback(hObject,eventdata,handles)
% hObject    handle to zMax_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of zMax_display as text
%        str2double(get(hObject,'String')) returns contents of zMax_display as a double
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function zMax_display_CreateFcn(hObject,eventdata,handles)
% hObject    handle to zMax_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function thetaMax_display_Callback(hObject,eventdata,handles)
% hObject    handle to thetaMax_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of thetaMax_display as text
%        str2double(get(hObject,'String')) returns contents of thetaMax_display as a double
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function thetaMax_display_CreateFcn(hObject,eventdata,handles)
% hObject    handle to thetaMax_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function num_xSteps_display_Callback(hObject,eventdata,handles)
% hObject    handle to num_xSteps_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of num_xSteps_display as text
%        str2double(get(hObject,'String')) returns contents of num_xSteps_display as a double
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function num_xSteps_display_CreateFcn(hObject,eventdata,handles)
% hObject    handle to num_xSteps_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function num_ySteps_display_Callback(hObject,eventdata,handles)
% hObject    handle to num_ySteps_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of num_ySteps_display as text
%        str2double(get(hObject,'String')) returns contents of num_ySteps_display as a double
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function num_ySteps_display_CreateFcn(hObject,eventdata,handles)
% hObject    handle to num_ySteps_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function num_zSteps_display_Callback(hObject,eventdata,handles)
% hObject    handle to num_zSteps_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of num_zSteps_display as text
%        str2double(get(hObject,'String')) returns contents of num_zSteps_display as a double
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function num_zSteps_display_CreateFcn(hObject,eventdata,handles)
% hObject    handle to num_zSteps_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function num_thetaSteps_display_Callback(hObject,eventdata,handles)
% hObject    handle to num_thetaSteps_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of num_thetaSteps_display as text
%        str2double(get(hObject,'String')) returns contents of num_thetaSteps_display as a double
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function num_thetaSteps_display_CreateFcn(hObject,eventdata,handles)
% hObject    handle to num_thetaSteps_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in xSpeed_popup.
function xSpeed_popup_Callback(hObject,eventdata,handles)
% hObject    handle to xSpeed_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns xSpeed_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from xSpeed_popup
val=get(hObject,'Value');
str=get(hObject,'String');
switch str{val}
    case '2.0 mm/sec'
        fprintf(handles.s,'1VA2');
    case '1.0 mm/sec'
        fprintf(handles.s,'1VA1');
    case '0.5 mm/sec'
        fprintf(handles.s,'1VA0.5');
    case '0.2 mm/sec'
        fprintf(handles.s,'1VA0.2');
    case '0.1 mm/sec'
        fprintf(handles.s,'1VA0.1');
    case '0.05 mm/sec'
        fprintf(handles.s,'1VA0.05');
end
guidata(hObject,handles);  % Update GUI handles

% --- Executes during object creation, after setting all properties.
function xSpeed_popup_CreateFcn(hObject,eventdata,handles)
% hObject    handle to xSpeed_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in ySpeed_popup.
function ySpeed_popup_Callback(hObject,eventdata,handles)
% hObject    handle to ySpeed_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns ySpeed_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ySpeed_popup
val=get(hObject,'Value');
str=get(hObject,'String');
switch str{val}
    case '2.0 mm/sec'
        fprintf(handles.s,'2VA2');
    case '1.0 mm/sec'
        fprintf(handles.s,'2VA1');
    case '0.5 mm/sec'
        fprintf(handles.s,'2VA0.5');
    case '0.2 mm/sec'
        fprintf(handles.s,'2VA0.2');
    case '0.1 mm/sec'
        fprintf(handles.s,'2VA0.1');
    case '0.05 mm/sec'
        fprintf(handles.s,'2VA0.05');
end
guidata(hObject,handles);  % Update GUI handles

% --- Executes during object creation, after setting all properties.
function ySpeed_popup_CreateFcn(hObject,eventdata,handles)
% hObject    handle to ySpeed_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in zSpeed_popup.
function zSpeed_popup_Callback(hObject,eventdata,handles)
% hObject    handle to zSpeed_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns zSpeed_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from zSpeed_popup
val=get(hObject,'Value');
str=get(hObject,'String');
switch str{val}
    case '3.0 mm/sec'
        fprintf(handles.s,'3VA3');
    case '2.0 mm/sec'
        fprintf(handles.s,'3VA2');
    case '1.0 mm/sec'
        fprintf(handles.s,'3VA1');
    case '0.5 mm/sec'
        fprintf(handles.s,'3VA0.5');
    case '0.2 mm/sec'
        fprintf(handles.s,'3VA0.2');
    case '0.1 mm/sec'
        fprintf(handles.s,'3VA0.1');
    case '0.05 mm/sec'
        fprintf(handles.s,'3VA0.05');
end
guidata(hObject,handles);  % Update GUI handles

% --- Executes during object creation, after setting all properties.
function zSpeed_popup_CreateFcn(hObject,eventdata,handles)
% hObject    handle to zSpeed_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in thetaSpeed_popup.
function thetaSpeed_popup_Callback(hObject,eventdata,handles)
% hObject    handle to thetaSpeed_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns thetaSpeed_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from thetaSpeed_popup
val=get(hObject,'Value');
str=get(hObject,'String');
switch str{val}
    case '1.0 deg/sec'
        handles.rot_stage_speed=1;
    case '2.0 deg/sec'
        handles.rot_stage_speed=2;
    case '4.0 deg/sec'
        handles.rot_stage_speed=4;
    case '6.0 deg/sec'
        handles.rot_stage_speed=6;
    case '10.0 deg/sec'
        handles.rot_stage_speed=10;
end
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function thetaSpeed_popup_CreateFcn(hObject,eventdata,handles)
% hObject    handle to thetaSpeed_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function xStepSize_display_Callback(hObject,eventdata,handles)
% hObject    handle to xStepSize_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of xStepSize_display as text
%        str2double(get(hObject,'String')) returns contents of xStepSize_display as a double
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function xStepSize_display_CreateFcn(hObject,eventdata,handles)
% hObject    handle to xStepSize_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function yStepSize_display_Callback(hObject,eventdata,handles)
% hObject    handle to yStepSize_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of yStepSize_display as text
%        str2double(get(hObject,'String')) returns contents of yStepSize_display as a double
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function yStepSize_display_CreateFcn(hObject,eventdata,handles)
% hObject    handle to yStepSize_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function zStepSize_display_Callback(hObject,eventdata,handles)
% hObject    handle to zStepSize_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of zStepSize_display as text
%        str2double(get(hObject,'String')) returns contents of zStepSize_display as a double
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function zStepSize_display_CreateFcn(hObject,eventdata,handles)
% hObject    handle to zStepSize_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function thetaStepSize_display_Callback(hObject,eventdata,handles)
% hObject    handle to thetaStepSize_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of thetaStepSize_display as text
%        str2double(get(hObject,'String')) returns contents of thetaStepSize_display as a double
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function thetaStepSize_display_CreateFcn(hObject,eventdata,handles)
% hObject    handle to thetaStepSize_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function x_mov_abs_display_Callback(hObject,eventdata,handles)
% hObject    handle to x_mov_abs_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of x_mov_abs_display as text
%        str2double(get(hObject,'String')) returns contents of x_mov_abs_display as a double
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function x_mov_abs_display_CreateFcn(hObject,eventdata,handles)
% hObject    handle to x_mov_abs_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function y_mov_abs_display_Callback(hObject,eventdata,handles)
% hObject    handle to y_mov_abs_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of y_mov_abs_display as text
%        str2double(get(hObject,'String')) returns contents of y_mov_abs_display as a double
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function y_mov_abs_display_CreateFcn(hObject,eventdata,handles)
% hObject    handle to y_mov_abs_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function z_mov_abs_display_Callback(hObject,eventdata,handles)
% hObject    handle to z_mov_abs_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of z_mov_abs_display as text
%        str2double(get(hObject,'String')) returns contents of z_mov_abs_display as a double
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function z_mov_abs_display_CreateFcn(hObject,eventdata,handles)
% hObject    handle to z_mov_abs_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function theta_mov_abs_display_Callback(hObject,eventdata,handles)
% hObject    handle to theta_mov_abs_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of theta_mov_abs_display as text
%        str2double(get(hObject,'String')) returns contents of theta_mov_abs_display as a double
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function theta_mov_abs_display_CreateFcn(hObject,eventdata,handles)
% hObject    handle to theta_mov_abs_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function x_pos_display_Callback(hObject,eventdata,handles)
% hObject    handle to x_pos_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of x_pos_display as text
%        str2double(get(hObject,'String')) returns contents of x_pos_display as a double
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function x_pos_display_CreateFcn(hObject,eventdata,handles)
% hObject    handle to x_pos_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function y_pos_display_Callback(hObject,eventdata,handles)
% hObject    handle to y_pos_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of y_pos_display as text
%        str2double(get(hObject,'String')) returns contents of y_pos_display as a double
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function y_pos_display_CreateFcn(hObject,eventdata,handles)
% hObject    handle to y_pos_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function z_pos_display_Callback(hObject,eventdata,handles)
% hObject    handle to z_pos_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of z_pos_display as text
%        str2double(get(hObject,'String')) returns contents of z_pos_display as a double
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function z_pos_display_CreateFcn(hObject,eventdata,handles)
% hObject    handle to z_pos_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function theta_pos_display_Callback(hObject,eventdata,handles)
% hObject    handle to theta_pos_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of theta_pos_display as text
%        str2double(get(hObject,'String')) returns contents of theta_pos_display as a double
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function theta_pos_display_CreateFcn(hObject,eventdata,handles)
% hObject    handle to theta_pos_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in xGo_pushbutton.
function xGo_pushbutton_Callback(hObject,eventdata,handles)
% hObject    handle to xGo_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Command Stage Motion
% Make certain to enter a move to ... position in "x_mov_abs_display"
x_mov_abs_cmd=str2double(get(handles.x_mov_abs_display,'String'));
x_stage_pos=move_abs_trans_stage(hObject,eventdata,handles,1,x_mov_abs_cmd);
set(handles.x_pos_display,'String',num2str(x_stage_pos));
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes on button press in yGo_pushbutton.
function yGo_pushbutton_Callback(hObject,eventdata,handles)
% hObject    handle to yGo_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Command Stage Movement
% Make certain to enter a move to ... position in "y_mov_abs_display"
y_mov_abs_cmd=str2double(get(handles.y_mov_abs_display,'String'));
y_stage_pos=move_abs_trans_stage(hObject,eventdata,handles,2,y_mov_abs_cmd);
set(handles.y_pos_display,'String',num2str(y_stage_pos));
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes on button press in zGo_pushbutton.
function zGo_pushbutton_Callback(hObject,eventdata,handles)
% hObject    handle to zGo_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Command Stage Translation
% Make certain to enter a move to ... position in "z_mov_abs_display"
z_mov_abs_cmd=str2double(get(handles.z_mov_abs_display,'String'));
z_stage_pos=move_abs_trans_stage(hObject,eventdata,handles,3,z_mov_abs_cmd);
set(handles.z_pos_display,'String',num2str(z_stage_pos));
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes on button press in thetaGo_pushbutton.
function thetaGo_pushbutton_Callback(hObject,eventdata,handles)
% hObject    handle to thetaGo_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Command Stage Rotation
% Make certain to enter a move to ... position in "theta_mov_abs_display"
% theta_mov_abs_cmd=str2double(get(handles.theta_mov_abs_display,'String'));
% theta_stage_pos=move_abs_rotation_stage(theta_mov_abs_cmd);
% set(handles.theta_pos_display,'String',num2str(theta_stage_pos));
guidata(hObject,handles);  % Update GUI handles structure

function readData_indicator_Callback(hObject,eventdata,handles)
% hObject    handle to readData_indicator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of readData_indicator as text
%        str2double(get(hObject,'String')) returns contents of readData_indicator as a double
set(handles.readData_indicator,'BackgroundColor','red');  % set(handles.readData_indicator,'String','on');
% pause(1);  % Blink cell or box (sec)
% set(handles.readData_indicator,'String','off');
set(handles.readData_indicator,'BackgroundColor','white');
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function readData_indicator_CreateFcn(hObject,eventdata,handles)
% hObject    handle to readData_indicator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function freq_disp_Callback(hObject,eventdata,handles)
% hObject    handle to freq_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of freq_disp as text
%        str2double(get(hObject,'String')) returns contents of freq_disp as a double
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes during object creation, after setting all properties.
function freq_disp_CreateFcn(hObject,eventdata,handles)
% hObject    handle to freq_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in disconnectCOMports_pushbutton.
function disconnectCOMports_pushbutton_Callback(hObject,eventdata,handles)
% hObject    handle to disconnectCOMports_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes on button press in pwr2sensorCkt_checkbox.
function pwr2sensorCkt_checkbox_Callback(hObject,eventdata,handles)
% hObject    handle to pwr2sensorCkt_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of pwr2sensorCkt_checkbox
guidata(hObject,handles);  % Update GUI handles structure

function stage_pos=move_abs_trans_stage(hObject,eventdata,handles,stage_num,cmd_pos)
% Move Translation Stage to Commanded Position:
% 1 - x-stage
% 2 - y-stage
% 3 - z-stage
switch stage_num
    case 1  % x-Stage
        pos_str=['1PA', num2str(cmd_pos)];  % Position String: pos_str
        fprintf(handles.s,pos_str);
        % fwrite(handles.s,13);
        a=0;
        while a==0
            fprintf(handles.s,'1MD?');  % Motion Done Status
            a=fscanf(handles.s,'%d');
            fprintf(handles.s,'1TP');  % Instantaneous Time Position
            stage_pos_str=fscanf(handles.s,'%f');  % Character String: stage_pos_str
            formatted_stage_pos_str=sprintf('%.4f',stage_pos_str);
            set(handles.x_pos_display,'String',formatted_stage_pos_str);
            pause(0.05);
        end
        fprintf(handles.s,'1PA?');
        stage_pos_str=fscanf(handles.s,'%f');  % Character String: stage_pos_str
        formatted_stage_pos_str=sprintf('%.4f',stage_pos_str);
        set(handles.x_pos_display,'String',formatted_stage_pos_str);
    case 2  % y-Stage
        pos_str=['2PA',num2str(cmd_pos)];  % Position String: pos_str
        fprintf(handles.s,pos_str);
        % fwrite(handles.s,13);
        a=0;
        while a==0
            fprintf(handles.s,'2MD?');  % Motion Done Status
            a=fscanf(handles.s,'%d');
            fprintf(handles.s,'2TP');  % Instantaneous Time Position
            stage_pos_str=fscanf(handles.s,'%f');  % Character String: stage_pos_str
            formatted_stage_pos_str=sprintf('%.4f',stage_pos_str);
            set(handles.y_pos_display,'String',formatted_stage_pos_str);
            pause(0.05);
        end
        fprintf(handles.s,'2PA?');
        stage_pos_str=fscanf(handles.s,'%f');  % Character String: stage_pos_str
        formatted_stage_pos_str=sprintf('%.4f',stage_pos_str);
        set(handles.y_pos_display,'String',formatted_stage_pos_str);
    case 3  % z-Stage
        pos_str=['3PA',num2str(cmd_pos)];  % Position String: pos_str
        fprintf(handles.s,pos_str);
        % fwrite(handles.s,13);
        a=0;
        while a==0
            fprintf(handles.s,'3MD?');  % Motion Done Status
            a=fscanf(handles.s,'%d');
            fprintf(handles.s,'3TP');  % Instantaneous Time Position
            stage_pos_str=fscanf(handles.s,'%f');  % Character String: stage_pos_str
            formatted_stage_pos_str=sprintf('%.4f',stage_pos_str);
            set(handles.z_pos_display,'String',formatted_stage_pos_str);
            % pause(0.05);
        end
        fprintf(handles.s,'3PA?');
        stage_pos_str=fscanf(handles.s,'%f');  % Character String: stage_pos_str
        formatted_stage_pos_str=sprintf('%.4f',stage_pos_str);
        set(handles.z_pos_display,'String',formatted_stage_pos_str);
end
stage_pos=str2double(formatted_stage_pos_str);
guidata(hObject,handles);  % Update GUI handles structure

function [stage_x_pos,stage_y_pos,stage_z_pos]=home_all_trans_stages(hObject,eventdata,handles,x_cmd,y_cmd,z_cmd)
% Home All Translation Stages Simultaneously: 1 - x-stage , 2 - y-stage , 3 - z-stage
% Command Positions for Stages to get to ...
% x_cmd , y_cmd , z_cmd
x_pos_str=['1PA',num2str(x_cmd)];  % Position String: x_pos_str
fprintf(handles.s,x_pos_str);
% fwrite(handles.s,13);
fprintf(handles.s,'1PA?');
x_stage_pos_str=fscanf(handles.s,'%f');  % Character String: x_stage_pos_str
formatted_x_stage_pos_str=sprintf('%.4f',x_stage_pos_str);
set(handles.x_pos_display,'String',formatted_x_stage_pos_str);

y_pos_str=['2PA',num2str(y_cmd)];  % Position String: y_pos_str
fprintf(handles.s,y_pos_str);
% fwrite(handles.s,13);
fprintf(handles.s,'2PA?');
y_stage_pos_str=fscanf(handles.s,'%f');  % Character String: y_stage_pos_str
formatted_y_stage_pos_str=sprintf('%.4f',y_stage_pos_str);
set(handles.y_pos_display,'String',formatted_y_stage_pos_str);

z_pos_str=['3PA',num2str(z_cmd)];  % Position String: z_pos_str
fprintf(handles.s,z_pos_str);
% fwrite(handles.s,13);
fprintf(handles.s,'3PA?');
z_stage_pos_str=fscanf(handles.s,'%f');  % Character String: z_stage_pos_str
formatted_z_stage_pos_str=sprintf('%.4f',z_stage_pos_str);
set(handles.z_pos_display,'String',formatted_z_stage_pos_str);
% Continuous Update to Display
a=0;
while a~=3
    fprintf(handles.s,'1MD?'); i=fscanf(handles.s,'%d');  % MD (Motion Done) Status
    fprintf(handles.s,'2MD?'); j=fscanf(handles.s,'%d');
    fprintf(handles.s,'3MD?'); k=fscanf(handles.s,'%d');
    fprintf(handles.s,'1TP');
    x_pos_str=fscanf(handles.s,'%f');
    formatted_x_pos_str=sprintf('%.4f',x_pos_str);
    set(handles.x_pos_display,'String',formatted_x_pos_str);
    fprintf(handles.s,'2TP');
    y_pos_str=fscanf(handles.s,'%f');
    formatted_y_pos_str=sprintf('%.4f',y_pos_str);
    set(handles.y_pos_display,'String',formatted_y_pos_str);
    fprintf(handles.s,'3TP');
    z_pos_str=fscanf(handles.s,'%f');
    formatted_z_pos_str=sprintf('%.4f',z_pos_str);
    set(handles.z_pos_display,'String',formatted_z_pos_str);
    % pause(0.05);
    a=i+j+k;
    if a>3, break,end
end
stage_x_pos=str2double(formatted_x_pos_str);
stage_y_pos=str2double(formatted_y_pos_str);
stage_z_pos=str2double(formatted_z_pos_str);
guidata(hObject,handles);  % Update GUI handles structure

function stage_pos=move_abs_rotation_stage(cmd_mode)
% Stepper Motor Command Mode (cmd_mod)
% 1) cmd_mode=0  CCW Single-Step
% 2) cmd_mode=1  CW  Single-Step
% 3) cmd_mode=2                                                                                                                                          Home (Return Stage to Start Position or Step Counter Origin)


% --- Executes when selected object is changed in scan_mode_radiobutton.
function scan_mode_radiobutton_SelectionChangeFcn(hObject,eventdata,handles)
% hObject    handle to the selected object in scan_mode_radiobutton
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
% switch get(eventdata.NewValue,'Tag')
%     case'xyz_theta_axis'
%         handles.scan_mode=0;
%     case'z_theta_axis'
%         handles.scan_mode=1;
% Select Scan Mode:
% 1) xyz,theta axis (or xp,zp,theta axis) Scanning
% 2) z,theta axis (or zp,theta axis) Scanning
% handles.scan_mode=0  for xyz,theta axis Scanning
% handles.scan_mode=1  for z,theta axis Scanning
guidata(hObject,handles);  % Update GUI handles structure


% --- Executes on button press in getData_pushbutton.
function getData_pushbutton_Callback(hObject,eventdata,handles)
% hObject    handle to getData_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
time_test=tic;
% Blink Read Data Indicator "red"
% readData_indicator_Callback(hObject,eventdata,handles);  % Blink/Flash "Reading Data" Display
% Data Sampling Rate: smpl_rate (samples/sec)
% smpl_tme=1/smpl_rate;  % Sample Time: smpl_tme (sec/sample)
% delta_tme=smpl_tme;  % Time Step Increment: delta_tme (sec/sample)
% tme_min=0; tme_max=tme_min+(num_data_smpls-1)*delta_tme;
% tme=tme_min:delta_tme:tme_max;  % tme=(0:num_data_smpls-1)*smpl_tme
% Probe (prb) Data
raw_data_prb_array=handles.rx;  % raw_data_prb_array(1,1:handles.num_data_smpls)  Complex Single Precision Col Vector Array
MaxOutputVoltage_prb=max(real(raw_data_prb_array));
set(handles.prbVolts_text,'String',num2str(1000*MaxOutputVoltage_prb));
mag_prb_array=abs(raw_data_prb_array);

% Fourier Transform (FT) Method

% if (length_tx < length_rx)
%     c = [ zeros(1,length_rx-1) handles.txChips' zeros(1,length_rx-length_tx) ];
%     d = [ raw_data_prb_array zeros(1,length_rx-1) ];
% else
%     c = [ zeros(1,length_tx-1) handles.txChips' ];
%     d = [ raw_data_prb_array zeros(1,length_tx-length_rx+length_tx-1) ];
% end
% 
% NFFT_c = 2^nextpow2(c);
% NFFT_d = 2^nextpow2(d);
% % calculate crosscorrelation
% e = fft(c,NFFT_c);
% f = fft(d,NFFT_d);
% fft_corr= e.*conj(f);
% ift_corr = ifft(fft_corr);
% [~,max_ift_indx] = max(abs(ift_corr));
% delay_fft = max_ift_indx/handles.smpl_rate;


length_rx = length(raw_data_prb_array);
NFFT_rx = 2^nextpow2(length_rx);
rx_win = hanning(length_rx,'periodic');
fft_rx = fft(raw_data_prb_array.*rx_win,NFFT_rx);

% phase difference calculation
[~,fftindx] = max(abs(fft_rx));
PhDiff_fft = angle(fft_rx(fftindx));

% data_rxx=handles.win.*handles.rxChips';  % Data Window: data_rxx(1:handles.num_data_smpls,1) complex col vector single precision (16-bit)
% num_fft_smpls=1024;  % num_fft_smpls=2*(handles.num_nyq_smpls-1);
% fft_daqData=fft(real(data_rxx),num_fft_smpls);  % FFT of Sensor Probe (prb) real part data
% % indx_prb=handles.indx_at_pwr_spec_pk_prb;  % Index Pointer
% indx_prb=16;
% Power Spectrum: pwr_spec
% Complex Variable: z = x + 1i*y
% pwr_spec = abs(z).^2
%          = real(z).^2 + imag(z).^2
%          = z .* conj(z)
% magnitude = sqrt(pwr_spec) = abs(z)

% At Peak in Power Spectrum: indx_at_pwr_spec_pk
% Magnitude: |Y| (rms volts)
mag_ft_scalar=abs(fft_rx(fftindx));  % mag_ft_scalar=sqrt(pwr_spec_prb(indx_prb))
phz_ft_scalar=PhDiff_fft;

% phz_ft_scalar=atan2(imag(fft_daqData(indx_prb)),real(fft_daqData(indx_prb)));  % Phase (rad)

% In-Phase Quadrature (IQ) Method
% phz_prb_array=atan2(imag(raw_data_prb_array),real(raw_data_prb_array));  % Phase (rad)
phz_iq_scalar = mean(angle(raw_data_prb_array));
% phz_prb_array=angle(raw_data_prb_array)-angle(handles.txChips);  % Phase (rad)
mag_iq_scalar=mean(mag_prb_array);
% phz_iq_scalar=mean(phz_prb_array);  % Phase (rad)


% Magnitude and Phase in GUI Display
% Placed in order of appearance
h_selected_processing_method=get(handles.daq_mode,'SelectedObject');
seletected_method_tag=get(h_selected_processing_method,'tag');
switch seletected_method_tag
    case 'ft_radiobutton'
        handles.data_processing_method=0;
        % Fourier Transform (FT) Method
        set(handles.phzPrb_text,'String',num2str(phz_ft_scalar));  % phz_prb (rad)
        set(handles.magPrb_text,'String',num2str(mag_ft_scalar*1000));  % mag_prb (mVolts)
    case 'iq_radiobutton'
        handles.data_processing_method=1;
        % In-Phase Quadrature (IQ) Method
        set(handles.phzPrb_text,'String',num2str(phz_iq_scalar));  % phz_prb (rad)
        set(handles.magPrb_text,'String',num2str(mag_iq_scalar*1000));  % mag_prb (mVolts)
end


% Evaluate Arrays from Workspace
% Scan
xp=evalin('base','Xp');
zp=evalin('base','Zp');
th=evalin('base','Theta');
% Outputs
mag_ft=evalin('base','Magn_FT');
phz_ft=evalin('base','Phze_FT');
mag_iq=evalin('base','Magn_IQ');
phz_iq=evalin('base','Phze_IQ');

% Get GUI Display Scan Positions/Displacements
x1_pos=str2double(get(handles.x_pos_display,'String'));
x2_pos=str2double(get(handles.y_pos_display,'String'));
zp_pos=str2double(get(handles.z_pos_display,'String'));
th_pos=str2double(get(handles.theta_pos_display,'String'));

% Scan Arrays
xp(handles.indx_xp,handles.indx_zp,handles.indx_theta)=x1_pos+x2_pos;
zp(handles.indx_xp,handles.indx_zp,handles.indx_theta)=zp_pos;
th(handles.indx_xp,handles.indx_zp,handles.indx_theta)=th_pos;

% Output Arrays
mag_ft(handles.indx_xp,handles.indx_zp,handles.indx_theta)=mag_ft_scalar;
phz_ft(handles.indx_xp,handles.indx_zp,handles.indx_theta)=phz_ft_scalar;
mag_iq(handles.indx_xp,handles.indx_zp,handles.indx_theta)=mag_iq_scalar;
phz_iq(handles.indx_xp,handles.indx_zp,handles.indx_theta)=phz_iq_scalar;


% Assign Arrays to Workspace
% Scan
assignin('base','Xp',xp);
assignin('base','Zp',zp);
assignin('base','Theta',th);
% Outputs
assignin('base','Magn_FT',mag_ft);
assignin('base','Phze_FT',phz_ft);
assignin('base','Magn_IQ',mag_iq);
assignin('base','Phze_IQ',phz_iq);

% Scrolling (streaming) Data Strip-Chart Plot
% Array Size: mag_prb_plt(1:handles.numWinPnts)
mag_ft_plt=evalin('base','Mag_FT_Plt');
phz_ft_plt=evalin('base','Phz_FT_Plt');
mag_iq_plt=evalin('base','Mag_IQ_Plt');
phz_iq_plt=evalin('base','Phz_IQ_Plt');

% Strip-Chart Plot Shift-Register of Data: Array Size (1:handles.numWinPnts)
% Probe Magnitude
mag_ft_plt(1:handles.numWinPnts-1)=mag_ft_plt(2:handles.numWinPnts);  % Shift-Register
mag_ft_plt(handles.numWinPnts)=mag_ft_scalar;
mag_iq_plt(1:handles.numWinPnts-1)=mag_iq_plt(2:handles.numWinPnts);  % Shift-Register
mag_iq_plt(handles.numWinPnts)=mag_iq_scalar;
% Probe Phase
phz_ft_plt(1:handles.numWinPnts-1)=phz_ft_plt(2:handles.numWinPnts);  % Shift-Register
phz_ft_plt(handles.numWinPnts)=phz_ft_scalar;
phz_iq_plt(1:handles.numWinPnts-1)=phz_iq_plt(2:handles.numWinPnts);  % Shift-Register
phz_iq_plt(handles.numWinPnts)=phz_iq_scalar;

% Assign Scrolling Data Strp-Chart Plot Arrays to Workspace
assignin('base','Mag_FT_Plt',mag_ft_plt);
assignin('base','Phz_FT_Plt',phz_ft_plt);
assignin('base','Mag_IQ_Plt',mag_iq_plt);
assignin('base','Phz_IQ_Plt',phz_iq_plt);


% Strip-Chart Scrolling Data Plot
% Toggle FT/IQ Data Processing Method for Strip-Chart Scrolling (or Streamed) Data Plots
% LHS Plot: Magnitude
% RHS Plot: Phase
if handles.data_processing_method==0
    % *** Fourier Transform (FT) Method ***
    % *** Magnitude ***
    % Plot (freq peak in power spectrum due to) Magnitude for Probe (prb)
    % axesHandle=findobj(gcf,'Tag','axes_magn_plt');  % set(axesHandle,'Tag','axes_magn_plt');
    plot(1000*handles.time_plt,1000*mag_ft_plt(1:handles.numWinPnts)/1,'Color','r','LineStyle','-','LineWidth',1.0,'Parent',handles.axes_magn_plt);  % LineWidth=0.5 (default)
    % hold(handles.axes_magn_plt,'on');  % set(handles.axes_magn_plt,'NextPlot','Add');  % hold(handles.axes_magn_plt,'off');
    title('\bfFourier Transform (FT) Method Magnitude','FontSize',14,'Parent',handles.axes_magn_plt);  % FontSize=12 (default)
    xlabel('Strip-Chart Time Tick per Sample-Step','FontSize',12,'Parent',handles.axes_magn_plt);  % FontSize=12 (default)
    ylabel('Magnitude (mVolts)','FontSize',12,'Parent',handles.axes_magn_plt);  % FontSize=12 (default)
    % legend(handles.axes_magn_plt,'FreqA','FreqB','Location','NorthWest');  % Location=NorthEast (default)
    % ylim(handles.axes_magn_plt,[0 +600]);
    % ylim(handles.axes_magn_plt,[0 +2000]);  % FT method high magnitude value
    % set(gca,'YLimMode','manual');
    axis([0 5000 0 1000]);
    % guidata(hObject,handles.axes_magn_plt);  % Update handles structure
    
    % *** Phase ***
    % Plot (freq peak in power spectrum due to) Phase for Probe (prb)
    plot(1000*handles.time_plt,phz_ft_plt(1:handles.numWinPnts)/1,'Color','b','LineStyle','-','LineWidth',1.0,'Parent',handles.axes_phze_plt);  % LineWidth=0.5 (default)
    % hold(handles.axes_phze_plt,'on');  % set(handles.axes_phze_plt,'NextPlot','Add');  % hold(handles.axes_phze_plt,'off');
    title('\bfFourier Transform (FT) Method Phase','FontSize',14,'Parent',handles.axes_phze_plt);  % FontSize=12 (default)
    xlabel('Strip-Chart Time Tick per Sample-Step','FontSize',12,'Parent',handles.axes_phze_plt);  % FontSize=12 (default)
    ylabel('Phase (rad)','FontSize',12,'Parent',handles.axes_phze_plt);  % FontSize=12 (default)
    % legend(handles.axes_phze_plt,'FreqA','FreqB','Location','NorthWest');  % Location=NorthEast (default)
    % ylim(handles.axes_phze_plt,[-2*pi +2*pi]);  % Note: atan2(imag(z),real(z)) = -pi to +pi range
    % set(gca,'YLimMode','manual');
    axis([0 5000 -2*pi +2*pi]);
    
    
    % guidata(hObject,handles.axes_phze_plt);  % Update handles structure
elseif handles.data_processing_method==1
    % *** In-Phase Quadrature (IQ) Method ***
    % *** Magnitude ***
    % Plot (freq peak in power spectrum due to) Magnitude for Probe (prb)
    % axesHandle=findobj(gcf,'Tag','axes_magn_plt');  % set(axesHandle,'Tag','axes_magn_plt');
    plot(1000*handles.time_plt,1000*mag_iq_plt(1:handles.numWinPnts)/1,'Color','r','LineStyle','-','LineWidth',1.0,'Parent',handles.axes_magn_plt);  % LineWidth=0.5 (default)
    % hold(handles.axes_magn_plt,'on');  % set(handles.axes_magn_plt,'NextPlot','Add');  % hold(handles.axes_magn_plt,'off');
    title('\bfIn-Phase Quadrature (IQ) Magnitude','FontSize',14,'Parent',handles.axes_magn_plt);  % FontSize=12 (default)
    xlabel('Strip-Chart Time Tick per Sample-Step','FontSize',12,'Parent',handles.axes_magn_plt);  % FontSize=12 (default)
    ylabel('Magnitude (mVolts)','FontSize',12,'Parent',handles.axes_magn_plt);  % FontSize=12 (default)
    % legend(handles.axes_magn_plt,'FreqA','FreqB','Location','NorthWest');  % Location=NorthEast (default)
    % ylim(handles.axes_magn_plt,[0 +600]);
    % set(gca,'YLimMode','manual');
    axis([0 5000 0 1000]);
    % guidata(hObject,handles.axes_magn_plt);  % Update handles structure
    
    % *** Phase ***
    % Plot (freq peak in power spectrum due to) Phase for Probe (prb)
    plot(1000*handles.time_plt,phz_iq_plt(1:handles.numWinPnts)/1,'Color','b','LineStyle','-','LineWidth',1.0,'Parent',handles.axes_phze_plt);  % LineWidth=0.5 (default)
    % hold(handles.axes_phze_plt,'on');  % set(handles.axes_phze_plt,'NextPlot','Add');  % hold(handles.axes_phze_plt,'off');
    title('\bfIn-Phase Quadrature (IQ) Phase','FontSize',14,'Parent',handles.axes_phze_plt);  % FontSize=12 (default)
    xlabel('Strip-Chart Time Tick per Sample-Step','FontSize',12,'Parent',handles.axes_phze_plt);  % FontSize=12 (default)
    ylabel('Phase (rad)','FontSize',12,'Parent',handles.axes_phze_plt);  % FontSize=12 (default)
    % legend(handles.axes_phze_plt,'FreqA','FreqB','Location','NorthWest');  % Location=NorthEast (default)
    % ylim(handles.axes_phze_plt,[-2*pi +2*pi]);  % Note: atan2(imag(z),real(z)) = -pi to +pi range
    % set(gca,'YLimMode','manual');
    axis([0 5000 -2*pi +2*pi]);
    % guidata(hObject,handles.axes_phze_plt);  % Update handles structure
end
drawnow;
% Close-Out Commands
time_test_done=toc(time_test);
disp(time_test_done);
guidata(hObject,handles);  % Update GUI handles structure
% end of getData_callback

% --- Executes on button press in startScan_pushbutton.
function startScan_pushbutton_Callback(hObject,eventdata,handles)
% hObject    handle to startScan_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Start Scan Steps
% ****************
% 1) Home All Stages (Re-Initialize or Reset in case of Re-Start or Multiple Scanning Start Intitiations)
% 2) Wait for Stages to arrive home or to min positions (motion done status)
% 3) Get Data at home position
% 4) Scanning Loop - Scan Mode: 1) xyz,theta axis or 2) z,theta axis
% 5)   Command: Send Stage to move to new position: x = x + delta_x
% 6)   Wait for Stage to arrive at new position (motion done status)
% 7)   Get Data at new scan position
% 8)   Evaluate new data by transferring from Workspace into Scanning function
% 9)   Save Data (Save data file within loop while scanning in order to recover partial data in case of crash)
%        a) Scan Position
%        b) Sensor Probe (prb) Signal Data
% 10) Loop Return

% Min Values
x1Min=str2double(get(handles.xMin_display,'String'));
x2Min=str2double(get(handles.yMin_display,'String'));
zpMin=str2double(get(handles.zMin_display,'String'));
thetaMin=str2double(get(handles.thetaMin_display,'String'));
% Max Values
x1Max=str2double(get(handles.xMax_display,'String'));
x2Max=str2double(get(handles.yMax_display,'String'));
zpMax=str2double(get(handles.zMax_display,'String'));
thetaMax=str2double(get(handles.thetaMax_display,'String'));
% Num of Steps
num_x1Steps=str2double(get(handles.num_xSteps_display,'String'));
num_x2Steps=str2double(get(handles.num_ySteps_display,'String'));
num_zpSteps=str2double(get(handles.num_zSteps_display,'String'));
num_thetaSteps=str2double(get(handles.num_thetaSteps_display,'String'));

% Compute Step-Size or delta_x
x1StepSize=(x1Max-x1Min)/num_x1Steps;  % x1StepSize=delta_x1
x2StepSize=(x2Max-x2Min)/num_x2Steps;  % x2StepSize=delta_x2
zpStepSize=(zpMax-zpMin)/num_zpSteps;  % zStepSize=delta_z
thetaStepSize=(thetaMax-thetaMin)/num_thetaSteps;  % thetaStepSize=delta_theta
% Note: Set x1StepSize=x2StepSize for Tandem x-axis and y-axis Stages

% Display Step-Size
set(handles.xStepSize_display,'String',num2str(x1StepSize));
set(handles.yStepSize_display,'String',num2str(x2StepSize));
set(handles.zStepSize_display,'String',num2str(zpStepSize));
set(handles.thetaStepSize_display,'String',num2str(thetaStepSize));

% Move All Translation Stages to Minimum Position
[x_start_pos,y_start_pos,z_start_pos]=home_all_trans_stages(hObject,eventdata,handles,x1Min,x2Min,zpMin);
set(handles.x_pos_display,'String',num2str(x_start_pos));
set(handles.y_pos_display,'String',num2str(y_start_pos));
set(handles.z_pos_display,'String',num2str(z_start_pos));

% Rotation Stage is Preset
% theta_min=str2num(get(handles.thetaMin_display,'String'));
% theta_max=str2num(get(handles.thetaMax_display,'String'));

% Scanning Cycle Loop for Translation and Rotation Stage Step-Stare Data Acquisition
% Scanning Loop Scanning Variables:
% 1) xyz,theta axis or xp,zp,theta axis Scanning (handles.scan_mode=0)
% 2) z,theta axis or zp,theta axis Scanning      (handles.scan_mode=1)

% Index Pointers to an Array Indicie Step Counter
handles.indx_xp=1;
handles.indx_zp=1;
handles.indx_theta=1;
theta_dir=+1;  % Rotation Stage angular direction

% Current Positions (pos) following step increment/decrement
% x1_pos=str2double(get(handles.x_pos_display,'String'));
% x2_pos=str2double(get(handles.y_pos_display,'String'));
% zp_pos=str2double(get(handles.z_pos_display,'String'));
% theta_pos=str2double(get(handles.theta_pos_display,'String'));

% Initialization
zp_pos=zpMin;
theta=thetaMin;
%% Generate Barker Code Sequence for spread
hBCode = comm.BarkerCode('Length', 13,...
               'SamplesPerFrame', 30);
seq_barker = step(hBCode);
data_tx = [seq_barker;ones(970,1)];

% Sending and Receving Signal using Ettus
handles.stopFlag=0;

while handles.stopFlag==0
    % Initialize Data Acquisition (DAQ) Stream Processing Loop
    radio=handles.radio;
    
    t_start_scan=tic;
    %%  Stream Processing Loop
        if (strcmp(radio.Status,'Success'))
        % disp('USRP/SDRu IP Address Connection Status: Successful');
        timeCounter=0;  % Initialization
        while timeCounter < handles.paramsTxStopTime  % Loop until transmitter stop time
            % --- Transmitter Stream Processing ---
            step(handles.hSDRuTx, data_tx); % Transmit to USRP(R) radio
            % Workspace Variable
            % data_tx=double(zeros(handles.num_data_smpls,1));  % Double Precision
            
            % --- Receiver Stream Processing ---
            % Get baseband samples from USRP(R) Board
            [data_rx, len_rx] = step(handles.hSDRuRx);
            % Workspace Variables
            % data_rx=single(zeros(handles.num_data_smpls,1) + 1i*zeros(handles.num_data_smpls,1));  % Complex Single Precision
            % data_rx(1:length_rx,1) complex single row vector
            % length_rx = handles.num_data_smpls = length(data_rx) = 1024 samples or data points
            % handles.data_rx=data_rx;  % only record here
            timeCounter=timeCounter+handles.paramsRxUSRPFrameTime;
        end
    else
        warning(message('USRP/SDRu IP Address Connection Status: Not Successful'));
    end % one signal streaming loop
    
    % set(handles.theta_pos_display,'String',num2str(theta));  % display theta position  
    
        %% Align the Received with transmitted signal
    [acor,Pos] = xcorr(data_rx,data_tx); % y slower than x >0
%     acor = acor(Pos>0);
%     Pos = Pos(Pos>0);
    [~,I] = max(abs(acor));
    lagDiff = Pos(I);
    delayDiff = lagDiff/(handles.smpl_rate);
    Synchro = sprintf('Uncertain delay between Tx and Rx is %d second',delayDiff);
    disp(Synchro);
%     acor_val= double(abs(acor));
%     figure;
%     plot(Pos*1e3/Fs,acor_val)
%     title('Delay in mircrosecond of Recevier side')
%     xlabel('Delay[\mus]')
%     ylabel('Cross correlation of received and transmitted sequence')
    
    % Synchronized received signal
    data_rxSync = data_rx(lagDiff+1+length(seq_barker):end);
    data_missing = ones((handles.num_data_smpls-length(data_rxSync)),1);
    
    % Send over to callback function
    handles.rx = [data_rxSync;data_missing];
    getData_pushbutton_Callback(hObject,eventdata,handles);  % Get Magnitude/Phase Data anyways
    
    
    %% Rotation Stage Motion Control after each scanning
    if (handles.indx_zp <= num_zpSteps+1)
        % Vertical Translation Loop: z(k)
        if (handles.indx_theta==1 && theta_dir==-1) || (handles.indx_theta==num_thetaSteps+1 && theta_dir==+1)
            % Stepper Motor Rotation Loop: theta(j)
            % getData_pushbutton_Callback(hObject,eventdata,handles);  % Get Magnitude/Phase Data
            
            % Scan Data Variables
            xp=evalin('base','Xp');
            zp=evalin('base','Zp');
            th=evalin('base','Theta');
            
            % Output Variables
            mag_ft=evalin('base','Magn_FT');
            phz_ft=evalin('base','Phze_FT');
            mag_iq=evalin('base','Magn_IQ');
            phz_iq=evalin('base','Phze_IQ');
            
            % Increment indx_zp and Change/Reverse theta Direction
            handles.indx_zp=handles.indx_zp+1;  % Increment Counter
            theta_dir=-theta_dir; % scan in the other direction
            zp_pos=zp_pos+zpStepSize;
            % Move z-axis Translation Stage to next position
            % stage_pos=move_abs_trans_stage(hObject,eventdata,handles,stage_num,cmd_pos);
            t_stage_move=tic;
            zp_pos_move=move_abs_trans_stage(hObject,eventdata,handles,3,zp_pos);
            set(handles.z_pos_display,'String',num2str(zp_pos_move));
            t_stage_move_done=toc(t_stage_move);
            disp(t_stage_move_done);
            
        else
            % Stepper Motor Rotation Loop: theta(j)
            % getData_pushbutton_Callback(hObject,eventdata,handles);  % Get Magnitude/Phase Data
            
            % Scan Data Variables
            xp=evalin('base','Xp');
            zp=evalin('base','Zp');
            th=evalin('base','Theta');
            
            % Output Variables
            mag_ft=evalin('base','Magn_FT');
            phz_ft=evalin('base','Phze_FT');
            mag_iq=evalin('base','Magn_IQ');
            phz_iq=evalin('base','Phze_IQ');
            
            % Increment/Decrement theta and Change/Reverse theta Direction
            theta=theta+theta_dir*thetaStepSize;
            handles.indx_theta=handles.indx_theta+theta_dir;  % Update (Increment/Decrement) Counter
            % pause(2);  % pause(10) (sec)
            % Move Rotation Stage - Stepper Motor to next position (Single-Step)
            % fprintf(handles.s,['DIR=%d ST0EP=%d'],[handles.dir 10],'mode','sync');
            if theta_dir==-1
                dir=0;  % CCW dir of rotation
            elseif theta_dir==+1
                dir=+1;  % CW dir of rotation
            end
            
            % Rotation Stage - Stepper Motor move command
            % fprintf(handles.ser,'DIR=%d',4,'sync');  % Scanning Mode
            fprintf(handles.ser,'DIR=%d',dir,'sync');
            set(handles.theta_pos_display,'String',num2str(theta));  % display theta position
%             if (strcmp(radio.Status,'Success'))
%                 % disp('USRP/SDRu IP Address Connection Status: Successful');
%                 timeCounter=0;  % Initialization
%                 while timeCounter < handles.paramsTxStopTime  % Loop until transmitter stop time
%                     % --- Transmitter Stream Processing ---
%                     step(handles.hSDRuTx,data_tx); % Transmit to USRP(R) Radio
%                     % Workspace Variable
%                     % data_tx=double(zeros(handles.num_data_smpls,1));  % Double Precision
%                     
%                     % --- Receiver Stream Processing ---
%                     [data_rrx,length_rrx]=step(handles.hSDRuRx);  % [data_rx,~]=step(handles.hSDRuRx)
%                     % Workspace Variables
%                     % data_rx=single(zeros(handles.num_data_smpls,1) + 1i*zeros(handles.num_data_smpls,1));  % Complex Single Precision
%                     % data_rx(1:length_rx,1) complex single row vector
%                     % length_rx = handles.num_data_smpls = length(data_rx) = 1024 samples or data points
%                     timeCounter=timeCounter+handles.paramsRxUSRPFrameTime;
%                 end
%             else
%                 warning(message('USRP/SDRu IP Address Connection Status: Not Successful'));
%             end % one signal streaming loop
            
            
            %             % fprintf(handles.ser,'DIR=%d',2,'sync'); % stop the motor and collect data
            %             if (strcmp(radio.Status,'Success'))
            %                 % disp('USRP/SDRu IP Address Connection Status: Successful');
            %                 timeCounter=0;  % Initialization
            %                 while timeCounter < handles.paramsTxStopTime  % Loop until transmitter stop time
            %                     % --- Transmitter Stream Processing ---
            %                     step(handles.hSDRuTx,data_tx); % Transmit to USRP(R) Radio
            %                     % Workspace Variable
            %                     % data_tx=double(zeros(handles.num_data_smpls,1));  % Double Precision
            %
            %                     % --- Receiver Stream Processing ---
            %                     [data_rrrx,length_rrrx]=step(handles.hSDRuRx);  % [data_rx,~]=step(handles.hSDRuRx)
            %                     % Workspace Variables
            %                     % data_rx=single(zeros(handles.num_data_smpls,1) + 1i*zeros(handles.num_data_smpls,1));  % Complex Single Precision
            %                     % data_rx(1:length_rx,1) complex single row vector
            %                     % length_rx = handles.num_data_smpls = length(data_rx) = 1024 samples or data points
            %                     timeCounter=timeCounter+handles.paramsRxUSRPFrameTime;
            %                 end
            %             else
            %                 warning(message('USRP/SDRu IP Address Connection Status: Not Successful'));
            %             end % one signal streaming loop
            
            % set(handles.theta_pos_display,'String',num2str(theta));
            drawnow;
            % pause(0.05);
        end
        save('-mat','Ultrasonic_Ettus_Bck_Scttrng_Scanner_Data_1Freq_1.mat','xp','zp','th', ...
            'mag_ft','phz_ft','mag_iq','phz_iq');  % Save after each horizontal scan (fwd & bck sweep)
        t_scan_time=toc(t_start_scan);
        disp(t_scan_time);
    else
        handles.stopFlag=1;
    end
    % stopScan_pushbutton_Callback(hObject,eventdata,handles);  % Returns rotation stage near home at origin
    
    if (handles.stopFlag==1)
        release(hSDRuTx);  % fclose(hSDRuTx);
        delete(hSDRuTx);
        clear hSDRuTx;
        release(hSDRuRx);  % fclose(hSDRuRx);
        delete(hSDRuRx);
        clear hSDRuRx;
        break
    end
end
disp('Mystery is Coming');
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes on button press in stopScan_pushbutton.
function stopScan_pushbutton_Callback(hObject,eventdata,handles)
% hObject    handle to stopScan_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Rotation Stage - Stepper Motor
fprintf(handles.ser,'DIR=%d',2,'sync');  % Stop
fprintf(handles.s,'1st;2st;3st');
handles.stopFlag=1;
fclose(handles.ser);
fclose(handles.s);
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes on button press in home_pushbutton.
function home_pushbutton_Callback(hObject,eventdata,handles)
% hObject    handle to home_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Rotation Stage - Stepper Motor
fprintf(handles.ser,'DIR=%d',3,'sync');  % Go to Home Position
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes when selected object is changed in daq_mode.
function daq_mode_SelectionChangeFcn(hObject,eventdata,handles)
% hObject    handle to the selected object in daq_mode
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue,'Tag')
    case 'ft_radiobutton'
        handles.data_processing_method=0;  % Fourier Transform (FT) Method
    case 'iq_radiobutton'
        handles.data_processing_method=1;  % In-Phase Quadrature (IQ) Method
end
guidata(hObject,handles);  % Update GUI handles structure

% --- Executes on button press in calib_pushbutton.
function calib_pushbutton_Callback(hObject,eventdata,handles)
% hObject    handle to calib_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Freq=str2double(get(handles.freq_disp,'String'));  % MHz (Typical: 1.5 MHz)

% Setup & Configure Ettus N210 USRP Data Acquisition Box
% Quadrature Output at 100 MSps = 100e6 samples/sec using LFTX and LFRX Daughterboards
% Center RF Frequency: Transmitter (Tx) and Receiver (Rx)
paramsTx.CenterFrequency=str2double(get(handles.freq_disp,'String'))*1e6;  % Hz (Typical: 1.5 MHz)
paramsRx.CenterFrequency=str2double(get(handles.freq_disp,'String'))*1e6;  % Hz (Typical: 1.5 MHz)

% Transmitter (Tx) Parameters
paramsTx.StopTime                 = 3;      % seconds (minimum)
paramsTx.USRPClockRate            = 100e6;  % Hz
paramsTx.USRPInterpolationFactor  = 400;    % xScale-Factor
paramsTx.USRPGain                 = 0;      % dB
paramsTx.USRPFrameLength          = 1000;   % samples or data points
paramsTx.FrontEndSampleRate       = paramsTx.USRPClockRate ...
    /paramsTx.USRPInterpolationFactor;

% Receiver (Rx) Parameters
paramsRx.USRPClockRate        = 100e6;  % ADC working sample in Hz
paramsRx.USRPDecimationFactor = 400;    % xScale-Factor (output sampling rate 200kHz)
paramsRx.USRPGain             = 0;      % dB
paramsRx.USRPFrameLength      = 1000;   % samples or data points
paramsRx.FrontEndSampleRate   = paramsRx.USRPClockRate...
    /paramsRx.USRPDecimationFactor;
paramsRx.USRPFrameTime        = paramsRx.USRPFrameLength ...
    /paramsRx.FrontEndSampleRate;  % seconds

handles.hSDRuTx=comm.SDRuTransmitter(...
    'IPAddress',              '192.168.10.2', ...
    'CenterFrequency',        paramsTx.CenterFrequency,...
    'Gain',                   paramsTx.USRPGain, ...
    'InterpolationFactor',    paramsTx.USRPInterpolationFactor);

handles.hSDRuRx=comm.SDRuReceiver(...
    'IPAddress',             '192.168.10.2', ...
    'CenterFrequency',       paramsRx.CenterFrequency,...
    'Gain',                  paramsRx.USRPGain, ...
    'DecimationFactor',      paramsRx.USRPDecimationFactor, ...
    'SampleRate',            paramsRx.FrontEndSampleRate, ...
    'FrameLength',           paramsRx.USRPFrameLength, ...
    'OutputDataType',        'single');

radio=findsdru(handles.hSDRuRx.IPAddress);  % Find USRP/SDRu Radio

% USRP/SDRu Radio Handles
handles.num_data_smpls=paramsRx.USRPFrameLength;  % handles.num_data_smpls=paramsTx.USRPFrameLength
handles.paramsTxStopTime=paramsTx.StopTime;
handles.paramsRxUSRPFrameTime=paramsRx.USRPFrameTime;
handles.smpl_rate=paramsRx.FrontEndSampleRate;
handles.radio=radio;

% Initialize Data Acquisition (DAQ) Stream Processing Loop
data_tx=ones(handles.num_data_smpls,1);  % Row Vector
% USRP/SDRu Radio Communication to IP Address Status
if (strcmp(radio.Status,'Success'))
    disp('USRP/SDRu IP Address Connection Status: Successful');
    timeCounter=0;  % Initialization
    while timeCounter < handles.paramsTxStopTime  % Loop until transmitter stop time
        % --- Transmitter Stream Processing ---
        step(handles.hSDRuTx,data_tx); % Transmit to USRP(R) Radio
        % Workspace Variable
        % data_tx=double(zeros(handles.num_data_smpls,1));  % Double Precision
        
        % --- Receiver Stream Processing ---
        [data_rx,len_rx]=step(handles.hSDRuRx);
        % Workspace Variable
        % data_rx=single(zeros(handles.num_data_smpls,1) + 1i*zeros(handles.num_data_smpls,1));  % Complex Single Precision
        % data_rx(1:len_rx,1) complex single row vector
        % len_rx = handles.num_data_smpls = 4000 amples or data points
        timeCounter=timeCounter+handles.paramsRxUSRPFrameTime;
    end
else
    warning(message('USRP/SDRu IP Address Connection Status: Not Successful'));
end

% Number of Data Samples per Trigger for Acquisition (samples or data points)
num_data_smpls=handles.num_data_smpls;
num_fft_smpls=2^(nextpow2(num_data_smpls));  % Number of FFT Samples: num_fft_smpls = 2^n = 64, 128, 256, 512, 1024 etc.
% num_fft_smpls = 2*(handles.num_nyq_smpls-1);
% Note: num_nyq_smpls = handles.num_nyq_smpls = num_fft_smpls/2+1
num_nyq_smpls=num_fft_smpls/2+1;  % Index Pointer to Nyquist Limit - Number of Nyquist Samples at Nyquist Limit: num_nyq_smpls
smpl_rate=handles.smpl_rate;  % Sample Rate: smpl_rate (samples/sec or Hz)
freq_nyq=smpl_rate/2;  % Nyquist Frequency: nyq_freq (cycles/sec or Hz)
% Note: Nyquist Frequency at 2 samples/cycle
% freq(1)=0;  % DC component
% delta_freq = (freq_nyq-freq(1))/(num_nyq_smpls-1) = freq_nyq/(num_nyq_smpls-1)
%            = (smpl_rate/2-0)/((num_fft_smpls/2+1)-1)
%            = (smpl_rate/2)/(num_fft_smpls/2)
%            = smpl_rate/num_fft_smpls
delta_freq=smpl_rate/num_fft_smpls;
% freq=freq(1):delta_freq:freq_nyq;  % freq=freq_nyq*linspace(0,1,num_nyq_smpls)
% freq(num_nyq_smpls)=freq_nyq
freq(1:num_nyq_smpls)=0;  % freq=zeros(1,num_nyq_smpls)
freq_sum=freq(1);  % DC Component: freq(1)
for i=1:num_nyq_smpls
    freq(i)=freq_sum;
    freq_sum=freq_sum+delta_freq;  % Update
end
% freq(1:num_nyq_smpls)
handles.num_nyq_smpls=num_nyq_smpls;

% Transmitter (Tx) Parameters for a Sinusoidal Wave
freq_tx=paramsTx.CenterFrequency;  % Transmitter Frequency: freq_tx (cycles/sec)
num_data_pnts=paramsTx.USRPFrameLength;  % Transmitter Data Points: num_data_pnts (samples or data points)
sample_rate=paramsTx.FrontEndSampleRate;  % Transmitter Sample Rate: sample_rate (samples/sec)
sample_tme=1/sample_rate;  % Sample Time or Time Step Interval: sample_tme = dt = delta_tme (sec/sample)
tme(1:num_data_pnts)=0;  % tme=zeros(1,num_data_pnts)
tme_sum=tme(1);  % Initialization
for i=1:num_data_pnts
    tme(i)=tme_sum;
    tme_sum=tme_sum+sample_tme;  % Update
end
% tme=tme_min:delta_tme:tme_max
raw_data_tx_array=sin(2*pi*freq_tx*tme);
% MaxOutputVoltage_tx=max(real(raw_data_tx_array));
% set(handles.prbVolts_text,'String',num2str(1000*MaxOutputVoltage_tx));

if handles.data_processing_method==0
    % Fourier Transform (FT) Method for Sensor Probe (prb)
    fft_daqData=fft(raw_data_tx_array,num_fft_smpls);  % FFT real part data
    % Power Spectrum: pwr_spec
    % Complex Variable: z = x + 1i*y
    % pwr_spec = abs(z).^2
    %          = real(z).^2 + imag(z).^2
    %          = z .* conj(z)
    % magnitude = sqrt(pwr_spec) = abs(z)
    pwr_spec=fft_daqData(1:num_nyq_smpls).*conj(fft_daqData(1:num_nyq_smpls));
    % Find Peak in Power Spectrum
    pwr_spec(1)=0;  % Set DC component (along with low frequencies) to zero
    %     if get(handles.pwr2sensorCkt_checkbox,'Value')
    %         % Label: "Pwr to Sensor Ckt"
    %         % Tag: handles.pwr2sensorCkt_checkbox
    %         % Checkbox is checked indicating Power is enabled to the Sensor Probe (default condition)
    %         % Find High Peak in Power Spectrum - Probe (pwr_spec_prb)
    %         [~,indx_tx]=max(pwr_spec);
    %         if indx_tx==0
    %             % Peak in Power Spectrum not found.
    %             indx_tx=16;  % Index at 1.5 MHz
    %         end
    %     else
    %         % Checkbox is not checked indicating Power is disabled to the Sensor Probe
    %         indx_tx=16;  % Index at 1.5 MHz
    %     end
    
    % *** Important ***
    % indx_rx = indx_tx = 16
    % indx=round(num_fft_smpls/sample_rate*freq_tx)+1;  % *** No Fourier Transform required ***
    % indx=round(1024/100M*1.5M)+1 = round(15.36)+1 = 15+1 = 16
    % Note: indx = indx_tx = 16  at Peak in Power Spectrum
    % *** end ***
    
    indx_tx=16;
    handles.indx_at_pwr_spec_pk_prb=indx_tx;  % Single Frequency Index Pointer
    
    % At Peak in Power Spectrum: indx_at_pwr_spec_pk
    % Magnitude: |Y| (rms volts)
    mag_ft_scalar=abs(fft_daqData(indx_tx));  % mag_prb_scalar=sqrt(pwr_spec_prb(indx_prb))
    phz_ft_scalar=atan2(imag(fft_daqData(indx_tx)),real(fft_daqData(indx_tx)));  % Phase (rad)
    
    % Display Magnitude and Phase in GUI
    % Placed in order of appearance in GUI display
    set(handles.phzPrb_text,'String',num2str(phz_ft_scalar));  % phz_prb (rad)
    set(handles.magPrb_text,'String',num2str(mag_ft_scalar*1000));  % mag_prb (mVolts)
    
    figure;  % Figure #1
    % Probe (prb)
    % Vertical Lines
    x_plt(1)=freq(1)/1000; y_plt(1)=0;
    x_plt(2)=freq(1)/1000; y_plt(2)=pwr_spec(1);
    plot(x_plt,y_plt,'Color','k','LineStyle','-','LineWidth',0.5);  % LineWidth=0.5 (default)
    hold on;
    plot(x_plt(2),y_plt(2),'Color','k','Marker','.','MarkerSize',8);  % MarkerSize=6 (default)
    for i=2:num_nyq_smpls
        x_plt(1)=freq(i)/1000; y_plt(1)=0;
        x_plt(2)=freq(i)/1000; y_plt(2)=pwr_spec(i);
        plot(x_plt,y_plt,'Color','k','LineStyle','-','LineWidth',0.5);  % LineWidth=0.5 (default)
        plot(x_plt(2),y_plt(2),'Color','k','Marker','.','MarkerSize',8);  % MarkerSize=6 (default)
    end
    % Peak in power spectrum
    x_plt(1)=freq(indx_tx)/1000; y_plt(1)=0;
    x_plt(2)=freq(indx_tx)/1000; y_plt(2)=pwr_spec(indx_tx);
    plot(x_plt,y_plt,'Color','r','LineStyle','-','LineWidth',1.0);  % LineWidth=0.5 (default)
    plot(x_plt(2),y_plt(2),'Color','r','Marker','.','MarkerSize',10);  % MarkerSize=6 (default)
    hold off;
    title('\bfCalibration Power Spectrum for "Sensor Probe"','FontSize',12);  % FontSize=12 (default)
    xlabel('Frequency (kHz)','FontSize',12);  % FontSize=12 (default)
    ylabel('|Y(f)|^2','FontSize',12);  % FontSize=12 (default)
    txt_str(1)={['Signal Generator ' num2str(Freq,'%6.3f') ' MHz']};
    txt_str(2)={['Index ' num2str(indx_tx,'%3.0f')]};
    text(10,1e-6,txt_str,'FontSize',12);  % FontSize=12 (default)
    xlim([0 freq_nyq/1000]);
elseif handles.data_processing_method==1
    % In-Phase Quadrature (IQ) Method
    mag_array=abs(raw_data_prb_array);
    phz_array=atan2(imag(raw_data_prb_array),real(raw_data_prb_array));
    mag_iq_scalar=mean(mag_array);
    phz_iq_scalar=mean(phz_array);
    
    % Display Magnitude and Phase in GUI
    % Placed in order of appearance in GUI display
    set(handles.phzPrb_text,'String',num2str(phz_iq_scalar));  % phz_prb (rad)
    set(handles.magPrb_text,'String',num2str(mag_iq_scalar*1000));  % mag_prb (mVolts)
    % No FFT or Power Spectrum
end
% Output: handles.indx_at_pwr_spec_pk_prb
guidata(hObject,handles);  % Update GUI handles structure
% *** end ***
