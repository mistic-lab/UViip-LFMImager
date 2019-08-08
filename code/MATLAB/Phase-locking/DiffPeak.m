function [FirstTwoPeaks,DiffPeakSample]=DiffPeak(Rx_Signal)

fprintf('DIFFPEAK: ');

Rx_Signal=abs(Rx_Signal);

% Difference between delayed and advance recieved signal

diff1=(Rx_Signal)-[(Rx_Signal(2:end)); 0];

diff2=(Rx_Signal)-[0; (Rx_Signal(1:end-1))];


% 1 means local peak

Binary=and(sign(diff1)>0,sign(diff2)>0);

% Determine the local maximum 

max1=find(Binary>0);

% Sorting for easy pick the max
ID=sort(Rx_Signal(max1),'descend');

% Result

DiffPeakSample=find(Rx_Signal==ID(2))-find(Rx_Signal==ID(1));
FirstTwoPeaks=ID(1:2);
fprintf('\n--> Difference between two peaks = %d Samples\n--> First peak value=%d\n--> Second peak value=%d\n',...
DiffPeakSample,ID(1),ID(2))
fprintf('done\n');
end