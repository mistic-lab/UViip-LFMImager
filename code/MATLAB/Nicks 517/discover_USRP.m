%% Header
%   USRP initialization file - this script communicates with an Ettus N2xx USRP
%   and is then called by either a Tx script or an Rx script.

%% Discover Radio
connectedRadios = findsdru;
for i = 1:length(connectedRadios)
    if connectedRadios(i).Status == 'Success' %these (1) need attention
        radioFound(i) = true;
        platform(i) = connectedRadios(i).Platform;
        address(i) = connectedRadios(i).IPAddress;
    else
        radioFound(i) = false;
        disp('No radio found');
        address(i) = '192.168.10.2';
        platform(i) = 'N200/N210/USRP2';
    end
    i = i+2;
end

