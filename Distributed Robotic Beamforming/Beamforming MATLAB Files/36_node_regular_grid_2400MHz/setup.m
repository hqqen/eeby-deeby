clc
clear all

%%
fid=fopen('Sensor_1_2400MHz_Tx_105_195_10.txt');
cdata=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','HeaderLines', 11 );
amat = cdata{8};
pmat = cdata{9};
save('amat_105_195.mat','amat')
save('pmat_105_195.mat','pmat')

fid=fopen('Sensor_1_2400MHz_Tx_195_200_10.txt');
cdata=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','HeaderLines', 11 );
amat = cdata{8};
pmat = cdata{9};
save('amat_195_200.mat','amat')
save('pmat_195_200.mat','pmat')

fid=fopen('Sensor_1_2400MHz_Tx_196_285_10.txt');
cdata=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','HeaderLines', 11 );
amat = cdata{8};
pmat = cdata{9};
save('amat_196_285.mat','amat')
save('pmat_196_285.mat','pmat')

fid=fopen('Sensor_1_2400MHz_Tx_285_195_10.txt');
cdata=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','HeaderLines', 11 );
amat = cdata{8};
pmat = cdata{9};
save('amat_285_195.mat','amat')
save('pmat_285_195.mat','pmat')

fid=fopen('Sensor_1_2400MHz_Tx_285_285_10.txt');
cdata=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','HeaderLines', 11 );
amat = cdata{8};
pmat = cdata{9};
save('amat_285_285.mat','amat')
save('pmat_285_285.mat','pmat')

% fid=fopen('Sensor_1_2400MHz_Tx_15_465_10.txt');
% cdata=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','HeaderLines', 11 );
% amat = cdata{8};
% pmat = cdata{9};
% save('amat_15_465.mat','amat')
% save('pmat_15_465.mat','pmat')