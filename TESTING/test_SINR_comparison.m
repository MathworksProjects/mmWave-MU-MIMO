clear all; close all; clc;

addpath('UTILITIES','-end');  % Add utilities folder at the end of search path
addpath('../data','-end');  % Add data folder at the end of search path

% load('TABLE-SNR-loc4.mat','SINRLCMV','SINRCBF','SINRHEU','nAntennasList','nUsersList');
load('TABLE-SNR-loc4-v2.mat','SINRLCMV','SINRCBF','SINRHEU','nAntennasList','nUsersList');

improve = SINRHEU - SINRLCMV;

figure; hold on;
leg = cell(10,1);
for i = 1:10
    plot(1:length(nAntennasList),improve(i,:),'linewidth',1.5);
    leg(i) = strcat(num2str(i),{' '},'users');
end
nAntennasList = {'12^2','14^2','16^2','18^2','20^2','22^2'};
legend(leg);
xticklabels(nAntennasList);
ylabel('improvement (dB)');
xlabel('Array size');
grid minor;