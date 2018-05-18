clear; close all; clc;
%% ======================== IMPORTANT INFO ============================== %
% This code downloads experimental data from the CAIDA servers to determine
% the packet size distribution for IPv4 and IPv6 encapsulation. The data
% serves for a better characterization to generate the input traffic on the
% main.m script.
% 
% The code runs wget to download the data from the server. At the beginning
% of the script, it ensures that wget is installed. If not, it installs it
% using brew. Afterwards, it uses it to download every txt file in the
% desired URL.
%
% If the DATA has already been  stored, then the code does not wget the
% files and jumps straight onto the execution.
%% ======================== CAIDA DB FORMAT ============================= %
% CAIDA Data Base format information (by column number):
% 1. SIZE: is packet size in bytes
% 2. =IPv4=: Total number of IPv4 packets for each packet size (3. SCTP, 4.
% IPv6, 5. ESP, 6. UDP, 7. GRE, 8. ICMP, 9. TCP, 10. UNKNOWN)
% 11. =IPv6=: Total number of IPv6 packets for each packet size ( 12.
% ICMP6, 13. UDP, 14. TCP)
% 15. =IPv6t=: Total number of tunneled IPv6 packets for each packet size
% (16. ICMP6, 17. UDP, 18. TCP
%% ====================================================================== %

%% DOWNLOAD DATA FROM CAIDA DATA BASE (DB)
% Check if we have already downloaded the data from CAIDA
list = dir('www.caida.org/data/passive/trace_stats');
% Download DATA if we haven't done so earier
if isempty(list)
    fprintf('********************************\n');
    fprintf('** Setting the correct environment path in MATLAB...\n');
    % Check if we already have wget
    [~,p] = unix('which wget');
    if isempty(p)
        % Modify path temporarily for current MATLAB session
        path1 = getenv('PATH');
        path1 = [path1 ':/usr/local/bin'];
        setenv('PATH', path1);
        fprintf('** New path configured -> OK\n');
        % Install brew
        fprintf('** Installing required brew to run wget...\n');
        unix('ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"');
        fprintf('** Brew installed correctly -> OK\n');
        % Install wget using brew
        fprintf('** Installing required wget...\n');
        [status,cmdout] = unix('brew install wget --with-libressl');
        if status==0
            fprintf('** Wget installed correctly -> OK\n');
        else
            fprintf('** Error while installing wget. Please check manually.\n');
            return;
        end
    else
        path = strrep(p,'wget','');
        path1 = getenv('PATH');
        path1 = [path1 ':' p];
        setenv('PATH', path1);
    end
    % Download CAIDA Database using wget
    fprintf('** The Following process will download Data from the CAIDA Server\n');
    fprintf('** The process may take a while, please be patient.\n');
    fprintf('** Downloading files from Server now...\n');
    [status,cmdout] = unix('wget --accept txt --mirror --page-requisites --adjust-extension --convert-links --backup-converted --no-parent http://www.caida.org/data/passive/trace_stats/');
    fprintf('** Wget installed correctly -> OK\n');
end

%% ANALYZE DATA
% Variables
minPktSize = 21;
maxPktSize = 1500;
varSpace = (minPktSize:1:maxPktSize);
% Choose Data Set from CAIDA Servers
location = 'chicago-A/';
year = '2014/';
rootFolder = strcat('www.caida.org/data/passive/trace_stats/',location,year);
list = dir(rootFolder);
filesInDir = list(~([list.isdir]));
% Global variables to store results in execution
leg = cell(length(filesInDir),1);  % Store the legend for plotting
resultIpv4BIG = cell(length(varSpace),1);
resultIpv6BIG = cell(length(varSpace),1);
for t = 1:length(varSpace); resultIpv4BIG{t} = []; resultIpv6BIG{t} = []; end;
for idx = 1:length(filesInDir)
    fullPath = strcat(rootFolder,filesInDir(idx).name);
    fileID = fopen(fullPath,'r');
    Intro = textscan(fileID,'%s',1e10,'Delimiter','\n');
    InputLines = Intro{1};
    firstLine = 23;  % Need to start reading from line 26. 
                     % The File first contain headers
    pktSizes = zeros(length(InputLines),1);
	resultIpv4TOT1 = zeros(length(InputLines),1);
    resultIpv6TOT1 = zeros(length(InputLines),1);
    for line = firstLine:length(InputLines)
        myString = strsplit(InputLines{line},' ');
        temp = str2double(myString{1});
        if (temp<=1500) && (temp>=21)
            % Copy x variable
            pktSizes(line) = temp;
            % Results for Ipv4
            resultIpv4TOT1(line) = str2double(myString{2});
            % Results for Ipv6
            resultIpv6TOT1(line) = str2double(myString{11});
            % Store results in final variable
            resultIpv4BIG{varSpace==temp} = [resultIpv4BIG{varSpace==temp} str2double(myString{2})];
            resultIpv6BIG{varSpace==temp} = [resultIpv6BIG{varSpace==temp} str2double(myString{11})];
        end
    end
    % Get CDFs
    totCumSum = cumsum(resultIpv4TOT1);
    resultIpv4CDF = totCumSum./totCumSum(end);
    totCumSum = cumsum(resultIpv6TOT1);
    resultIpv6CDF = totCumSum./totCumSum(end);
    % Concatenate Results for averaging
    % TODO
    %% Plot Results
    % Parse name
    myName = strsplit(filesInDir(idx).name,'.');
    myName = strsplit(myName{3},'-');
    myName = myName{1};
    % Individual set analysis
    figure; subplot(1,2,1); hold on;
    semilogy(pktSizes,resultIpv4TOT1,'lineStyle','none','MarkerSize',4,'Marker','.','MarkerEdgeColor','r');
    semilogy(pktSizes,resultIpv6TOT1,'lineStyle','none','MarkerSize',4,'Marker','.','MarkerEdgeColor','b');
    set(gca,'yscale','log')
    xlim([0 1500]);
    xlabel('Packet Size (bytes)','FontSize',12);
    tit = strcat('Data from Date:',{' '},myName);
    title(tit,'FontSize',14);
    lg = legend('IPV4','IPv6');
    set(lg,'FontSize',12,'Location','NorthWest');
    grid minor;
    subplot(1,2,2); hold on;
    plot(pktSizes,resultIpv4CDF,'lineStyle','-','Color','r','LineWidth',3);
    plot(pktSizes,resultIpv6CDF,'lineStyle','-','Color','b','LineWidth',3);
    xlabel('Packet Size (bytes)','FontSize',12);
    tit = strcat('CDF from Date:',{' '},myName);
    title(tit,'FontSize',14);
    lg = legend('IPV4','IPv6');
    set(lg,'FontSize',12,'Location','NorthWest');
    grid minor;
end
for line = 1:length(varSpace)
    resultIpv4BIG{line} = mean(resultIpv4BIG{line});
    resultIpv6BIG{line} = mean(resultIpv6BIG{line});
end
resultIpv4BIG_2 = cell2mat(resultIpv4BIG);
resultIpv6BIG_2 = cell2mat(resultIpv6BIG);
totCumSum = cumsum(resultIpv4BIG_2);
resultIpv4CDF = totCumSum./totCumSum(end);
totCumSum = cumsum(resultIpv6BIG_2);
resultIpv6CDF = totCumSum./totCumSum(end);
%% Plot Final Result
figure;
subplot(1,2,1); hold on;
title('Average Packet Size','FontSize',14);
semilogy(varSpace,resultIpv4BIG_2,'lineStyle','none','MarkerSize',4,'Marker','.','MarkerEdgeColor','r');
semilogy(varSpace,resultIpv6BIG_2,'lineStyle','none','MarkerSize',4,'Marker','.','MarkerEdgeColor','b');
xlabel('Packet Size (bytes)','FontSize',12);
[lg,icons,plots,legend_text] = legend('IPV4','IPv6');
set(lg,'FontSize',12,'Location','NorthWest');
icons(1).FontSize = 12;
icons(2).FontSize = 12;
icons(4).MarkerSize = 10;
icons(6).MarkerSize = 10;
set(gca,'yscale','log');
xlim([0 1500]);
grid minor;
subplot(1,2,2); hold on;
title('Average CDF','FontSize',14);
plot(varSpace,resultIpv4CDF,'lineStyle','-','Color','r','LineWidth',3);
plot(varSpace,resultIpv6CDF,'lineStyle','-','Color','b','LineWidth',3);
xlabel('Packet Size (bytes)','FontSize',12);
lg = legend('IPV4','IPv6');
set(lg,'FontSize',12,'Location','NorthWest');
grid minor;