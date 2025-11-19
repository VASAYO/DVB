clc;
clear;
% close all;

addpath("Lib\");
addpath("Signals\");

File.Name = 'Rus1_small';
File.HeadLenInBytes = 0;
File.NumOfChannels = 1;
File.ChanNum = 0;
File.DataType = 'int16';
File.Fs0 = 64/7*10^6;
File.dF = 0;
File.FsDown = 1;
File.FsUp = 1;

NumOfShiftedSamples = 0;
NumOfNeededSamples = (8192+2048)*10;
[Signal, ~] = ReadSignalFromFile(File, NumOfShiftedSamples, NumOfNeededSamples);


% Цифровая фильтрация
    load("LPF_Rus1.mat", "LPF_Rus1");
    FSignal = conv(Signal, LPF_Rus1);

% Символьная синхронизация по циклическому префиксу
NumReps = 4;
CorrLen = NumReps*(8192+256);
CPLen = 256;

CorrRes = zeros(1, CorrLen);

for k = 1:CorrLen
    CorrRes(k) = sum(FSignal((1:CPLen) + (k-1)) .* conj(FSignal((1:CPLen) + 8192 + (k-1))));
end

plot(abs(CorrRes)); hold on;