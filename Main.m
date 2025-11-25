clc;
clear;
close all;

addpath("Lib\");
addpath("Signals\");

% Параметры обработки
    % Структура с параметрами сигнала
        File.HeadLenInBytes = 0;
        File.NumOfChannels = 1;
        File.ChanNum = 0;
        File.DataType = 'int16';
        File.dF = 0;
        % Запись: 'Rus1_small' | 'Rus2_small' | 'Fin_small'
            File.Name = 'Rus2_small';
        % Частота дискретизации
            File.Fs0 = 64/7*10^6;
        % Коэффициенты передискретизации
            File.FsDown = 1;
            File.FsUp = 1;

    % Длительность OFDM символа в отсчётах
        SamplesPerSymbol = 8192;

% Вычисляемые параметры
    % Длина циклического префикса
        if strcmp(File.Name, 'Rus1_small') || ...
                strcmp(File.Name, 'Rus2_small')
            CPLen = 256;

        elseif strcmp(File.Name, 'Fin_small')
            CPLen = 1024;
        end

% Загрузка сигнала из файла
    NumOfShiftedSamples = 0;
    NumOfNeededSamples = (8192+2048)*10;
    [Signal, ~] = ReadSignalFromFile( ...
        File, NumOfShiftedSamples, NumOfNeededSamples);


% Цифровая фильтрация сигнала
    load("LPF_Rus1.mat", "LPF_Rus1");
    FSignal = conv(Signal, LPF_Rus1);

% Символьная синхронизация по циклическому префиксу
    % Число используемых для синхронизации символов
        NumReps = 10;
    % Длина результата корреляции в отсчётах
        CorrLen = NumReps*(8192+CPLen);
    % Память под результат
        CorrRes = zeros(1, CorrLen);

    for k = 1:CorrLen
        CorrRes(k) = FSignal((1:CPLen) + (k-1)) * ...
            conj(FSignal((1:CPLen) + SamplesPerSymbol + (k-1)).');
    end
    % Накопление результата
        Buf = reshape(CorrRes, SamplesPerSymbol+CPLen, []);
        NonKohAcc = sum(abs(Buf), 2);
    % Сдвиг в FSignal до начала циклического префикса
        [~, Ind] = max(NonKohAcc);
    
%% Рисунки
    f1 = figure;
    SPDEstPlotFun(Signal, File.Fs0, 5e3); hold on; grid on;
    SPDEstPlotFun(FSignal, File.Fs0, 5e3); hold on;
    legend('До фильтрации', 'После фильтрации');

    f2 = figure;
    plot(abs(CorrRes));  grid on;

    f3 = figure;
    plot(NonKohAcc);  grid on;