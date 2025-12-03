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
            File.Name = 'Fin_small';
        % Частота дискретизации
            File.Fs0 = 64/7*10^6;
        % Коэффициенты передискретизации
            File.FsDown = 1;
            File.FsUp = 1;

% Загрузка сигнала из файла
    NumOfShiftedSamples = 0;
    NumOfNeededSamples = (8192+2048)*10;
    [Signal, ~] = ReadSignalFromFile( ...
        File, NumOfShiftedSamples, NumOfNeededSamples);

% Цифровая фильтрация сигнала
    load("LPF_Rus1.mat", "LPF_Rus1");
    FSignal = conv(Signal, LPF_Rus1);

%% Символьная синхронизация по циклическому префиксу
% Параметры
    SymLen = 8192;
    NumCorPers = 8;
    % Порог при определении длины ЦП
        CPTreshold = 6;

% cell-массив для результатов корреляций для различных вариантов длин ЦП
    NCorsCell = cell(1, 4);
    NCorsAccumulate = cell(1, 4);

% Массив различных длин ЦП
    CPLenVals = SymLen * [1/4 1/8 1/16 1/32];

for CPLenIdx = 1:length(CPLenVals) % Цикл по различным длинам ЦП

    CPLen = CPLenVals(CPLenIdx);
    CorPer = SymLen + CPLen;
    CorLen = CorPer * NumCorPers;

    ShortSig = FSignal( 1:(CorPer + CorLen - 1) );
    Cors = conv( ...
        ShortSig( 1:CPLen+CorLen-1 ) .* ...
        conj( ShortSig( (1:CPLen+CorLen-1)+SymLen ) ), ...
        ones( 1, CPLen ), "valid" ...
    );

    En  = conv( ShortSig .* conj(ShortSig), ones(1, CPLen), "valid" );
    En1 = En(  1:CorLen ); 
    En2 = En( (1:CorLen) + SymLen ); 
    NCors = Cors ./ sqrt( En1 .* En2 );

    NCorsCell{CPLenIdx} = NCors;

    % Накопление 
        Buf = reshape(NCors, [], NumCorPers);
        NCorsAccumulate{CPLenIdx} = abs( sum(Buf, 2) );
end

% Определение длины ЦП
    isTresholdExceeded = zeros(1, 4);

    for CPLenIdx = 1:length(CPLenVals)
        isTresholdExceeded(CPLenIdx) = ...
            sum( NCorsAccumulate{CPLenIdx} >= CPTreshold);
    end
    CPLenIdx = find( isTresholdExceeded );
    CPLen    = CPLenVals( CPLenIdx );
    % Грубая оценка временной синхронизации
        Buf = NCorsAccumulate{CPLenIdx};
        [~, Symbol_Offset] = max(Buf);

%% ДЗ от 03.12.2025:
% - Задаёмся различными сдвигами относительно грубой оценки начала ОФДМ 
%   символа (-40:1:40, 81 значение);
% - Задаёмся возможными значениями кратного сдвига частоты 
%   ( (-3:3)*1/(2T), 7 значений );
% - Обрабатываем один ОФДМ-символ;
% - Двумерная корреляционная функция:
%   * По первому измерению кратный сдвиг частоты;
%   * По второму измерению отсчёт, соответствующий предположительному
%     началу сигнала;
%   * По ординате значение корреляционной функции;
%
% - Последовательность действий для каждой точки корреляционной функции:
%   1. Сдвигаемся по времени на определенный отчёт (смещение от -40 до 40
%      относительно грубой оценки синхронизации);
%   2. Определяем и компенсируем дробный сдвиг по частоте;
%   3. Берём ДПФ и имплементируем кратный сдвиг по частоте;
%   4. Извлекаем пилотные поднесущие;
%   5. Коррелируем набор пилотов с эталонной последовательностью.
%   * Результирующая картинка - картина многолучевости в канале, с пиками 
%     в различных сдвигах по времени и частоте.
%
% - bar3 для рисования трёхмерной картинки столбцами.

%% Точная символьная и частотная синхронизация
% Массив сдвигов относительно грубой временной синхронизации
    ShiftSamps = ( -40:40 );
% Массив кратных частотных сдвигов в единицах 1/T
    ShiftFreqs = ( -3:3 );

% Массив значений корреляционной функции
    CorrVals = zeros(length(ShiftSamps), length(ShiftFreqs));

% Цикл по сдвигам во времени
for shIdx = 1:length(ShiftSamps)
    % Текущий сдвиг до начала символа
        TOffset = Symbol_Offset+ShiftSamps(shIdx);

    % Цикл по кратным сдвигам частоты
    for frIdx = 1:length(ShiftFreqs)
        % Текущий кратный частотный сдвиг
            FOffset = ShiftFreqs(frIdx);

        % Сдвигаемся к началу ОФДМ-символа и выбираем его отсчёты 
            CPSymbol = FSignal( ( 1:SymLen+CPLen )-1 + TOffset );

        % Значения корреляции отсчётов ЦП и символа
            % CorrCloud = 

    end
end

% Обязательно доделать к следующей паре!
