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
            File.Name = 'Rus1_small';
        % Частота дискретизации
            File.Fs0 = 64/7*10^6;
        % Коэффициенты передискретизации
            File.FsDown = 1;
            File.FsUp = 1;
    % Длина ОФДМ символа в отсчётах
        SymLen = 8192;

% Загрузка сигнала из файла
    NumOfShiftedSamples = 0;
    NumOfNeededSamples = (8192+2048)*300;
    [Signal, ~] = ReadSignalFromFile( ...
        File, NumOfShiftedSamples, NumOfNeededSamples);

% Инверсия спектра
    Signal = conj(Signal);

% Цифровая фильтрация сигнала
    load("LPF_Rus1.mat", "LPF_Rus1");
    FSignal = conv(Signal, LPF_Rus1);

%% Символьная синхронизация по циклическому префиксу
% Определение длины ЦП и построение КФ
    [CPLen, CorrFun] = Cycle_Prefix_Length_Determination(FSignal);

% Оценка сдвига до начала ЦП
    [~, Symbol_Offset] = max(CorrFun);

%% ДЗ от 03.12.2025:
% - Задаёмся различными сдвигами относительно грубой оценки начала ОФДМ 
%   символа ( -40:1:40, 81 значение );
% - Задаёмся возможными значениями кратного сдвига частоты 
%   ( (-3:3)*1/T, 7 значений );
% - Обрабатываем один ОФДМ-символ;
% - Двумерная корреляционная функция:
%   * По первому измерению кратный сдвиг частоты;
%   * По второму измерению смещение до соответствующего луча;
%   * По аппликате значение корреляционной функции;
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
% Индексы пилотных поднесущих
    load("ContPilotsInds.mat", "ContPilotsInds");
    [ ~, ContPilotsInds ] = GetPoses();

% Генерация PRS
    PRBS = GenPRBS(6817);
    PRS  = 4/3 * 2 * ( 1/2 - PRBS );

% Массив сдвигов по времени
    tOffsets = Symbol_Offset + ( -40:40 );
% Массив кратных частотных отстроек
% в единицах расстояния между поднесущими
    fOffsets = ( -3:3 );

% Инициализация матрицы со значениями двумерной КФ
    CorrVals = zeros( length( fOffsets ), length( tOffsets ) );

% Цикл по временным сдвигам
for tIdx = 1:length( tOffsets )
    % Выбор отсчётов ЦП и полезной части символа
        Buf = FSignal( (1:CPLen+SymLen)-1 + tOffsets(tIdx) );
        Prefix = Buf( 1:CPLen );
        UsefulPart = Buf( CPLen+1:end );

    % Определение и компенсация дробной частотной отстройки
        dphi = angle( sum( UsefulPart( end-CPLen+1:end ) ) / sum(Prefix) );
        fFrac = dphi / ( 2*pi * SymLen/File.Fs0 );
        UPFrecComp1 = UsefulPart .* ...
            exp( -1j*2*pi*fFrac * ( 0:SymLen-1 )/File.Fs0 );

    % Цикл по грубым частотным сдвигам
    for fIdx = 1:length( fOffsets )
        % Взятие БПФ и грубая частотная подстройка
            FFTVals = fftshift(fft(UPFrecComp1));
            FrecComp2 = circshift(FFTVals, fOffsets(fIdx));

        % Выбор отсчётов, соответствующих поднесущим символа
            SCs = FrecComp2( (689:7505) );

        % Выбор пилотных поднесущих
            Pilots = SCs( ContPilotsInds +1 );

        % Умножение пилотных поднесущих на элементы PRS и сложение
            Pilots = Pilots .* PRS( ContPilotsInds +1 );

        % Вычисление метрики
            CorrVals( fIdx, tIdx ) = abs( sum( Pilots ) );
    end
end

% Построение дифференциального сигнального созвездия для TPS поднесущей для
% 300 ОФДМ символов

% Построение двумерной КФ
    bar3(CorrVals)
