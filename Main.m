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
    NumOfNeededSamples = (8192+2048)*10;
    [Signal, ~] = ReadSignalFromFile( ...
        File, NumOfShiftedSamples, NumOfNeededSamples);

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
% Номера поднесущих, на которых располагаются непрерывные пилоты
    load("ContPilotsInds.mat", "ContPilotsInds");

% Сдвиги до моментов, определяющих начало ОФДМ символа
    TOffsets = Symbol_Offset + 0;%( -40:40 );
    TOffsets(TOffsets < 1) = [];

% Массив грубых частотных сдвигов в единицах 1/T
    CoarseFreqShifts = 0;%( -3:3 );

% Память под значения корреляционной функции
    CorrVals = zeros(length(CoarseFreqShifts), length(TOffsets));

% Цикл по сдвигам до лучей
for tIdx = 1:length(TOffsets)
    RayOffset = TOffsets(tIdx);

    % Цикл по кратным сдвигам частоты
    for fIdx = 1:length(CoarseFreqShifts)
        CoarseFreqShift = CoarseFreqShifts(fIdx);

        % Сдвигаемся до луча и выбирает отсчёты ЦП и символа
            SymWithCP = FSignal( (1:CPLen+SymLen)-1 + RayOffset );

        % Определяем и компенсируем дробный частотный сдвиг
            % Оценка набега фазы
            Ksi = angle( ...
                sum( SymWithCP( end-CPLen+1:end ) ) / ...
                    sum( SymWithCP( 1:CPLen ) ) ...
            );

            % Дробный частотный сдвиг
                FracFreqShift = Ksi / ( 2*pi * SymLen/File.Fs0 );
            % Компенсация
                NoFraqShift = SymWithCP .* ...
                    exp( -1j*2*pi*FracFreqShift * ( 0:length(SymWithCP)-1 ) / File.Fs0 );

        % Взятие ДПФ и имплементация кратного сдвига частоты
            SCs = fftshift( fft( NoFraqShift( CPLen+1:end ) ) );
            SCsShift = circshift(SCs, CoarseFreqShift);

        figure; semilogy(abs(SCsShift))
    end
end

% Обязательно доделать к следующей паре!
