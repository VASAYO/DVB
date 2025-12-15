clc;
clear;
close all;

addpath("Lib\");
addpath("Signals\");

% Параметры
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

    % Структура с параметрами ОФДМ-символа, содержащая следующие поля:
    % 
    %   Fs         - частота дискретизации;
    %   Nfft       - длина полезной части в отсчётах;
    %   Tu         - длительность полезной части в секундах;
    %   SCsSpacing - разнос между поднесущими в Гц;
    %   CPLen      - длина циклического префикса в отсчётах;
    %   Nscs       - число активных поднесущих символа;
    %   N0         - число нулевых поднесущих;
    %   SCsInds    - индексы отсчётов преобразования Фурье от полезной
    %                части символа, на которых расположены активные
    %                поднесущие.
        OFDM.Fs    = File.Fs0;
        OFDM.Nfft  = 8192;
        OFDM.Nscs  = 6817;
        OFDM.CPLen = [];

        OFDM.Tu         = OFDM.Nfft / OFDM.Fs;
        OFDM.SCsSpacing = 1 / OFDM.Tu;
        OFDM.N0         = OFDM.Nfft - OFDM.Nscs;
        OFDM.SCsInds    = ( 0 : OFDM.Nscs-1 ) + ( OFDM.N0 + 1 ) / 2 + 1;

    % Число обрабатываемых ОФДМ-символов
        NumProcSymbs = 500;

% Вычисляемые параметры
    % Индексы непрерывных пилотов
        [ TPSInds, ContPilotsInds ] = GetPoses();

% Загрузка сигнала из файла
    NumOfShiftedSamples = 0;
    NumOfNeededSamples = (8192+2048)*NumProcSymbs;
    [Signal, ~] = ReadSignalFromFile( ...
        File, NumOfShiftedSamples, NumOfNeededSamples);

% Инверсия спектра сигнала для корректной дальнейшей обработки
    Signal = conj(Signal);

% Цифровая фильтрация сигнала
    load("LPF_FIR_Coeffs.mat", "LPF_FIR_Coeffs");
    FSignal = conv(Signal, LPF_FIR_Coeffs);

%% Символьная синхронизация по циклическому префиксу
% Определение длины ЦП и построение КФ
    [OFDM.CPLen, CorrFun] = Cycle_Prefix_Length_Determination(FSignal);

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

%% Обработка 300 последовательных ОФДМ-символов
% Грубые оценки сдвигов до начала каждого символа
    Coarse_Offsets = Symbol_Offset + ...
        ( 0:NumProcSymbs-1 ) * ( OFDM.CPLen+OFDM.Nfft );

% Точная временная и частотная синхронизация для каждого символа
    Presize_Offsets = zeros( 1, NumProcSymbs);
    Freq_Offsets = zeros( 1, NumProcSymbs);
    % Цикл по символам
        for symIdx = 1:NumProcSymbs
            [ Presize_Offsets(symIdx), Freq_Offsets(symIdx) ] = ...
                OFDM_Syncronization(FSignal, OFDM, ...
                    Coarse_Offsets( symIdx ), ContPilotsInds, 0 ...
                );
        end

% Полезные части всех символов
    UsefulParts = zeros( OFDM.Nfft, NumProcSymbs );

    for symIdx = 1:NumProcSymbs
        UsefulParts( :, symIdx ) = ...
            FSignal( ( 0:OFDM.Nfft-1 ) + Presize_Offsets( symIdx ) + ...
                OFDM.CPLen ...
            );
    end

% Компенсация частотных отстроек
    expVals = exp( ...
        -1j * 2 * pi * ...
        repmat( Freq_Offsets, size(UsefulParts, 1), 1 ) .* ...
        repmat( ( 0:size(UsefulParts, 1)-1 )' / OFDM.Fs, 1, NumProcSymbs ) ...
    );
    
    UsefulParts = UsefulParts .* expVals;

% Преобразование Фурье
    FFTVals = fftshift( fft( UsefulParts, [], 1 ) );

% Выбор отсчётов, на которых находятся активные поднесущие
    SCs = FFTVals( OFDM.SCsInds, : );

% Оценка канала
    % Извлечение пилотов
        RxContPilots = SCs( ContPilotsInds +1, :);
    % Генерация опорных пилотов
        PRS = GenPRBS( OFDM.Nscs )';
        RefContPilots = 4/3 * 2 * ( 1/2 - PRS( ContPilotsInds +1 ) );
    % Коэффициенты передачи канала
        ChEst = RxContPilots ./ repmat( RefContPilots, 1, NumProcSymbs );
    % Интерполяция на все поднесущие 
        ChEstInterp = interp1( ContPilotsInds+1, ChEst, 1:OFDM.Nscs );

% ZF-эквалайзер
    SCsEq = SCs ./ ChEstInterp;

% Обработка TPS сигнала
    % Извлечение поднесущих
        SCsTPS = SCsEq( TPSInds +1, : );
    % Построение дифференциального созвездия
        SCsTPSDiff = SCsTPS( :, 2:end ) .* SCsTPS( :, 1:end-1 );
    % Демодуляция
        TPSBits = pskdemod( SCsTPSDiff, 2, 0 );

%% Прорисовка результатов
% Сигнальные созвездия до и после эквалайзинга
    figure("Name", "Принятые сигнальные созвездия");
    subplot( 1, 2, 1 );
    plot( SCs( :, 1 ), '.' ); axis equal; grid on;
    title('До эквалайзинга');
    xlabel('I');
    ylabel('Q');

    subplot( 1, 2, 2 );
    plot( SCsEq( :, 1 ), '.' ); axis equal; grid on;
    title( 'После эквалайзинга при' );
    xlabel('I');
    ylabel('Q');

% Созвездие TPS сигнала
    figure("Name", "Созвездия TPS сигнала")
    subplot( 1, 2, 1 );
    plot( SCsTPS( 10, : ), '.' ); axis equal; grid on;
    title('Недифференциальное');
    xlabel('I');
    ylabel('Q');

    subplot( 1, 2, 2 );
    plot( SCsTPSDiff( 10, : ), '.' ); axis equal; grid on;
    title('Дифференциальное');
    xlabel('I');
    ylabel('Q');

% Корреляция бит TPS с синхрословом
    figure(Name="Корреляция TPS с синхрословом");
    plot( ...
        conv( 1 - 2*TPSBits(1, :), ...
            fliplr( 1 - 2*[0 0 1 1 0 1 0 1 1 1 1 0 1 1 1 0] ), "valid" ...
        ) ...
    ); grid on;
    xlabel('Сдвиг до начала синхрослова')
    ylabel('Корреляция')