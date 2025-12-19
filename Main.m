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
        NumProcSymbs = 499;

    % Общее число пилотов в одном символе
        TotNumPilots = 701;

% Вычисляемые параметры
    % Позиции непрерывных пилотов и TPS-сигналов
        [ TPSInds, ContPilotsInds ] = GetPoses();
    % ПСП для модуляции пилотных поднесущих
        PRBS = GenPRBS( OFDM.Nscs )';

% Загрузка сигнала из файла
    NumOfShiftedSamples = 0;
    NumOfNeededSamples = (8192+2048)*NumProcSymbs;
    [Signal, ~] = ReadSignalFromFile( ...
        File, NumOfShiftedSamples, NumOfNeededSamples);

% Инверсия спектра сигнала для корректной дальнейшей обработки
    Signal = conj( Signal );

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
% Определение расположения распределённых пилотов первого OFDM-символа
    % Определение начала и отстройки символа
        [ dt1, df1 ] = OFDM_Syncronization( FSignal, OFDM, ...
            Symbol_Offset, ContPilotsInds, false );

    % Извлечение поднесущих
        UP1 = FSignal( ( 0:OFDM.Nfft-1 ) + dt1 + OFDM.CPLen );
        UP1df = UP1 .* ...
            exp( -1j*2*pi*df1 * (0:OFDM.Nfft-1) / OFDM.Fs );
        FFTVals1 = fftshift( fft( UP1df ) );
        SCs1 = FFTVals1( OFDM.SCsInds );

    % Построение КФ 
        ScatPilotsCorrVal = zeros( 1, 4 );
        lVals = 0 : 3;
        for i = 1 : length( lVals )
            [ ~, ScatPilotPoses ] = GetAllPilotPoses( lVals( i ) );
            ScatPilots = SCs1( ScatPilotPoses +1 );

            RefScatPilots = 4/3 * 2 * ( 1/2 - PRBS( ScatPilotPoses +1 ) );

            ScatPilotsCorrVal( i ) = ScatPilots * conj( RefScatPilots );
        end

    % Определение значения mod( n, 4 ), где n - номер первого ОФДМ-символа
        [ ~, Ind] = max( abs(ScatPilotsCorrVal ) );
        l = lVals( Ind );


% Точная временная и частотная синхронизация для каждого символа
    Coarse_Offsets = zeros( 1, NumProcSymbs );
    Presize_Offsets = zeros( 1, NumProcSymbs );
    Freq_Offsets = zeros( 1, NumProcSymbs );

    Coarse_Offsets ( 1 ) = Symbol_Offset;
    Presize_Offsets( 1 ) = dt1;
    Freq_Offsets   ( 1 ) = df1;
    clear df1 dt1 UP1df UP1 SCs1 FFTVals1;
    % Цикл по символам
        for symIdx = 2:NumProcSymbs
            Coarse_Offsets(symIdx) = Presize_Offsets( symIdx - 1 ) + ...
                OFDM.CPLen + OFDM.Nfft;

            [ Presize_Offsets( symIdx ), Freq_Offsets( symIdx ) ] = ...
                OFDM_Syncronization( FSignal, OFDM, ...
                    Coarse_Offsets( symIdx ), ...
                    GetAllPilotPoses( l + symIdx-1 ), false ...
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
        repmat( ...
            ( 0:size(UsefulParts, 1)-1 )' / OFDM.Fs, 1, NumProcSymbs ...
        ) ...
    );
    
    UsefulParts = UsefulParts .* expVals;

% Преобразование Фурье
    FFTVals = fftshift( fft( UsefulParts, [], 1 ) );

% Выбор отсчётов, на которых находятся активные поднесущие
    SCs = FFTVals( OFDM.SCsInds, : );

% Оценка канала
    % Извлечение пилотов
        RxPilots = SCs( ContPilotsInds +1, :);
    % Генерация опорных пилотов
        RefContPilots = 4/3 * 2 * ( 1/2 - PRBS( ContPilotsInds +1 ) );
    % Коэффициенты передачи канала
        ChEst = RxPilots ./ repmat( RefContPilots, 1, NumProcSymbs );
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
        TPSBits = round( mean( TPSBits, 1 ) );

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
    title( 'После эквалайзинга' );
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
        conv( 1 - 2*TPSBits, ...
            fliplr( 1 - 2*[0 0 1 1 0 1 0 1 1 1 1 0 1 1 1 0] ), "valid" ...
        ) ...
    ); grid on;
    xlabel('Сдвиг до начала синхрослова')
    ylabel('Корреляция')
