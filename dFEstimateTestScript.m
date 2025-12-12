% Тестирование оценки и компенсации частотной отстройки в системах CP-OFDM
    clc; clear; 
    close all;

% Параметры ОФДМ символа
    % Разнос между поднесущими в Гц
        OFDM.SCsSpacing = 20e3;
    % Размер БПФ
        OFDM.Nfft  = 1024;
    % Число нулевых поднесущих по бокам сигнала
        OFDM.No    = 256;
    % Размер ЦП в единицах длины полезной части символа
        OFDM.CPLen = 1/4;
    % Порядок модуляции (ФМ)
        OFDM.ModOrder = 2;

% Параметры канала
    % Дробная частотная отстройка в Гц
        dfFrac = 9000;
    % Кратная частотная отстройка в единицах разноса между поднесущими
        dfMult = 8;

% Вычисляемые параметры
    % Число информационных поднесущих
        OFDM.Ndata = OFDM.Nfft - OFDM.No;
    % Длительность полезной части символа в секундах
        OFDM.Ts = 1 / OFDM.SCsSpacing;
    % Частота дискретизации
        Fs = OFDM.Nfft * OFDM.SCsSpacing;
    % Индексы, соответствующие активным поднесущим
        OFDM.ActiveInds = OFDM.No / 2 + 1 : OFDM.Nfft - OFDM.No / 2;
    % Часто используемое значение
        OFDM.log2M = log2(OFDM.ModOrder);
    % Размер ЦП в отсчётах
        OFDM.CPLenSamples = OFDM.CPLen * OFDM.Nfft;
    % Итоговая частотная отстройка в канале 
        df = dfFrac + dfMult * OFDM.SCsSpacing;

% Генерация ОФДМ символа
    % Генерация данных
        TxData = randi([0 1], 1, OFDM.Ndata * OFDM.log2M);
    % Маппинг
        TxSymbols = pskmod(TxData, OFDM.ModOrder, 0, "gray", ...
            "InputType", "bit");

    % Инициализация массива поднесущих
        FFTVals = zeros(1, OFDM.Nfft);
    % Заполнение активных поднесущих
        FFTVals(OFDM.ActiveInds) = TxSymbols;
    % Обратное преобразование Фурье
        UsefulPart = ifft( ifftshift( FFTVals ) );
    % Циклический префикс
        CP = UsefulPart( end-OFDM.CPLenSamples+1:end );
    % ОФДМ символ
        TxWaveform = [CP, UsefulPart];

% Добавление частотной отстройки
    RxWaveform = TxWaveform .* ...
        exp( 1j*2*pi*df * ( 0:length(TxWaveform)-1 ) / Fs );

%% Приёмник
% Оценка дробной частотной отстройки
    % Оценка набега фазы за длительность полезной части символа
        dPhi = angle( ...
            sum( RxWaveform( end-OFDM.CPLenSamples+1:end ) ) / ...
            sum( RxWaveform( 1:OFDM.CPLenSamples ) ) ...
        );

    % Оценка дробного частотного сдвига
        dfFracEst = dPhi / ( 2*pi * OFDM.Nfft / Fs );

    % Компенсация дробного частотного сдвига
        FreqComp = RxWaveform .* ...
            exp( -1j*2*pi* dfFracEst * ( 0:length(RxWaveform)-1 ) / Fs );

% Компенсация кратной частотной отстройки
    RxSCs = fftshift( fft( FreqComp( OFDM.CPLenSamples+1:end ) ) );
    FreqComp2 = circshift( RxSCs, -dfMult );

% Извлечение данных из информационных поднесущих
    RxSymbols = FreqComp2( OFDM.ActiveInds );

% Демаппинг 
    RxData = pskdemod( RxSymbols, OFDM.ModOrder, 0, "gray", ...
        "OutputType", "bit");

%% Верификация резульатов
    numErr = sum( RxData ~= TxData );

%% Логгирование 
    fprintf( '%s Error Rate: %d / %d\r\n', ...
        datestr(now), numErr, length(TxData) );