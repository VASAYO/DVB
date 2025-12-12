clc; clear;

% Параметры 
    % Частота дискретизации
        Fs = 64/7 * 1e6;
    
    % Размер преобразования Фурье
        OFDM.Nfft = 8192;
    % Число активных поднесущих
        OFDM.Nactive = 6817;

% Вычисляемые параметры
    % Длительность полезной части символа в секундах
        OFDM.Ts = OFDM.Nfft / Fs;
    % Разнос поднесущих в Гц
        OFDM.SCsSpacing = 1 / OFDM.Ts;
    % Число нулевых поднесущих
        OFDM.N0 = OFDM.Nfft - OFDM.Nactive;
    % Индексы активных поднесущих
        OFDM.ActiveInds = ( 0:OFDM.Nactive-1 ) + ( OFDM.N0 + 1 ) / 2 + 1
        