function l = DetermineL( Signal, Symbol_Offset, Freq_Offset, OFDM, ...
    Flag_Draw )
% Функция выполняет определение остатка от деления номера ОФДМ-символа в
% кадре (0, ..., 67) на 4 путём анализа расположения распределённых пилотов
% символа
%
% Входные аргументы:
%   Signal        - массив-строка с отсчётами сигнала;
%   Symbol_Offset - сдвиг до начала символа в Signal;
%   Freq_Offset   - частотная отстройка символа;
%   OFDM          - структура параметров ОФДМ-символа, определённая в 
%                   Main.m;
%   Flag_Draw     - флаг, указывающий на необходимость прорисовки
%                   результатов;
%
% Выходные аргументы:
%   l - остаток от деления номера ОФДМ-символа в кадре на 4;

% Генерация ПСП
    PRBS = GenPRBS( OFDM.Nscs )';

% Взятие отсчётов полезной части символа
    UsefulPart = Signal( (0:OFDM.Nfft-1) + Symbol_Offset + OFDM.CPLen );

% Компенсация частотной отстройки
    FCompUP = UsefulPart .* ...
        exp( -1j*2*pi*Freq_Offset * (0:OFDM.Nfft-1) / OFDM.Fs );

% Извлечение поднесущих
    FFTVals = fftshift( fft( FCompUP ) );
    SCs     = FFTVals( OFDM.SCsInds );

% Построение КФ 
    ScatPilotsCorrVals = zeros( 1, 4 );
    lVals = 0 : 3;
    for i = 1 : length( lVals )
        [ ~, ScatPilotPoses ] = GetAllPilotPoses( lVals( i ) );
        ScatPilots = SCs( ScatPilotPoses +1 );

        RefScatPilots = 4/3 * 2 * ( 1/2 - PRBS( ScatPilotPoses +1 ) );

        ScatPilotsCorrVals( i ) = ScatPilots * conj( RefScatPilots );
    end

% Определение выходного аргумента
    [ ~, Ind] = max( abs(ScatPilotsCorrVals ) );
    l = lVals( Ind );

if exist( "Flag_Draw", "var" )
    if Flag_Draw
        stem( lVals, abs( ScatPilotsCorrVals ) ); grid on;
        xlabel( 'l' );
        ylabel( 'Корреляция' );
    end
end
