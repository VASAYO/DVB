function [tOffset, fOffset] = OFDM_Syncronization( ...
    Signal, OFDM, Symbol_Offset, Flag_isDraw)
% Функция выполняет точную частотную и временную синхронизацию ОФДМ
% символа стандарта DVB-T используя непрерывные пилотные поднесущие
% 
% Входные аргументы:
%   Signal        - отсчёты сигнала, взятые с частотой дискретизации 
%                   64/7 МГц;
%   OFDM          - структура с параметрами ОФДМ-символа, содержащая 
%                   следующие поля:
%       Fs         - частота дискретизации;
%       Nfft       - длина полезной части в отсчётах;
%       Tu         - длительность полезной части в секундах;
%       SCsSpacing - разнос между поднесущими в Гц;
%       CPLen      - длина циклического префикса в отсчётах;
%       Nscs       - число активных поднесущих символа;
%       N0         - число нулевых поднесущих;
%       SCsInds    - индексы отсчётов преобразования Фурье от полезной
%                    части символа, на которых расположены активные
%                    поднесущие.
%   Symbol_Offset - грубая оценка значения сдвига до начала ОФДМ-символа.
% 
% Выходные аргументы:
%   tOffset     - точное значение сдвига до начала ОФДМ-символа, пришедшего 
%                 по главному лучу;
%   fOffset     - точное значение частотной отстройки ОФДМ-символа, 
%                 пришедшего по главному лучу;
%   Flag_isDraw - флаг, указывающий на необходимость прорисовки
%                 результатов.

    % Параметры
        % Частота дискретизации
            Fs         = OFDM.Fs;
        % Длина циклического префикса
            CPLen      = OFDM.CPLen;
        % Длина полезной части
            Nfft       = OFDM.Nfft;
        % Длительность полезной части в секундах
            Tu         = OFDM.Tu;
        % Разнос между поднесущими
            SCsSpacing = OFDM.SCsSpacing;
        % Индексы отсчётов преобразования Фурье от полезной части символа, 
        % на которых расположены активные поднесущие
            SCsInds = OFDM.SCsInds;
    
    % Индексы пилотных поднесущих
        [~, ContPilotPoses] = GetPoses();
    % Генерация PRS
        PRS = GenPRBS( OFDM.Nscs );
        PRS = 4/3 * 2 * ( 1/2 - PRS );
    
    % Значения сдвигов до начала символа
        tOffsetVals = Symbol_Offset + 0%( -40 : 40 );
        tOffsetVals( tOffsetVals < 1 ) = [];
    
    % Значение кратных частотных отстроек в единицах разноса между 
    % поднесущими
        fMultVals = 0%( -3 : 3 );
    
    % Память под значения двумерной корреляционной функции
        CorrVals = zeros( length( fMultVals ), length( tOffsetVals ) );

    % Цикл по сдвигам до начала символа
    for tIdx = 1:length( tOffsetVals )
        tOffset = tOffsetVals( tIdx );

        % Выбор отсчётов символа
            Symbol = Signal( ( 0:Nfft+CPLen-1 ) + tOffset );

        % Определение и компенсация дробной частотной отстройки
            dPhi = angle( ...
                sum( Symbol( end-CPLen+1:end ) ) / ...
                sum( Symbol( 1:CPLen ) ) ...
            );
            % Дробная частотная отстройка в рад/с
                fFrac = dPhi / Tu;

            fFracComp = Symbol .* ...
                exp( -1j*fFrac * (0:Nfft+CPLen-1) / Fs );

            % Цикл по кратным частотным отстройкам
            for fIdx = 1 : length( fMultVals )
                fMult = fMultVals( fIdx );

                % Компенсация кратной частотной отстройки
                    fMultComp = fFracComp .* ...
                        exp( -1j * 2*pi * (fMult*SCsSpacing) * ...
                            (0:Nfft+CPLen-1) / Fs ...
                        );

                % Взятие БПФ
                    FFTVals = fftshift( fft( fMultComp( CPLen+1:end ) ) );

                % Выбор отсчётов, соответствующих активным поднесущим
                    SCs = FFTVals( SCsInds );
                % Выбор пилотных поднесущих
                    Pilots = SCs( ContPilotPoses +1 );
                % Умножение пилотных поднесущих на элементы PRS и сложение
                    Pilots = Pilots .* PRS( ContPilotPoses +1 );

                % Вычисление метрики
                    CorrVals( fIdx, tIdx ) = abs( sum( Pilots ) );
            end
    end

% Прорисовка результатов
    if Flag_isDraw
        bar3(CorrVals);
    end

tOffset = [];
fOffset = [];
