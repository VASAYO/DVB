function [CPLen, CorrFun] = Cycle_Prefix_Length_Determination(Signal)
% Функция выполняет определение длины циклического префикса сигнала,
% соответствующего стандарту DVB-T
%
% Входные аргументы:
%   Signal - массив отсчётов сигнала, записанного на частоте дискретизации
%            64/7 * 10^6 Гц.
%   
% Выходные аргументы:
%   CPLen   - длина циклического префикса в отсчётах Signal;
%   CorrFun - массив отсчётов нормированной по Коши-Буняковского 
%             автокорреляционной функции циклического префикса.

% Параметры
    % Длина ОФДМ символа отсчётах
        SymLen = 8192;
    % Число периодов КФ
        NumCorPers = 8;
    % Порог при определении длины ЦП
        CPTreshold = NumCorPers-1;

% cell-массив для результатов корреляций для различных вариантов длин ЦП
    NCorsCell = cell(1, 4);
% cell-массив для результатов корреляций после накопления
    NCorsAccumulate = cell(1, 4);

% Массив различных длин ЦП
    CPLenVals = SymLen * [1/4 1/8 1/16 1/32];

% Цикл по различным длинам ЦП
for CPLenIdx = 1:length(CPLenVals)
    % Расчёт параметров на текущей итерации на основе длины ЦП
        CPLen = CPLenVals(CPLenIdx);
        CorPer = SymLen + CPLen;
        CorLen = CorPer * NumCorPers;

    % Выбор отсчётов сигнала для корреляции
        ShortSig = Signal( 1:(CorPer + CorLen - 1) );

    % Корреляция
    Cors = conv( ...
        ShortSig( 1:CPLen+CorLen-1 ) .* ...
        conj( ShortSig( (1:CPLen+CorLen-1)+SymLen ) ), ...
        ones( 1, CPLen ), "valid" ...
    );

    % Нормировка корреляции по Коши-Буняковского
        En  = conv( ...
            ShortSig .* conj(ShortSig), ones(1, CPLen), "valid" ...
        );
        En1 = En(  1:CorLen ); 
        En2 = En( (1:CorLen) + SymLen ); 
        NCors = Cors ./ sqrt( En1 .* En2 );

    % Сохранение результата для данной итерации в cell-массив
        NCorsCell{CPLenIdx} = NCors;

    % Накопление и сохранение результата
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

CorrFun = NCorsAccumulate{CPLenIdx};
