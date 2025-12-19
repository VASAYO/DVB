function [AllPilotPoses, ScatPilotPoses] = GetAllPilotPoses( n )
% Функция возвращает массив-столбец позиций всех (непрерывных и 
% распределённых) пилотных поднесущих n-го (n = 0, ..., 67) OFDM-символа в 
% кадре.
%
% Выходные аргументы:
%   AllPilotPoses  - массив-столбец позиций всех поднесущих OFDM-символа;
%   ScatPilotPoses - массив-столбец позиций распределённых поднесущих 
%                    OFDM-символа;

% Позиции непрерывных поднесущих
    [~, ContPilotPoses] = GetPoses();

% Позиции распределённых поднесущих
    ScatPilotPoses = 3 * mod( n, 4 ) + 12 * ( 0 : 568 )';
    ScatPilotPoses( ScatPilotPoses > 6816 ) = [];

% Объединение позиций пилотов
    AllPilotPoses = unique( [ ContPilotPoses'; ScatPilotPoses ] );
