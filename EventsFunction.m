% Questa funzione, usata in abbinamento con ode45, permette di concludere
% l'integrazione del sistema di equazioni differenziali quando la variabile
% position vale zero. Con questa funzione l'integrazione si ferma quando 
% la quota del meteoroide in caduta è zero. 
%
% Albino Carbognani, INAF-OAS
% Versione del 12 giugno 2023

function [position,isterminal,direction] = EventsFunction(t,y)
    
    % Reference ellipsoid for World Geodetic System of 1984
    wgs84 = wgs84Ellipsoid;
    % Meteoroid height, m
    [~,~,quota_meteoroide]=ecef2geodetic(wgs84, y(4), y(5), y(6));
    
    % Condition to exit from integration 
    position = quota_meteoroide;    % The value that we want to be zero
    isterminal = 1;                 % Halt integration 
    direction = 0;                  % The zero can be approached from either direction
end