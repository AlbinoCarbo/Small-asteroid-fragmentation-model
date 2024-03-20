% Questa funzione, usata in abbinamento con ode45, permette di concludere
% l'integrazione del sistema di equazioni differenziali quando position vale zero.
% Usata per terminare l'integrazione della fase di espansione del pancake.
%
% Albino Carbognani, INAF-OAS
% Versione del 26 giugno 2023

function [position,isterminal,direction] = EventsFunction_L(t,y)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read meteoroid's radius %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    met_dat='meteoroid_data.txt';
    Dat1(:,:)=load(met_dat, '-ascii');
    R0=Dat1(:, 5);  % Meteoroid's original radius, m

    % Read reduced expansion factor for asteroid's radius
    wholefile_set = fileread('.\Settings.txt');
    set = regexp(wholefile_set,'\$+','split');
    nn=str2double(strtrim(set{55}));
    
    % Condition to exit from integration 
    position = y(8)-nn*R0;              % The value that we want to be zero
    isterminal = 1;                     % Halt integration 
    direction = 0;                      % The zero can be approached from either direction
end