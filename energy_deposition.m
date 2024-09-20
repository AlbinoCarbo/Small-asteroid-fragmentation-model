% Script to compute energy deposition in atmosphere of a meteoroid
%
% Input:
% Time, s
% Speed, km/s
% Height, km
% Relative_mass,
% M0, original mass, kg
%
% Output:
% Energy_dep1, kt/km
% Height1, km
% Power, total power lost, Watt
%
% Albino Carbognani, INAF-OAS
% Versione 10 settembre 2024

function [Energy_dep1, Height1, Time1, Power]=energy_deposition(Time, Speed, Height, Relative_mass, M0)

% Kinetic energy, kt
K=(0.5*M0*Relative_mass.*(1000*Speed).^2)/(4.184*10^12); 

% Gradient of kinetic energy, kt
% N.B. Al di sopra di circa 70 km di quota dK>0 per effetto dell'accelerazione di gravità,
% che aumenta la velocità di caduta, mentre il drag con l'atmosfera è trascurabile.
dK=gradient(K);

% Gradient of height, km
dHeight=gradient(Height);

% Compute energy deposition array, kt/km
Energy_dep=dK./dHeight; 
ix = find(Energy_dep>0,1); % Indice del primo elemento di Energy_dep > 0

% Array dell'energy deposition con tutti gli elementi > 0, 
% così si può usare la scala logaritmica per il plot della figura
Energy_dep1=Energy_dep(ix:1:end); 
Height1=Height(ix:1:end);
Time1=Time(ix:1:end);

% Compute fireball's absolute magnitude
dT1=gradient(Time1);
dK1=gradient(K(ix:1:end));
Power=(4.184*10^12)*abs(dK1)./dT1;      % Fireball's total power, Watt

end
