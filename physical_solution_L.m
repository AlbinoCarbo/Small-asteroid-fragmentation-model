% Function for extracting the physical solution for a meteoroid fall from the results of the
% numerical integration during pancake phase. The solution is physically acceptable if elevation >= final elevation 
% and expansion radius <= meteoroid radius.
%
% Input:
% working_path, path of saving model file
% Qfin, final elevation, m
% M0, starting mass, kg
% y, matrix with the solution of the numerical integration:
% y(1) = array speed x, m/s (ECECF)
% y(2) = array speed y, m/s (ECECF)
% y(3) = array speed x, m/s (ECECF)
% y(4) = array position x, m (ECECF)
% y(5) = array position y, m (ECECF)
% y(6) = array position z, m (ECECF)
% y(7) = array mass, kg
% y(8) = expansion radius, m 
% t, array integration time, s
%
% Output:
% N = numero di elementi della soluzione fisica 
% Array quota_meteoroide, m
% Array speed, m/s
% Array lat, long (deg) 
% Array height, m
% Array Dynamic_pressure, MPa
% Array residual fraction mass (adimensional)
% Array diameter, m
%
% Albino Carbognani
% Versione del 26 giugno 2023

function [N, quota_meteoroide, speed, lat, long, height, Dynamic_pressure, residual_mass, diameter]=physical_solution_L(working_path, Qfin, M0, y, t)

wgs84 = wgs84Ellipsoid;
% Polynomial function for the drag coefficient in m^4/kg^2 vs Mach number x (for Mach <= 4).
G_S=@(x)((13.2566820069031e-003)*x^4-(101.740672494020e-003)*x^3+(186.705667397016e-003)*x^2+(104.950251795627e-003)*x+291.373797865474e-003);       

% Read reduced expansion factor for asteroid's radius
wholefile_set = fileread('.\Settings.txt');
set = regexp(wholefile_set,'\$+','split');
nn=str2double(strtrim(set{55}));       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read starting meteoroid's data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

met_dat='meteoroid_data.txt';
Dat1(:,:)=load(met_dat, '-ascii');
Ga0=Dat1(:, 2);         % Drag coefficient
R0=Dat1(:, 5);          % Meteoroid radius, m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the atmospheric profile %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

atmospheric_data='atmospheric_data.txt';
s0=atmospheric_data;
Dat(:,:)=load(s0, '-ascii');

quota_vento=Dat(:, 1);    % Height wind, m 
temperatura=Dat(:, 2);    % Temperature, K
density=Dat(:, 3);        % Air density, kg/m^3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select physical solution, i.e. Strength >= Dynamip pressure and height >= Qfin %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('Select physical solution for pancake phase')
disp('   ')

N=0;

for i=1:length(y(:,4))
     
    % Meteoroid height, m
    [~,~,quota_meteoroide]=ecef2geodetic(wgs84, y(i, 4), y(i, 5), y(i, 6));
       
   if quota_meteoroide >= Qfin && y(i,8)<=nn*R0
       N=N+1;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute dynamic pressure %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('PSL: compute dynamic pressure')
disp('   ')

Dynamic_pressure=zeros(1, N); % MPa

for i=1:N
    % Meteoroid height, m
    [~,~,quota_meteoroide]=ecef2geodetic(wgs84, y(i, 4), y(i, 5), y(i, 6));
    
    % Meteoroid speed respect to soil, m/s
    v=sqrt(y(i, 1)^2+y(i, 2)^2+y(i, 3)^2); 
    
    % Interpolate temperature, K
    T_a=interp1(quota_vento,temperatura,quota_meteoroide);
    
    % Speed sound nearest quota_meteoroide, m/s
    vs=20.0468*sqrt(T_a); 
    
    if v/vs <=4
       Ga=G_S(v/vs);     % Drag for low Mach numbers
    end
    
    if v/vs > 4
       Ga=Ga0;           % Drag for high Mach numbers
    end
    
    % Interpolate air density (kg/m^3)
    rho_a=interp1(quota_vento,density,quota_meteoroide);
    
    % Dynamic pressure, MPa
    Dynamic_pressure(i)=Ga*rho_a*(y(i, 1)^2+y(i, 2)^2+y(i, 3)^2)/10^6;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute meteoroid speed %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('PSL: compute meteoroid speed')
disp('   ')

speed=zeros(1, N); % m/s

for i=1:N
    speed(i)=sqrt(y(i, 1)^2+y(i, 2)^2+y(i, 3)^2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute meteoroid height, lat and long %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('PSL: compute meteoroid height, lat and long (deg)')
disp('   ')

height=zeros(1, N);
lat=zeros(1, N);
long=zeros(1, N);

for i=1:N
    [lat(i), long(i), height(i)]=ecef2geodetic(wgs84, y(i, 4), y(i, 5), y(i, 6));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute residual fraction mass %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
residual_mass=zeros(1, N); 

for i=1:N
    residual_mass(i)=y(i, 7)/M0;
end

%%%%%%%%%%%%%%%%%%%%
% Compute diameter %
%%%%%%%%%%%%%%%%%%%%
diameter=zeros(1, N);

for i=1:N
    diameter(i)=2*(R0+y(i, 8));
end

fid1 = fopen(strcat(working_path, '/Asteroid_model.txt'),'a+');
% Salvataggio dei valori del modello nel file di output
for i=1:N
     fprintf(fid1,'%10.6f \t  %10.6f\t\t %10.4f  \t\t %10.4f \t\t %10.3f \t\t %10.7f \t\t %10.7f \t\t %10.7f\n', height(i)/1000, speed(i)/1000, lat(i), long(i), Dynamic_pressure(i), residual_mass(i), diameter(i), t(i)); 
end

end
