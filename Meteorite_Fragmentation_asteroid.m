% Meteoroid fragmentation for little asteroids (SPECIAL EDITION WITH STRENGTH AND WITHOUT PANCAKE)
% Matlab script for the computation of an atmospheric meteoroid path in a 3D ECEF system
% with main explosion and fragment fall. 
% Input data with initial conditions are read from the Settings.txt file. It needs 
% the atmospheric profile from the starting height to the ground.
%
% ALGORITHM
% 1-Read settings file
% 2-Copy file "Settings.txt" in fireball's folder
% 3-Compute and save file "meteoroid_data.txt" in main folder to be read from ODEsystem.m, 
%   EventsFunction.m and EventsFunction_L.m
% 4-Define working path and atmospheric profile location
% 5-Open the output model file
% 5-Compute ECEF atmospheric profile from geopotential atmospheric profile 
%   and save ECEF "atmospheric_data.txt" in main folder to be read from ODEsystem.m script
% 6-Copy file "atmospheric_data.txt" also in fireball's folder
% 7-Compute meteoroid starting speed in ECEF reference system
% 8-Compute meteoroid starting position in ECEF reference system
% 9-Define starting conditions vector
% 10-Integrate 3D motion equations until main fragmentation with Runge-Kutta 4th/5th order ODE
%    solver taking into account Coriolis and centrifugal forces (ODEsystem.m script)
% 11-Select physical solution (only points above Earth's surface and with strength superior to dynamic pressure)
% 12-
% 13-
%
% 20-End script
%
% Matlab's functions:
%
% [U,V,W] = ned2ecefv(uNorth,vEast,wDown,lat0,lon0) 
% returns vector components U, V, and W in a geocentric Earth-centered Earth-fixed (ECEF) system 
% corresponding to vector components uNorth, vEast, and wDown in a local north-east-down (NED) system. 
% Specify the origin of the system with the geodetic coordinates lat0 and lon0. 
% Each coordinate input argument must match the others in size or be scalar.
%
% [X,Y,Z] = geodetic2ecef(spheroid,lat,lon,h) transforms the geodetic coordinates specified by lat, lon, and h 
% to the geocentric Earth-Centered Earth-Fixed (ECEF) Cartesian coordinates specified by X, Y, and Z. 
% Specify spheroid as the reference spheroid for the geodetic coordinates.
%
% interpolation.m, function for the linear interpolation of wind speed, air density and temperature 
% between two different altitudes of the atmospheric profile.
%
% Albino Carbognani, INAF-OAS
% Version Sep 19, 2024

clear all

disp('                                                                       ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%               METEORITE FALL MODEL FOR LITTLE ASTEROIDS             %')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%              ATMOSPHERIC PATH AND FRAGMENTS STREWN FIELD            %')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('                                                                       ')
disp('                    by Albino Carbognani (INAF-OAS)                    ')
disp('                              Sep 2024                                 ')
disp('                                                                       ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('   ')

% Reference ellipsoid for World Geodetic System of 1984.
% NOTE: by default, the lengths of the semimajor axis and semiminor axis are in
% meters.
wgs84 = wgs84Ellipsoid;

%%%%%%%%%%%%%%%%%%%%%%
% Read settings file %
%%%%%%%%%%%%%%%%%%%%%%
disp('Read settings file')
disp('   ')

wholefile_set = fileread('.\Asteroid_2023CX1\Settings_2023CX1.txt');
% Split of the Settings and initialization of the string variables
set1 = regexp(wholefile_set,'\$+','split');

fireball_results=strtrim(set1{4});      % Fireball's folder results
working_path=strtrim(set1{7});          % Path of the atmospheric profile
ATMO=strtrim(set1{10});                 % Atmospheric profile name
E0=str2double(strtrim(set1{13}));       % Meteoroid kinetic energy, Mt
rho_M=str2double(strtrim(set1{16}));    % Mean density meteoroid, kg/m^3
V0=str2double(strtrim(set1{19}));       % Starting speed, m/s
I0=str2double(strtrim(set1{22}));       % Trajectory inclination at starting speed, degrees
A0=str2double(strtrim(set1{25}));       % Trajectory azimut at starting speed, degrees
Q0=str2double(strtrim(set1{28}));       % Starting height above Earth surface, m
LAT=str2double(strtrim(set1{31}));      % Latitude starting point, degrees
LONG=str2double(strtrim(set1{34}));     % Longitude starting point, degrees
Qfin=str2double(strtrim(set1{37}));     % Mean elevation of the ground (estimated) on which the meteoroid falls, m
T0=str2double(strtrim(set1{40}));       % Maximum integration time, s
Ga=str2double(strtrim(set1{43}));       % Asintotic drag coefficient
F=str2double(strtrim(set1{46}));        % Dimensionless spherical form factor
Strength=str2double(strtrim(set1{49})); % Mean body strength (MPa)
M1=str2num(strtrim(set1{52}));          % Final fragments mass, kg
f=str2double(strtrim(set1{55}));        % Fragments mass correction factor (to compensate ablation)
fireball_name=strtrim(set1{58});        % Fireball's name
LE=str2double(strtrim(set1{61}));       % Luminous Efficiency (adimensional)

% Fragments starting mass, kg
M1=M1/f;

% Compute meteoroid starting mass, kg
M0=2*E0*(4.1840e+15)/(V0^2);

% Compute meteoroid starting effective radius, m
R0=(3*M0/(4*pi*rho_M))^(1/3);

% Compute starting cross section meteoroid, m^2
S0=F*(M0/rho_M)^(2/3);

%%%%%%%%%%%%%%%%%%%%%%%
% Save meteoroid data %
%%%%%%%%%%%%%%%%%%%%%%%
disp('Save meteoroid data for ODEsystem.m script')
disp('   ')

% Save data in "meteoroid_data.txt", file that will be read by ODEsystem.m,
% and EventsFunction.m 
fid0 = fopen('meteoroid_data.txt','w');
fprintf(fid0,'%10.5f \t %10.5f \t %10.5f \t %10.5f \t %10.5f\n', F, Ga, rho_M, Strength, R0);   
fclose(fid0);

% Copy the current meteoroid file in the current fireball's folder
copyfile("meteoroid_data.txt", fireball_results)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path of the atmospheric data file %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Atmospheric data file')
disp('   ')
 
atmospheric_profile=strcat(working_path, ATMO);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open the output model file %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid1 = fopen(strcat(fireball_results, 'Asteroid_model.txt'),'w'); 
fprintf(fid1, '%%    \n');
fprintf(fid1, strcat('%%', " ", fireball_name," ", 'FRAGMENTATION WITH NO PANCAKE PHASE \n'));
fprintf(fid1, '%%    \n');
fprintf(fid1, '%% Starting Mass (kg)                       %10.2f \n', M0);
fprintf(fid1, '%% Diameter (m)                             %10.2f \n', 2*R0);
fprintf(fid1, '%% Mean density (kg/m^3)                    %10.2f \n', rho_M);
fprintf(fid1, '%% M/A (kg/m^2)                             %10.2f \n', M0/S0);
fprintf(fid1, '%% Asintotic drag coefficient               %10.2f \n', Ga);
fprintf(fid1, '%% Starting speed (km/s)                    %10.2f \n', V0/1000);
fprintf(fid1, '%% Starting inclination (degree)            %10.2f \n', I0);
fprintf(fid1, '%% Starting height (km)                     %10.2f \n', Q0/1000);
fprintf(fid1, '%% Starting strength (MPa)                  %10.2f \n', Strength);
fprintf(fid1, '%% Fragments Mass (kg)                      %10.2f \n', M1(1));
for J=2:length(M1)
fprintf(fid1, '%%                                          %10.3f \n', M1(J));
end
fprintf(fid1, '    \n');
fprintf(fid1, '%% Height (km)     Speed (km/s)       Lat (degrees)     Long (degrees)    Dynamic pressure (MPa)    Relative mass    Diameter (m)   Time (s)\n');
fprintf(fid1, '    \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open the strewn field file %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid5 = fopen(strcat(fireball_results, 'Asteroid_strewn_field.txt'),'w');
fprintf(fid5, strcat('%%', " ", fireball_name," ", 'FRAGMENTATION WITH NO PANCAKE PHASE \n'));
fprintf(fid5, '%% Starting strength (MPa)                  %10.2f \n', Strength);
fprintf(fid5, '    \n');
fprintf(fid5, '%%  Lat (degrees)     Long (degrees)    Mass (kg)\n');
fprintf(fid5, '    \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute and save "atmospheric_data.txt" from ATMO %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
string1=strcat('Compute and save ECEF atmospheric data');
disp(string1)
disp('   ')

s0=atmospheric_profile;
Dat(:,:)=load(s0, '-ascii');
pres=Dat(:, 1)';                            % Pressure, hPa
ugr=Dat(:, 2)';                             % Wind u speed, m/s
vgr=Dat(:, 3)';                             % Wind v speed, m/s
temp=Dat(:, 4)';                            % Temperature, K
quota_vento=Dat(:, 6)';                     % Quote, m

density=0.1*(3.483676*pres./temp);          % Air density, kg/m^3
Vw_north=vgr;                               % Wind speed north component, m/s
Vw_east=ugr;                                % Wind speed east component, m/s

% Computation of wind speed components in the ECEF system, m/s
[VX, VY, VZ]=ned2ecefv(Vw_north, Vw_east, 0, LAT, LONG);

% Save atmospheric ECEF profile in current "atmospheric_data.txt", file that will be read by ODEsystem.m
fid2 = fopen('atmospheric_data.txt','w');
for i=1:length(quota_vento)
fprintf(fid2,'%10.5f \t %10.5f \t %10.6f  %10.5f  %10.5f \t %10.5f\n', quota_vento(i), temp(i), density(i), VX(i), VY(i), VZ(i)); 
end

% Close the file atmospheric_profile.txt
fclose(fid2);

% Copy the current atmospheric data file in the current fireball's folder
copyfile("atmospheric_data.txt", fireball_results)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute meteoroid starting speed in ECEF reference system %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Compute meteoroid starting speed in ECEF reference system')
disp('   ')

V_north=-V0*cosd(I0)*cosd(A0); % Speed component towards north, m/s
V_east=-V0*cosd(I0)*sind(A0);  % Speed component towards east, m/s
V_down=V0*sind(I0);            % Speed component down, m/s

[VX0, VY0, VZ0]=ned2ecefv(V_north, V_east, V_down, LAT, LONG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute meteoroid starting position in ECEF reference system %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Compute meteoroid starting position in ECEF reference system')
disp('   ')

% Compute ECEF coordinates of starting point, m 
[X0, Y0, Z0] = geodetic2ecef(wgs84,LAT,LONG,Q0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vector of starting conditions: speed (m/s) and position (m) in ECEF system, mass (kg) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Define starting conditions vector')
disp('   ')

y0=[VX0 VY0 VZ0 X0 Y0 Z0 M0];

% Integration time vector, s
tspan=0:0.01:T0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical integration of the equations of motion %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Integrate motion equations original meteoroid')
disp('   ')

% Runge-Kutta 4th/5th order ODE solver
opt = odeset('Events',@EventsFunction); % Condition to exit the integration
[t, y]=ode45(@ODEsystem, tspan, y0, opt);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract and save physical solution until main fragmentation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('Extract physical solution until main fragmentation')
disp('   ')
[N, quota_meteoroide, speed, lat, long, height, Dynamic_pressure, residual_mass, diameter]=physical_solution2(fireball_results, Qfin, M0, y, t);
t=t(1:N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical integration for all fragments  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Starting position and speed for fragments
VX1=y(N,1);
VY1=y(N,2);
VZ1=y(N,3);
X1=y(N,4);
Y1=y(N,5);
Z1=y(N,6);

for J=1:length(M1)

disp(strcat('Integrate motion equations of the fragment n.', num2str(J)))
disp('   ')    
    
% Output file with single fragment's path 
fid4 = fopen(strcat(fireball_results, 'Asteroid_fragment_n', num2str(J), '.txt'),'w'); 
fprintf(fid4, '%%    \n');
fprintf(fid4, strcat('%%', " ", fireball_name," ", 'FRAGMENTATION WITH NO PANCAKE PHASE \n'));
fprintf(fid4, '%% Starting strength (MPa)                  %10.2f \n', Strength);
fprintf(fid4, '%% Starting Mass (kg)                       %10.3f \n', M1(J));

% Compute fragment's diameter and strength
D1=2*(3*M1(J)/(4*pi*rho_M))^(1/3);  % Fragment's diameter, m
Strength1=Strength*(M0/M1(J))^0.20; % Fragment's strength with alpha=0.2, MPa (Weibull's law)

% Superior limit to fragment's strength (MPa)
if Strength1 > 100
    Strength1=100;
end

y1=[VX1 VY1 VZ1 X1 Y1 Z1 M1(J)];

fid0 = fopen('meteoroid_data.txt','w');
fprintf(fid0,'%10.5f \t %10.5f \t %10.5f \t %10.5f \t %10.5f\n', F, Ga, rho_M, Strength1, D1);   
fclose(fid0);

% Runge-Kutta 4th/5th order ODE solver
opt = odeset('Events',@EventsFunction);
[t2, y2]=ode45(@ODEsystem, tspan, y1, opt);

[N2, quota_meteoroide2, speed2, lat2, long2, height2, Dynamic_pressure2, residual_mass2, diameter2]=physical_solution(fireball_results, Qfin, M0, y2, t2+t(N), J);

fclose(fid1);

end

% Load model whole data
model_data=strcat(fireball_results, '/Asteroid_model.txt');
s1=model_data;
Dat1(:,:)=load(s1, '-ascii');

Height_tot=Dat1(:,1);           % km     
Speed_tot=Dat1(:,2);            % km/s       
Lat_tot=Dat1(:,3);              % deg 
Long_tot=Dat1(:,4);             % deg     
Dynamic_pressure_tot=Dat1(:,5); % MPa    
Relative_mass_tot=Dat1(:,6);    % Adimensional 
Diameter_tot=Dat1(:,7);         % m
Time_tot=Dat1(:,8);             % s
N_tot=length(Time_tot);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute energy deposition in atmosphere %
% and fireball's absolute magnitude       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Salvataggio dei valori dell'energy deposition e della magnitudine assoluta in un file di output
% fid3 = fopen(strcat(fireball_results, 'Asteroid_energy_deposition_absolute_mag.txt'),'w');
% fprintf(fid3, strcat('%%', " ", fireball_name," ", 'METEOROID FRAGMENTATION WITH NO PANCAKE MODEL \n'));
% fprintf(fid3, '%% Starting Mass (kg)                       %10.2f \n', M0);
% fprintf(fid3, '%% Diameter (m)                             %10.2f \n', 2*R0);
% fprintf(fid3, '%% Mean density (kg/m^3)                    %10.2f \n', rho_M);
% fprintf(fid3, '%% M/A (kg/m^2)                             %10.2f \n', M0/S0);
% fprintf(fid3, '%% Asintotic drag coefficient               %10.2f \n', Ga);
% fprintf(fid3, '%% Starting speed (km/s)                    %10.2f \n', V0/1000);
% fprintf(fid3, '%% Starting inclination (degree)            %10.2f \n', I0);
% fprintf(fid3, '%% Starting height (km)                     %10.2f \n', Q0/1000);
% fprintf(fid3, '%% Starting strength (MPa)                  %10.2f \n', Strength);
% fprintf(fid3, '%% Luminous Efficiency                      %10.2f \n', LE);
% fprintf(fid3, '%%    \n');
% fprintf(fid3, '%% Height (km)     Energy dep (kt/km)       Absolute mag     Time (s)\n');
% fprintf(fid3, '    \n');
% 
% % Original body energy deposition
% [Energy_dep_main, Height_main, Time_main, Power_main]=energy_deposition(Time_tot(1:1:N), Speed_tot(1:1:N), Height_tot(1:1:N), Relative_mass_tot(1:1:N), M0);
% 
% % Efficienza luminosa fissa e pari alla frazione LE della perdita di energia cinetica
% % 1500 W è la costante di normalizzazione per la mag assoluta pari a zero.
% MAG_main=-2.5*log10(LE*Power_main/1500);
% 
% for i=1:length(Time_main)
%      fprintf(fid3,'%10.6f \t  %10.6f\t\t %10.6f  \t\t %10.6f \n', Height_main(i), Energy_dep_main(i), MAG_main(i), Time_main(i)); 
% end
% 
% % Fragments energy deposition
% Energy_dep_frag_tot=0;
% Power_tot_frag=0;
% % Select minimum array length between fragments
% s0=strcat(fireball_results,'Asteroid_fragment_n', '1', '.txt');
% Dat_frag=load(s0, '-ascii');
% Time_frag=Dat_frag(:,8);   
% MM=length(Time_frag); % Minimum array length
% 
% for J=1:length(M1)
%     s0=strcat(fireball_results,'Asteroid_fragment_n', num2str(J), '.txt');
%     Dat_frag=load(s0, '-ascii');
%     Time_frag=Dat_frag(:,8);   % s
%     Speed_frag=Dat_frag(:,2);  % km/s
%     Height_frag=Dat_frag(:,1); % km
%     Relative_mass_frag=Dat_frag(:,6);
%     [Energy_dep_frag, Height_frag, Time_frag1, Power_frag]=energy_deposition(Time_frag(2:1:MM), Speed_frag(2:1:MM), Height_frag(2:1:MM), Relative_mass_frag(2:1:MM), M0);
%     Energy_dep_frag_tot=Energy_dep_frag_tot+Energy_dep_frag;
%     Power_tot_frag=Power_tot_frag+Power_frag;
% end
% 
% % Efficienza luminosa fissa e pari alla frazione LE della perdita di energia cinetica
% % 1500 W è la costante di normalizzazione per la mag assoluta pari a zero.
% MAG_frag=-2.5*log10(LE*Power_frag/1500); 
% 
% % MAG2=MAG(ix:1:end);
% % Time_tot2=Time_tot1(ix:1:end);
%   
% for i=1:length(Time_frag1)
%      fprintf(fid3,'%10.6f \t  %10.6f\t\t %10.6f  \t\t %10.6f \n', Height_frag(i), Energy_dep_frag_tot(i), MAG_frag(i), Time_frag1(i)); 
% end

%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
% Plot model results %
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%

disp(strcat('Plot and save', " ", fireball_name," ", 'model results'))
disp('   ')

% Compute distance range from starting point, km
[arc,az] = distance(LAT, LONG, Lat_tot, Long_tot);
range_tot=111.19*arc;

%========================================================================

figure % Plot generale del modello

% Up figure

subplot(2,3,1)
% Speed vs time
semilogy(Time_tot, Speed_tot, 'k.', 'MarkerSize', 5)
hold on
semilogy(Time_tot(N), Speed_tot(N), 'r.', 'MarkerSize', 30) % Start fragmentation
grid
legend("On air", "Fragmentation",'interpreter','latex')
title('Speed vs time', 'FontSize',14) 
xlabel('Time (s)','FontSize',14)
ylabel('Speed (km/s)','Fontsize',14)
hold off

subplot(2,3,2)
% Speed vs height
loglog(Speed_tot, Height_tot, 'k.', 'MarkerSize', 5)
hold on
loglog(Speed_tot(N), Height_tot (N), 'r.', 'MarkerSize', 30) % Start fragmentation
grid
legend("On air", "Fragmentation",'interpreter','latex')
title('Speed vs height', 'FontSize',14) 
xlabel('Speed (km/s)','FontSize',14)
ylabel('Height (km)','Fontsize',14)
hold off


subplot(2,3,3)
% Height vs time
semilogx(Time_tot, Height_tot, 'k.', 'MarkerSize', 5)
hold on
semilogx(Time_tot(N), Height_tot(N), 'r.', 'MarkerSize', 30) % Start fragmentation
grid
legend("On air", "Fragmentation",'interpreter','latex')
title('Height vs time', 'FontSize',14)
xlabel('Time (s)','FontSize',14)
ylabel('Height (km)','Fontsize',14)
hold off


% Down figure

subplot(2,3,4)
% Dynamic pressure vs height
plot(Height_tot, Dynamic_pressure_tot, 'k.', 'MarkerSize', 5)
hold on
plot(Height_tot(N), Dynamic_pressure_tot(N), 'r.', 'MarkerSize', 30) % Start fragmentation
hold on
grid
legend("On air", "Fragmentation",'interpreter','latex')
title('Dynamic pressure vs height', 'FontSize',14)
xlabel('Height (km)','FontSize',14)
ylabel('Dynamic pressure (MPa)','Fontsize',14)
hold off

subplot(2,3,5)
% Height vs range
axis equal
plot(range_tot, Height_tot, 'k.', 'MarkerSize', 5)
hold on
plot(range_tot(N), Height_tot(N), 'r.', 'MarkerSize', 30) % Start fragmentation
hold on
grid
legend("On air", "Fragmentation",'interpreter','latex')
title('Height vs range', 'FontSize',14)
xlabel('Range (km)','FontSize',14)
ylabel('Height (km)','Fontsize',14)
hold off

subplot(2,3,6)
% Dynamic pressure vs speed
plot(Speed_tot, Dynamic_pressure_tot, 'k.', 'MarkerSize', 5)
hold on
plot(Speed_tot(N), Dynamic_pressure_tot(N), 'r.', 'MarkerSize', 30) % Start fragmentation
hold on
grid
legend("On air", "Fragmentation",'interpreter','latex')
title('Dynamic pressure vs speed', 'FontSize',14)
xlabel('Speed (km/s)','FontSize',14)
ylabel('Dynamic pressure (MPa)','Fontsize',14)
hold off

suptitle('Asteroid fall model with meteoroid disintegration')

hold off

file_output=strcat(fireball_results,'Asteroid_model_',fireball_name, '_', num2str(Strength), 'MPa','.fig');
saveas(gcf, file_output, 'fig')   % Save the figure in .fig format

%========================================================================

figure % Figure of the paper

% Compute local minimum of the total height.
% Each minimum correspond to a meteorite on the ground.
TF1 = islocalmin(Height_tot,'FlatSelection','first');
IN = find(TF1==1);
IN(end+1) = length(Height_tot); % Add last minimum

% Restore final mass value to plot in the graph
M1=f*M1;
nn=1;

subplot(2,2,1)
% Height vs range
axis equal
plot(range_tot(1:1:N), Height_tot(1:1:N), 'k-', 'LineWidth', 2, 'DisplayName', sprintf('Original mass'))
hold on
plot(range_tot(N), Height_tot(N), 'r.', 'MarkerSize', 30, 'DisplayName', sprintf('Fragmentation')) % Start fragmentation
plot(range_tot(N:IN(1)), Height_tot(N:IN(1)), '-', 'LineWidth', 2, 'DisplayName', sprintf(strcat('Final mass'," ", num2str(M1(1)), 'kg'))) % Plot first meteorite
for J=1:length(M1)-1 % plot of subsequent meteorites
    plot(range_tot(IN(J)+nn:IN(J+1)-nn), Height_tot(IN(J)+nn:IN(J+1)-nn), '-', 'LineWidth', 2, 'DisplayName', sprintf(strcat('Final mass'," ", num2str(M1(J+1)), 'kg')))
end
plot(range_tot(N), Height_tot(N), 'r.', 'MarkerSize', 30, 'HandleVisibility','off') % Start fragmentation
grid
axis equal
legend('-DynamicLegend');
legend('show');
ax = gca;
ax.FontSize = 10;  % Font Size of 10
ylim([0 Height_tot(N)+5])
%title('Height vs range', 'FontSize',14)
xlabel('Range (km)','FontSize',20)
ylabel('Height (km)','Fontsize',20)
hold off


subplot(2,2,2)
% Height vs speed
semilogx(Speed_tot(1:1:N), Height_tot(1:1:N),  'k-', 'LineWidth', 2, 'DisplayName', sprintf('Original mass'))
hold on
semilogx(Speed_tot(N), Height_tot(N), 'r.', 'MarkerSize', 30, 'DisplayName', sprintf('Fragmentation')) % Start fragmentation
semilogx(Speed_tot(N:IN(1)), Height_tot(N:IN(1)), '-', 'LineWidth', 2, 'DisplayName', sprintf(strcat('Final mass'," ", num2str(M1(1)), 'kg'))) % Plot first meteorite
for J=1:length(M1)-1 % plot of subsequent meteorites
    semilogx(Speed_tot(IN(J)+nn:IN(J+1)-nn), Height_tot(IN(J)+nn:IN(J+1)-nn), '-', 'LineWidth', 2, 'DisplayName', sprintf(strcat('Final mass'," ", num2str(M1(J+1)), 'kg')))
end
semilogx(Speed_tot(N), Height_tot(N), 'r.', 'MarkerSize', 30, 'HandleVisibility','off') % Start fragmentation
grid
%legend('-DynamicLegend');
%legend('show');
ax = gca;
ax.FontSize = 10;  % Font Size of 10
%title('Height vs speed', 'FontSize',14)
ylim([0 Height_tot(N)+5])
xlabel('Speed (km/s)','FontSize',20)
ylabel('Height (km)','Fontsize',20)
hold off


subplot(2,2,3)
% Dynamic pressure vs height
semilogx(Dynamic_pressure_tot(1:1:N), Height_tot(1:1:N), 'k-', 'LineWidth', 2, 'DisplayName', sprintf('Original mass'))
hold on
semilogx(Dynamic_pressure_tot(N), Height_tot(N), 'r.', 'MarkerSize', 30, 'DisplayName', sprintf('Fragmentation')) % Start fragmentation
semilogx(Dynamic_pressure_tot(N:IN(1)), Height_tot(N:IN(1)), '-', 'LineWidth', 2, 'DisplayName', sprintf(strcat('Final mass'," ", num2str(M1(1)), 'kg'))) % Plot first meteorite
for J=1:length(M1)-1 % plot of subsequent meteorites
    semilogx(Dynamic_pressure_tot(IN(J)+nn:IN(J+1)-nn), Height_tot(IN(J)+nn:IN(J+1)-nn), '-', 'LineWidth', 2, 'DisplayName', sprintf(strcat('Final mass'," ", num2str(M1(J+1)), 'kg')))
end
semilogx(Dynamic_pressure_tot(N), Height_tot(N), 'r.', 'MarkerSize', 30, 'HandleVisibility','off') % Start fragmentation
grid
%legend('-DynamicLegend');
%legend('show');
ax = gca;
ax.FontSize = 10;  % Font Size of 10
ylim([0 Height_tot(N)+5])
%title('Height vs Dynamic pressure', 'FontSize',14)
xlabel('Dynamic pressure (MPa)','FontSize',20)
ylabel('Height (km)','Fontsize',20)
hold off


subplot(2,2,4)
% Height vs mass
semilogx(M0*Relative_mass_tot(1:1:N), Height_tot(1:1:N), 'k-', 'LineWidth', 2, 'DisplayName', sprintf('Original mass'))
hold on
semilogx(M0*Relative_mass_tot(N), Height_tot(N), 'r.', 'MarkerSize', 30, 'DisplayName', sprintf('Fragmentation')) % Start fragmentation
semilogx(M0*Relative_mass_tot(N:IN(1)), Height_tot(N:IN(1)), '-', 'LineWidth', 2, 'DisplayName', sprintf(strcat('Final mass'," ", num2str(M1(1)), 'kg'))) % Plot first meteorite
for J=1:length(M1)-1 % plot of subsequent meteorites
    semilogx(M0*Relative_mass_tot(IN(J)+nn:IN(J+1)-nn), Height_tot(IN(J)+nn:IN(J+1)-nn), '-', 'LineWidth', 2, 'DisplayName', sprintf(strcat('Final mass'," ", num2str(M1(J+1)), 'kg')))
end
semilogx(M0*Relative_mass_tot(N), Height_tot(N), 'r.', 'MarkerSize', 30, 'HandleVisibility','off') % Start fragmentation
grid
%legend('-DynamicLegend');
%legend('show');
ax = gca;
ax.FontSize = 10;  % Font Size of 10
ylim([0 Height_tot(N)+5])
%title('Height vs mass', 'FontSize',14)
xlabel('Mass (kg)','FontSize',20)
ylabel('Height (km)','Fontsize',20)
hold off

file_output1=strcat(fireball_results,'Asteroid_model_for_paper_',fireball_name, '_', num2str(Strength), 'MPa','.fig');
saveas(gcf, file_output1, 'fig')   % Save the figure in .fig format

%========================================================================

% Figure for the paper
% 
% % Load energy data from file
% model_data1=strcat(fireball_results, '/Asteroid_energy_deposition_absolute_mag.txt');
% s2=model_data1;
% Dat2(:,:)=load(s2, '-ascii');
% 
% Height_energy=Dat2(:,1);     % km     
% Energy_dep_tot=Dat2(:,2);    % kt/km       
% Absolute_mag_tot=Dat2(:,3);  % mag 
% Time_energy=Dat2(:,4);       % s    
% 
% subplot(1,2,1)
% 
% % Absolute mag vs time
% semilogx(Time_energy, Absolute_mag_tot, 'k.', 'MarkerSize', 5)
% hold on
% set(gca,'Ydir','reverse')
% hold on
% grid
% title('Absolute mag vs time', 'FontSize',14)
% xlabel('Time (s)','FontSize',14)
% ylabel('Absolute mag','Fontsize',14)
% hold off
% 
% subplot(1,2,2)
% 
% % Height vs energy deposition
% semilogx(Energy_dep_tot, Height_energy, 'k.', 'MarkerSize', 5)
% hold on
% grid
% axis([Energy_dep_tot(1) max(Energy_dep_tot) 10 75])
% title('Height vs energy deposition', 'FontSize',14)
% xlabel('Energy deposition (kt/km)','FontSize',14)
% ylabel('Height (km)','Fontsize',14)
% hold off

%========================================================================

disp(strcat('End'," ", fireball_name," ", 'model computation'))
disp('   ')
fclose('all');