% Meteorite fragmentation for 2023 CX1 (SPECIAL EDITION WITH STRENGTH)
% Matlab script for the computation of an atmospheric meteoroid path in a 3D ECEF system
% with main explosion, pancake phase end fragment fall. 
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
% Version Feb 23, 2024

clear all

disp('                                                                       ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%                      METEORITE MODEL FOR 2023 CX1                   %')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%                     PATH AND FRAGMENT COMPUTATION                   %')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('                                                                       ')
disp('                    by Albino Carbognani (INAF-OAS)                    ')
disp('                              Feb 2024                                 ')
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

wholefile_set = fileread('.\Settings.txt');
% Split of the Settings and initialization of the string variables
set = regexp(wholefile_set,'\$+','split');

fireball_name=strtrim(set{4});         % Fireball's folder name
WP=strtrim(set{7});                    % Working path of the script
ATMO=strtrim(set{10});                 % Atmospheric profile name
E0=str2double(strtrim(set{13}));       % Meteoroid kinetic energy, Mt
rho_M=str2double(strtrim(set{16}));    % Mean density meteoroid, kg/m^3
V0=str2double(strtrim(set{19}));       % Starting speed, m/s
I0=str2double(strtrim(set{22}));       % Trajectory inclination at starting speed, degrees
A0=str2double(strtrim(set{25}));       % Trajectory azimut at starting speed, degrees
Q0=str2double(strtrim(set{28}));       % Starting height above Earth surface, m
LAT=str2double(strtrim(set{31}));      % Latitude starting point, degrees
LONG=str2double(strtrim(set{34}));     % Longitude starting point, degrees
Qfin=str2double(strtrim(set{37}));     % Mean elevation of the ground (estimated) on which the meteoroid falls, m
T0=str2double(strtrim(set{40}));       % Maximum integration time, s
Ga=str2double(strtrim(set{43}));       % Asintotic drag coefficient
F=str2double(strtrim(set{46}));        % Dimensionless spherical form factor
Strength=str2double(strtrim(set{49})); % Mean body strength (MPa)
M1=str2double(strtrim(set{52}));       % Fragment's mass, kg

% Copy the current settings file in the current fireball's folder
copyfile("Settings.txt", fireball_name)

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
% EventsFunction.m and EventsFunction_L.m
fid0 = fopen('meteoroid_data.txt','w');
fprintf(fid0,'%10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f\n', F, Ga, rho_M, Strength, R0);   
fclose(fid0);

% Copy the current meteoroid file in the current fireball's folder
copyfile("meteoroid_data.txt", fireball_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define path and atmospheric data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Define working path and atmospheric data location')
disp('   ')

working_path=strcat(WP, fireball_name);  
atmospheric_profile=strcat(working_path, '/', ATMO);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open the output model file %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid1 = fopen(strcat(working_path, '/Asteroid_model.txt'),'w'); 
fprintf(fid1, '%%    \n');
fprintf(fid1, '%% 2023 CX1 SINGLE BODY MODEL \n');
fprintf(fid1, '%%    \n');
fprintf(fid1, '%% Fireball name                            %s \n', '2023 CX1');
fprintf(fid1, '%% Starting Mass (kg)                       %10.2f \n', M0);
fprintf(fid1, '%% Diameter (m)                             %10.2f \n', 2*R0);
fprintf(fid1, '%% Mean density (kg/m^3)                    %10.0f \n', rho_M);
fprintf(fid1, '%% M/A (kg/m^2)                             %10.2f \n', M0/S0);
fprintf(fid1, '%% Asintotic drag coefficient               %10.2f \n', Ga);
fprintf(fid1, '%% Starting speed (km/s)                    %10.2f \n', V0/1000);
fprintf(fid1, '%% Starting inclination (degree)            %10.2f \n', I0);
fprintf(fid1, '%% Starting height (km)                     %10.2f \n', Q0/1000);
fprintf(fid1, '%% Starting strength (MPa)                  %10.2f \n', Strength);
fprintf(fid1, '    \n');
fprintf(fid1, '%% Height (km)     Speed (km/s)       Lat (degrees)     Long (degrees)    Dynamic pressure (MPa)    Relative mass    Diameter (m)   Time (s)\n');
fprintf(fid1, '    \n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute and save "atmospheric_data.txt" from ATMO %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
string1=strcat('Compute and save ECEF atmospheric data from', " ", ATMO);
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
fprintf(fid2,'%10.3f \t %10.3f \t %10.6f  %10.3f  %10.3f \t %10.3f\n', quota_vento(i), temp(i), density(i), VX(i), VY(i), VZ(i)); 
end

% Close the file atmospheric_profile.txt
fclose(fid2);

% Copy the current atmospheric data file in the current fireball's folder
copyfile("atmospheric_data.txt", fireball_name)

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
[N, quota_meteoroide, speed, lat, long, height, Dynamic_pressure, residual_mass, diameter]=physical_solution2(working_path, Qfin, M0, y, t);
t=t(1:N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical integration for lateral expansion (pancake phase) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Integrate motion equations lateral expansion (pancake phase)')
disp('   ')

VXL=y(N,1);
VYL=y(N,2);
VZL=y(N,3);
XL=y(N,4);
YL=y(N,5);
ZL=y(N,6);

yL=[VXL VYL VZL XL YL ZL M0*residual_mass(N) 0];

% Runge-Kutta 4th/5th order ODE solver
optL = odeset('Events',@EventsFunction_L); % Condition to exit the integration
[t_lat, y_lat]=ode45(@ODEsystem_L, tspan, yL, optL);

[NL, quota_meteoroideL, speedL, latL, longL, heightL, Dynamic_pressureL, residual_massL, diameterL]=physical_solution_L(working_path, Qfin, M0, y_lat, t_lat+t(N));
t_lat=t_lat(1:NL);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical integration for fragment path %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Integrate fragment''s motion equations ')
disp('   ')

D1=2*(3*M1/(4*pi*rho_M))^(1/3);  % Fragment's diameter, m
Strength1=Strength*(M0/M1)^0.20; % Fragment's strength with alpha=0.2, MPa (Weibull's law)

% Superior limit to fragment's strength (MPa)
if Strength1 > 100
    Strength1=100;
end

VX1=y_lat(NL,1);
VY1=y_lat(NL,2);
VZ1=y_lat(NL,3);
X1=y_lat(NL,4);
Y1=y_lat(NL,5);
Z1=y_lat(NL,6);

y1=[VX1 VY1 VZ1 X1 Y1 Z1 M1];

fid0 = fopen('meteoroid_data.txt','w');
fprintf(fid0,'%10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f\n', F, Ga, rho_M, Strength1, D1);   
fclose(fid0);

% Runge-Kutta 4th/5th order ODE solver
opt = odeset('Events',@EventsFunction);
[t2, y2]=ode45(@ODEsystem, tspan, y1, opt);

[N2, quota_meteoroide2, speed2, lat2, long2, height2, Dynamic_pressure2, residual_mass2, diameter2]=physical_solution(working_path, Qfin, M0, y2, t2+t_lat(NL)+t(N));

fclose(fid1);

%%%%%%%%%%%%%%%%
% Plot results %
%%%%%%%%%%%%%%%%

% Load model whole data
model_data=strcat(working_path, '/Asteroid_model.txt');
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

% Compute energy deposition, Mt/km
% Kinetic energy, Mt
K=(0.5*M0*Relative_mass_tot.*(1000*Speed_tot).^2)/(4.184*10^15); 
dK=gradient(K); % Mt
dHeight=gradient(Height_tot); % km

Energy_dep=abs(dK./dHeight); % Mt/km

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save Google Earth files, atmospheric path, dark flight and strewn field %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Save atmospheric path in .kml files for Google Earth')
disp('   ')

% Syntax: kmlwriteline(filename,latitude,longitude)
% kmlwriteline(filename,latitude,longitude,altitude)
iconFilename1 = 'http://maps.google.com/mapfiles/kml/shapes/placemark_circle_highlight.png';
iconFilename2 = 'http://maps.google.com/mapfiles/kml/shapes/star.png';
iconFilename3 = 'http://maps.google.com/mapfiles/kml/paddle/ylw-blank-lv.png';
iconFilename4 = 'https://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png';
% Array of white characters to not show the name near the icons in the kml file
Name_blank = blanks(length(Lat_tot));
filename1 = strcat(working_path, '/', '2023CX1_atmospheric_path_', num2str(M1), 'kg_', num2str(Strength), 'MPa', '.kml');
kmlwritepoint(filename1, Lat_tot(1:10:length(Lat_tot))', Long_tot(1:10:length(Long_tot))', 1000*Height_tot(1:10:length(Height_tot))', 'Icon', iconFilename1, 'IconScale', 1, 'Name', Name_blank);

% Reads filename1 and adds the line that reaches the ground into the trajectory kml file
filename2 = [working_path '/','2023CX1_atmospheric_path_and_height_', num2str(M1), 'kg_', num2str(Strength), 'MPa','.kml'];
fidi=fopen(filename1,'r');
fido=fopen(filename2,'w');
while ~feof(fidi)
  l=fgetl(fidi);           % read line
  if strfind(l,'</Point>') % Stringa da sostituire
    % modify line here
    l='<extrude>1</extrude></Point>';
  end
  fprintf(fido,'%s',l);  % 'fgetl returns \n so it's embedded
end
fidi=fclose(fidi);
fido=fclose(fido);

% Save strewn field data for Google Earth
filename3 = strcat(working_path, '/', '2023CX1_strewn_field_', num2str(M1), 'kg_', num2str(Strength), 'MPa.kml');
Name_strewn=strcat('2023CX1,', " ",'Lat'," ", num2str(Lat_tot(length(Lat_tot))),", ", 'Long'," ", num2str(Long_tot(length(Long_tot))), " ", num2str(M1), 'kg'," ", '(', num2str(Strength), 'MPa)'); 
kmlwritepoint(filename3, Lat_tot(length(Lat_tot)), Long_tot(length(Long_tot)), 'Icon', iconFilename4, 'IconScale', 1.5, 'Name', Name_strewn);


%%%%%%%%%%%%%%%%%%%%%%%%%%
% End Google Earth files %
%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Plot and save 2023 CX1 model results')
disp('   ')

% Compute distance range from starting point, km
[arc,az] = distance(LAT, LONG, Lat_tot, Long_tot);
range_tot=111.19*arc;

figure

% Up figure

subplot(2,4,1)
% Speed vs time
plot(Time_tot, Speed_tot, 'k.', 'MarkerSize', 10)
hold on
plot(Time_tot(N), Speed_tot(N), 'r.', 'MarkerSize', 30) % Start fragmentation
hold on
plot(Time_tot(N+NL), Speed_tot(N+NL), 'g*', 'MarkerSize', 20) % End pancake phase
grid
legend("On air", "Start fragmentation","Airburst",'interpreter','latex')
title('Speed vs time', 'FontSize',14) 
xlabel('Time (s)','FontSize',14)
ylabel('Velocity (km/s)','Fontsize',14)
hold off

subplot(2,4,2)
% Residual mass vs time
plot(Time_tot, Relative_mass_tot, 'k.', 'MarkerSize', 10)
hold on
plot(Time_tot(N), Relative_mass_tot(N), 'r.', 'MarkerSize', 30) % Start fragmentation
hold on
plot(Time_tot(N+NL), Relative_mass_tot(N+NL), 'g*', 'MarkerSize', 20) % End pancake phase
hold on
grid
legend("On air", "Start fragmentation","Airburst",'interpreter','latex')
title('Mass vs time', 'FontSize',14)
xlabel('Time (s)','FontSize',14)
ylabel('Relative mass','Fontsize',14)
hold off

subplot(2,4,3)
% Height vs time
plot(Time_tot, Height_tot, 'k.', 'MarkerSize', 10)
hold on
plot(Time_tot(N), Height_tot(N), 'r.', 'MarkerSize', 30) % Start fragmentation
hold on
plot(Time_tot(N+NL), Height_tot(N+NL), 'g*', 'MarkerSize', 20) % End pancake phase
grid
legend("On air", "Start fragmentation","Airburst",'interpreter','latex')
title('Height vs time', 'FontSize',14)
xlabel('Time (s)','FontSize',14)
ylabel('Height (km)','Fontsize',14)
hold off

subplot(2,4,4)
% Diameter vs time
plot(Time_tot, Diameter_tot, 'k.', 'MarkerSize', 10)
hold on
plot(Time_tot(N), Diameter_tot(N), 'r.', 'MarkerSize', 30) % Start fragmentation
hold on
plot(Time_tot(N+NL), Diameter_tot(N+NL), 'g*', 'MarkerSize', 20) % End pancake phase
grid
legend("On air", "Start fragmentation","Airburst",'interpreter','latex')
title('Diameter vs time', 'FontSize',14)
xlabel('Time (s)','FontSize',14)
ylabel('Diameter (m)','Fontsize',14)
hold off

% Down figure

subplot(2,4,5)
% Dynamic pressure vs height
plot(Height_tot, Dynamic_pressure_tot, 'k.', 'MarkerSize', 10)
hold on
plot(Height_tot(N), Dynamic_pressure_tot(N), 'r.', 'MarkerSize', 30) % Start fragmentation
hold on
plot(Height_tot(N+NL), Dynamic_pressure_tot(N+NL), 'g*', 'MarkerSize', 20) % End pancake phase
hold on
grid
legend("On air", "Start fragmentation","Airburst",'interpreter','latex')
title('Dynamic pressure vs height', 'FontSize',14)
xlabel('Height (km)','FontSize',14)
ylabel('Dynamic pressure (MPa)','Fontsize',14)
hold off

subplot(2,4,6)
% Height vs range
plot(range_tot, Height_tot, 'k.', 'MarkerSize', 10)
hold on
plot(range_tot(N), Height_tot(N), 'r.', 'MarkerSize', 30) % Start fragmentation
hold on
plot(range_tot(N+NL), Height_tot(N+NL), 'g*', 'MarkerSize', 20) % End pancake phase
hold on
grid
legend("On air", "Start fragmentation","Airburst",'interpreter','latex')
title('Height vs range', 'FontSize',14)
xlabel('Range (km)','FontSize',14)
ylabel('Height (km)','Fontsize',14)
hold off

subplot(2,4,7)
% % Path latitude vs longitude
% plot(Long_tot, Lat_tot, 'k.', 'MarkerSize', 10)
% hold on
% plot(Long_tot(N), Lat_tot(N), 'r.', 'MarkerSize', 30) % Start fragmentation
% hold on
% plot(Long_tot(N+NL), Lat_tot(N+NL), 'g*', 'MarkerSize', 20) % End pancake phase
% grid
% legend("On air", "Start fragmentation","Airburst",'interpreter','latex')
% title('Path on the ground', 'FontSize',14)
% xlabel('Long (deg)','FontSize',14)
% ylabel('Lat (deg)','Fontsize',14)
% hold off

% Height vs energy deposition
%semilogx(Energy_dep, Height_tot, '--', 'MarkerSize', 10)
plot(Energy_dep, Height_tot, '--', 'MarkerSize', 10)
hold on
grid
title('Height vs energy deposition', 'FontSize',14)
xlabel('Energy deposition (Mt/km)','FontSize',14)
ylabel('Height (km)','Fontsize',14)
hold off


subplot(2,4,8)
% Dynamic pressure vs speed
plot(Speed_tot, Dynamic_pressure_tot, 'k.', 'MarkerSize', 10)
hold on
plot(Speed_tot(N), Dynamic_pressure_tot(N), 'r.', 'MarkerSize', 30) % Start fragmentation
hold on
plot(Speed_tot(N+NL), Dynamic_pressure_tot(N+NL), 'g*', 'MarkerSize', 20) % End pancake phase
hold on
grid
legend("On air", "Start fragmentation","Airburst",'interpreter','latex')
title('Dynamic pressure vs speed', 'FontSize',14)
xlabel('Speed (km/s)','FontSize',14)
ylabel('Dynamic pressure (MPa)','Fontsize',14)
hold off

suptitle('Asteroid fall model')

hold off

cd(working_path)
file_output=strcat('Asteroid_model_', num2str(M1), 'kg_', num2str(Strength), 'MPa','.fig');
saveas(gcf, file_output, 'fig')   % Save the figure in .fig format
cd(WP)

figure % Per il paper

subplot(2,2,1)
% Speed vs time
plot(Time_tot, Speed_tot, 'k.', 'MarkerSize', 10)
hold on
plot(Time_tot(N), Speed_tot(N), 'r.', 'MarkerSize', 30) % Start fragmentation
hold on
plot(Time_tot(N+NL), Speed_tot(N+NL), 'g*', 'MarkerSize', 20) % End pancake phase
grid
legend("On air", "Start fragmentation","Airburst",'interpreter','latex')
title('Speed vs time', 'FontSize',14) 
xlabel('Time (s)','FontSize',14)
ylabel('Velocity (km/s)','Fontsize',14)
hold off

subplot(2,2,2)
% Height vs range
plot(range_tot, Height_tot, 'k.', 'MarkerSize', 10)
hold on
plot(range_tot(N), Height_tot(N), 'r.', 'MarkerSize', 30) % Start fragmentation
hold on
plot(range_tot(N+NL), Height_tot(N+NL), 'g*', 'MarkerSize', 20) % End pancake phase
hold on
grid
legend("On air", "Start fragmentation","Airburst",'interpreter','latex')
title('Height vs range', 'FontSize',14)
xlabel('Range (km)','FontSize',14)
ylabel('Height (km)','Fontsize',14)
hold off

subplot(2,2,3)
% Dynamic pressure vs height
plot(Height_tot, Dynamic_pressure_tot, 'k.', 'MarkerSize', 10)
hold on
plot(Height_tot(N), Dynamic_pressure_tot(N), 'r.', 'MarkerSize', 30) % Start fragmentation
hold on
plot(Height_tot(N+NL), Dynamic_pressure_tot(N+NL), 'g*', 'MarkerSize', 20) % End pancake phase
hold on
grid
legend("On air", "Start fragmentation","Airburst",'interpreter','latex')
title('Dynamic pressure vs height', 'FontSize',14)
xlabel('Height (km)','FontSize',14)
ylabel('Dynamic pressure (MPa)','Fontsize',14)
hold off

subplot(2,2,4)
% Dynamic pressure vs speed
plot(Speed_tot, Dynamic_pressure_tot, 'k.', 'MarkerSize', 10)
hold on
plot(Speed_tot(N), Dynamic_pressure_tot(N), 'r.', 'MarkerSize', 30) % Start fragmentation
hold on
plot(Speed_tot(N+NL), Dynamic_pressure_tot(N+NL), 'g*', 'MarkerSize', 20) % End pancake phase
hold on
grid
legend("On air", "Start fragmentation","Airburst",'interpreter','latex')
title('Dynamic pressure vs speed', 'FontSize',14)
xlabel('Speed (km/s)','FontSize',14)
ylabel('Dynamic pressure (MPa)','Fontsize',14)
hold off

disp('End 2023 CX1 model computation')
disp('   ')