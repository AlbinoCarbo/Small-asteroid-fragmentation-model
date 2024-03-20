% System of differential equations for the description of the path of a
% single-body meteoroid in atmosphere with WGS84 Earth model in a
% geocentric Cartesian reference system Earth-centered Earth-fixed (ECECF).
% It needs the atmospheric profile, from the starting height to the ground.
% 
% y(1) = vector speed x, m/s (ECECF)
% y(2) = vector speed y, m/s (ECECF)
% y(3) = vector speed x, m/s (ECECF)
% y(4) = vector position x, m (ECECF)
% y(5) = vector position y, m (ECECF)
% y(6) = vector position z, m (ECECF)
% y(7) = vector mass, kg
%
% The initialization of the constants describing the system takes place
% by reading the "meteoroid_data.txt" and "atmospheric_data.txt" files, 
% created by the main script which must be in the "ODEsystema.m" folder.
%
% Albino Carbognani, INAF-OAS
% Version Feb 7, 2024

function dydt=ODEsystem(t, y)

dydt=zeros(7, 1);

M=5.972*10^(24);  % Earth mass, kg
G=6.672*10^(-11); % Gravitational constant, MKS
Q=8*10^6;         % Heat of ablation for stony asteroid, J/kg
CH=0.1;           % Heat transfer coefficient (adimensional)

% Reference ellipsoid for World Geodetic System of 1984
wgs84 = wgs84Ellipsoid;

    %%%%%%%%%%%%%%%%%%%%%%%
    % Read meteoroid data %
    %%%%%%%%%%%%%%%%%%%%%%%

    met_dat='meteoroid_data.txt';
    s1=met_dat;
    Dat1(:,:)=load(s1, '-ascii');
    F=Dat1(:, 1);         % Dimensionless spherical form factor
    Ga0=Dat1(:, 2);       % Meteoroid's asintotic drag coefficient (adimensional)
    rho_M=Dat1(:, 3);     % Mean density, kg/m^3
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read the atmospheric profile %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    atmospheric_data='atmospheric_data.txt';
    s0=atmospheric_data;
    Dat(:,:)=load(s0, '-ascii');

    quota_vento=Dat(:, 1);    % Height wind, m 
    temperatura=Dat(:, 2);    % Temperature, K 
    density=Dat(:, 3);        % Air density, kg/m^3
    Vx=Dat(:, 4);             % x component wind speed in ECEF, m/s
    Vy=Dat(:, 5)';            % y component wind speed in ECEF, m/s
    Vz=Dat(:, 6)';            % z component wind speed in ECEF, m/s 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute values of temperature, air density, speed components wind %
    % and sound speed for the integration height                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Compute scalar value of meteoroid height, m
    [~,~,quota_meteoroide]=ecef2geodetic(wgs84, y(4), y(5), y(6));
    [~, I]=min(abs(quota_vento-quota_meteoroide));
    
    % Shape-preserving piecewise cubic interpolation of temperature, air density and speed winds
    % component vs "quota_meteoroide".
    [Ti, rho_a, VX, VY, VZ]=interpolation(I, quota_meteoroide, quota_vento, temperatura, density, Vx, Vy, Vz);
    vs=20.0468*sqrt(Ti);   % Speed sound, m/s
    
    % Warning: according to the modern theory of sound, its speed in the Earth’s 
    % atmosphere depends only on temperature and does not depend on its density (height).
    % V. Kirtskhalia 2021 IOP Conf. Ser.: Mater. Sci. Eng. 1024 012037
    
    % OLD ALGORITHM
    % rho_a=density(I);
    % VX=Vx(I);
    % VY=Vy(I);
    % VZ=Vz(I); 
    % vs=20.0468*sqrt(temperatura(I));
             
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COEFFICIENT FOR ATMOSPHERIC DRAG MODEL %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Table Ceplecha, 1987
    % M=[4 3 2 1.5 1.2 1 0.8 0.6 0.4 0.2];                             % Mach number
    % G=[0.580 0.618 0.632 0.596 0.552 0.504 0.441 0.389 0.351 0.328]; % Drag
    %
    % Polynomial function for the drag coefficient in 
    % m^4/kg^2 vs Mach number x (for Mach <= 4).

    G_S=@(x)((13.2566820069031e-003)*x^4-(101.740672494020e-003)*x^3+(186.705667397016e-003)*x^2+(104.950251795627e-003)*x+291.373797865474e-003);       

    v=sqrt(y(1)^2+y(2)^2+y(3)^2); % Meteoroid speed respect to soil, m/s
    
    if v/vs <=4
    Ga=G_S(v/vs);     % Drag for low Mach numbers
    end
    
    if v/vs > 4
    Ga=Ga0;           % Drag for high Mach numbers
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differential equations of motion %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OM=2*pi/86400;                                 % Earth's angular speed, s^-1
r=sqrt(y(4)^2+y(5)^2+y(6)^2);                  % Meteoroid radial distance, m
Vr=sqrt((y(1)-VX)^2+(y(2)-VY)^2+(y(3)-VZ)^2);  % Modulus of relative velocity with air, m/s

% Compute effective cross section meteoroid, m^2
S0=F*(y(7)/rho_M)^(2/3);

% Acceleration along x, m/s^2 (include Coriolis acceleration)
dydt(1)=-G*M*y(4)/(r^3)-(Ga*rho_a*Vr*S0*(y(1)-VX))/y(7)+(OM^2)*y(4)+2*OM*y(2);

% Acceleration along y, m/s^2 (include Coriolis acceleration)
dydt(2)=-G*M*y(5)/(r^3)-(Ga*rho_a*Vr*S0*(y(2)-VY))/y(7)+(OM^2)*y(5)-2*OM*y(1);

% Acceleration along z, m/s^2
dydt(3)=-G*M*y(6)/(r^3)-(Ga*rho_a*Vr*S0*(y(3)-VZ))/y(7);

% Speed x, m/s
dydt(4)=y(1);

% Speed y, m/s
dydt(5)=y(2);

% Speed z, m/s
dydt(6)=y(3);

% Mass ablation, kg/s
dydt(7)=-Ga*(CH/Q)*rho_a*S0*Vr*Vr*Vr;

% Condition for ablation end
if Vr <= 3000
    dydt(7)=0;
end

end