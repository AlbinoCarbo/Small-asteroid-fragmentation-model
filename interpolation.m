% Function for the interpolation of wind speed, air density and temperature 
% between two different altitudes of the atmospheric profile.
%
% FUNCTION USED:
% vq = interp1(x,v,xq,method,extrapolation) specifies a strategy for evaluating points that lie outside the domain of x. 
% Set extrapolation to 'extrap' when you want to use the method algorithm for extrapolation. 
% Alternatively, you can specify a scalar value, in which case, interp1 returns that value for all points 
% outside the domain of x.
%
% 'pchip' = Shape-preserving piecewise cubic interpolation. 
% The interpolated value at a query point is based on a shape-preserving piecewise cubic interpolation 
% of the values at neighboring grid points.
%
% Input
% I = index of Q nearest to Qmet
% Qmet = meteoroid height (scalar), m
% Q = array winds speed, m
% T = array temperature, K
% rho = array air density, kg/m^3
% vx = Array x component wind speed in ECEF, m/s
% vy = Array y component wind speed in ECEF, m/s
% vz = Array z component wind speed in ECEF, m/s
%
% Output:
% Interpolated values of T (K), rho (kg/m^3), vx (m/s), vy (m/s), vz (m/s)
% at meteoroid height.
%
% Albino Carbognani, INAF-OAS
% Version Feb 7, 2024

function [Ti, rho_i, vxi, vyi, vzi]=interpolation(I, Qmet, Qw, T, rho, vx, vy, vz)

Ti=interp1(Qw, T, Qmet, 'pchip', T(I));         % Interpolate temperature, K
rho_i=interp1(Qw, rho, Qmet, 'pchip', rho(I));  % Interpolate density, kg/m^3
vxi=interp1(Qw, vx, Qmet, 'pchip', vx(I));   % Interpolate Vx, m/s
vyi=interp1(Qw, vy, Qmet, 'pchip', vy(I));   % Interpolate Vy, m/s
vzi=interp1(Qw, vz, Qmet, 'pchip', vz(I));   % Interpolate Vz, m/s

end