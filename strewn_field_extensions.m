% Matlab function to compute strewn field extension
%
% Input
% LAT, LONG = coordinates strewn field nominal solution (degrees)
% a_R = fireball's azimut (degrees)
% Lat_impact, Long_impact = coordinates of possible impact points (degrees)
%
% Output
% strewn_field_X = strewn field in orthogonal direction coordinates, relative to LAT, LONG point, km
% sigma_X = standard deviation in orthogonal direction, km 
% strewn_field_Y = strewn field in parallel direction coordinates, relative to LAT, LONG point, km 
% sigma_Y = standard deviation in parallel direction, km
%
% Albino Carboganni, INAF-OAS
% Version Mar 09, 2023


function [strewn_field_X, sigma_X, strewn_field_Y, sigma_Y]=strewn_field_extensions(LAT, LONG, a_R, Lat_impact, Long_impact)

% Position of strewn field points relative to the nominal centre, km
grado_Lat=111.195;
grado_Long=(2*pi*6371*cosd(LAT))/360;
relative_Lat_impact=(Lat_impact-LAT)*grado_Lat;     

relative_Long_impact=(Long_impact-LONG)*grado_Long;

% Rotation of strewn field points to place them in a frame of reference where the
% direction of motion of the fireball is the y axis and the orthogonal direction is
% the x-axis

theta=(a_R-180); % Rotation angle (degrees)

% Counterclockwise rotation matrix

R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];

% Rotation of strewn field points

N=length(Lat_impact);
strewn_field_X=zeros(1, N); strewn_field_Y=zeros(1, N); 

for i=1:N
point = [relative_Long_impact(i) relative_Lat_impact(i)]';
rotpoint = R*point;
strewn_field_X(i)=rotpoint(1);
strewn_field_Y(i)=rotpoint(2);
end

% Compute standar deviation strewn field in orthogonal and parallel
% direction fireball's motion, km
sigma_X=std(strewn_field_X);
sigma_Y=std(strewn_field_Y);

end
