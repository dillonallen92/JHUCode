% LLA2ECEF - convert latitude, longitude, and altitude to
%            earth-centered, earth-fixed (ECEF) cartesian
% 
% USAGE:
% [x,y,z] = lla2ecef(lat,lon,alt)
% 
% x = ECEF X-coordinate (m)
% y = ECEF Y-coordinate (m)
% z = ECEF Z-coordinate (m)
% lat = geodetic latitude (radians)
% lon = longitude (radians)
% alt = height above WGS84 ellipsoid (m)
% 
% Notes: This function assumes the WGS84 model.
%        Latitude is customary geodetic (not geocentric).
% 
% Source: "Department of Defense World Geodetic System 1984"
%         Page 4-4
%         National Imagery and Mapping Agency
%         Last updated June, 2004
%         NIMA TR8350.2
% 
% Michael Kleder, July 2005
% Modified by Dillon Allen to match notation in TSB pg. 79

function [x,y,z] = lla2ecef(lat,lon,alt)

% WGS84 ellipsoid constants:
Re = 6378137;
e = 8.1819190842622e-2;

% intermediate calculation
% (prime vertical radius of curvature)
Nh = Re ./ sqrt(1 - e^2 .* sin(lat).^2);

% results:
x = (Nh+alt) .* cos(lat) .* cos(lon);
y = (Nh+alt) .* cos(lat) .* sin(lon);
z = ((1-e^2) .* Nh + alt) .* sin(lat);

end