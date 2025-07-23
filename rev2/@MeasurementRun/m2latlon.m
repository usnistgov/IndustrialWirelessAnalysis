function [ latO, lonO ] = m2latlon( lat, dn, de )
%M2LATLON Convert meters to lat-lon offset
%   distance offset in meters must be under 1 km
%
%   lat: origin latitude in degrees
%   dn: distance north from origin in meters
%   de: distance east from origin in meters
%
% Author: Rick Candell
% National Institute of Standards and Technology
% Source: http://gis.stackexchange.com/questions/2951/algorithm-for-offsetting-a-latitude-longitude-by-some-amount-of-meters

 %Earth’s radius, sphere
 R=6378137;
 lon = 0;


 %Coordinate offsets in radians
 dLat = dn/R;
 dLon = de/(R*cos(pi*lat/180));

 %OffsetPosition, decimal degrees
 latO = lat + dLat * 180/pi;
 lonO = lon + dLon * 180/pi; 

end

