function trueDistance = distanceFromCenterOfCloud(locations)
% The locations are obained from the names of the TDMS files.
% locations should be a n by 2 array. The first column is the
% x coordinates, the second is the y coordinates.


    xPositions = locations(:,1);
    trueXPosition = xPositions - median(xPositions);
    
    yPositions = locations(:,2);
    trueYPosition = yPositions - median(yPositions);
    
    distance = sqrt(trueXPosition.^2 + trueYPosition.^2);
    
    trueDistance = distance * (1/4000) * (1/39.3701);
end