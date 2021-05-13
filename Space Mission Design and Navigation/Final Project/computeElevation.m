function [elevation_deg] = computeElevation(rVec,r_IVec)
    r_tVec = rVec - r_IVec;
    % Compute azimuth and elevation
    r_t    = norm(r_tVec);
    % r_t_xy = norm(r_tVec(1:2));
    %
    arg_el = r_tVec(3)/r_t;
    %
    % az_rad = asin(arg_az);% not complete
    % az_rad = atan2(r_tVec(1),r_tVec(2));
    elevation_deg = asin(arg_el) * 180/pi;
end