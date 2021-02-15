function [a_au,e_norm,i,omega,Omega,theta] = problem1Function(r,v)

    % Constants
    mu = 132712440041.94; % km^3/s^2
    auToKm = 149597870.7; % km
    dayToSeconds = 1/(24*3600); % 1/day to 1/s
    I_hat = [1; 0; 0];
    K_hat = [0; 0; 1];
    
    % first convert r and v into desired units
    r_km = r.*auToKm; % km
    v_kms = v.*(auToKm)*(dayToSeconds); % km/s
    
    % Following the algorithm in lecture
    E = dot(v_kms,v_kms)/2 - mu/norm(r_km); %km^2/s^2
    a = -mu/(2*E); % km
    a_au = a / auToKm; % AU
    
    % now find h_vec and e_vec and e_norm
    h_vec = cross(r_km,v_kms); % km^2/s
    e_vec = cross(v_kms,h_vec)/mu - r_km/norm(r_km); % unitless
    e_norm = norm(e_vec); % unitless
    
    % Now find the inclination
    i = acos(dot(h_vec,K_hat)/norm(h_vec)); % rads
    
    % Now find the right ascension of the ascending node Omega
    N_hat = cross(K_hat,h_vec)/norm(cross(K_hat,h_vec));
    Omega = acos(dot(N_hat, I_hat)); % radians
    if N_hat(2) < 0
        Omega = 2*pi - Omega; % radians
    end
    
    % Calculate the argument of periapsis (omega)
    omega = acos(dot(N_hat, e_vec)/e_norm); % radians
    if e_vec(3) < 0
        omega = 2*pi - omega;
    end
    
    % Finally, find the true anomaly (theta)
    theta = acos(dot(e_vec, r_km)/(e_norm * norm(r_km)));
    if dot(r_km, v_kms) < 0
        theta = 2*pi - theta;
    end
    
end