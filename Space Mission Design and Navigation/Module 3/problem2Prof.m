function [rV,vV] = problem2Prof(a,e,i,Omega,omega,theta)

    mu = 132712440041.94; % km^3/s^2
    KmInAU = 149597870.7; % AU/ km
    secondsInDay = 24*3600; % days / second
    
    % Using Prussing as a guide
    a = a * KmInAU;

    C=(C3(omega)*C1(i)*C3(Omega))';

    p=a*(1 - e^2);

    rVperifocal=p/(1+e*cos(theta))*[cos(theta);sin(theta);0];

    vVperifocal=sqrt(mu/p)*[-sin(theta); e + cos(theta); 0];

    rV=C*rVperifocal / KmInAU;

    vV=C*vVperifocal * secondsInDay / KmInAU;


 

    function rot=C1(x)

    rot=[1,0,0;0,cos(x),sin(x);0,-sin(x),cos(x)];

    end



    function rot=C3(x)

    rot=[cos(x),sin(x),0;-sin(x),cos(x),0;0,0,1];

    end


end