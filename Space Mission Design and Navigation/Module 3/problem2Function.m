function [r,v] = problem2Function(a,e,i,Omega,omega,theta)

    mu = 132712440041.94; % km^3/s^2
    KmInAU = 149597870.7; % AU/ km
    secondsInDay = 24*3600; % days / second
    

    a = a * KmInAU;

    C_313=(C3(omega)*C1(i)*C3(Omega))';

    p=a*(1 - e^2);
    
    % First transform into the perifocal frame

    rPerifocal=p/(1+e*cos(theta))*[cos(theta);sin(theta);0];

    vPerifocal=sqrt(mu/p)*[-sin(theta); e + cos(theta); 0];

    % Perform 3-1-3 Rotation on rPerifocal and vPerifocal
    r=C_313*rPerifocal / KmInAU;

    v=C_313*vPerifocal * secondsInDay / KmInAU;


 
    % Rotation Inner Functions
    function C1=C1(x)

    C1=[1,0,0;0,cos(x),sin(x);0,-sin(x),cos(x)];

    end



    function C3=C3(x)

    C3=[cos(x),sin(x),0;-sin(x),cos(x),0;0,0,1];

    end
    
    
end