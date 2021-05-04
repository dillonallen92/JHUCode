% Frame rotation matrices:
function [R] = rot(xyzOpt, ang_Rad)
Ca = cos(ang_Rad);
Sa = sin(ang_Rad);
if (xyzOpt==1) % Rot about x
    R = [1,0,0;
        0,Ca,Sa;
        0,-Sa,Ca];
elseif (xyzOpt==2) % Rot about y
    R = [Ca,0,-Sa;
        0,1,0;
        Sa,0,Ca];
elseif (xyzOpt==3) %Rot about z
    R = [Ca,Sa,0;
        -Sa,Ca,0;
        0,0,1];
else
    error('Rotation options are only in 1-3')
end