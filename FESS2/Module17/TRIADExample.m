clear, clc, close all;

%%
b1 = [1;0;0];
b2 = [0;0;1];
i1 = [sqrt(2)/2; sqrt(2)/2; 0];
i2 = [0;0;1];

ub = b1/norm(b1);
vb = cross(ub,b2)/norm(cross(ub,b2));
wb = cross(ub,vb)/norm(cross(ub,vb));

ui = i1 / norm(i1);
vi = cross(ui, i2) / norm(cross(ui,i2));
wi = cross(ui, vi) / norm(cross(ui,vi));
A = [ub vb wb] * transpose([ui vi wi]);

theta = acos(.5*(trace(A) - 1));
theta_deg = theta * 180 / pi;

% axis of rotation
v_rot = (1/(2*sin(theta))) * ...
    [ A(2,3) - A(3,2), A(3,1) - A(1,3), A(1,2) - A(2,1)];

%%

b1 = [1;0;0];
b2 = [0;0;1];
i1 = [sqrt(2)/2; sqrt(2)/2; 0];
i2 = [0;0;1];

[attitude_mat, eigen_axis, theta_deg] = TRIAD_Algo(b1,b2,i1,i2)
