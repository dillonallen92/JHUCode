function [attitude_mat, eigen_axis, theta_deg] =  TRIAD_Algo(b1,b2,i1,i2)

ub = b1/norm(b1);
vb = cross(ub,b2)/norm(cross(ub,b2));
wb = cross(ub,vb)/norm(cross(ub,vb));

ui = i1 / norm(i1);
vi = cross(ui, i2) / norm(cross(ui,i2));
wi = cross(ui, vi) / norm(cross(ui,vi));
attitude_mat = [ub vb wb] * transpose([ui vi wi])

theta = acos(.5*(trace(attitude_mat) - 1));
theta_deg = theta * 180 / pi

% axis of rotation
eigen_axis = (1/(2*sin(theta))) * ...
    [ attitude_mat(2,3) - attitude_mat(3,2);...
        attitude_mat(3,1) - attitude_mat(1,3);...
            attitude_mat(1,2) - attitude_mat(2,1)]

end