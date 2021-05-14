function C=fromInToFixed(nB,t)
theta=nB*t;
C=[cos(theta),sin(theta),0;-sin(theta),cos(theta),0;0,0,1];
end
