function [xDot,t]=twoBody(~,x,mu)
xDot=zeros(6,1);
r=norm(x(1:3,1));
c1=1/(r^3);
xDot(1:3,1)=x(4:6,1);
xDot(4:6,1)=-mu*c1*x(1:3,1);
end