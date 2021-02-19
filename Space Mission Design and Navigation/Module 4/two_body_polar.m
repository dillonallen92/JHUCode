function dx = two_body_polar(~,x,mu_in)

dx=zeros(4,1);

dx(1)=x(3);
dx(2)=x(4);
dx(3)=x(1)*x(4)^2-mu_in/(x(1)^2);
dx(4)=-2*x(3)*x(4)/x(1);

end