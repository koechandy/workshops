function Malaria()
clear all; clc;
H=100;V=1000;a=1/3;bv=0.05;muv=0.1;bh=0.2;r=1/50 ;
options=odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-4 1e-4]);
[T,Y]=ode45(@MALARIAmodel,[0 100],[H-1 1 V 0],options);
plot(T,Y(:,1),'r',T,Y(:,2),'b',T,Y(:,3),'k',T,Y(:,4),'g','Linewidth',2)
xlabel('time')
ylabel('population')
function dy=MALARIAmodel(t,y)
dy=zeros(4,1);
dy(1)=-(V/H)*a*bh*(y(4)/V)*y(1)+r*y(2);
dy(2)=(V/H)*a*bh*(y(4)/V)*y(1)-r*y(2);
dy(3)=muv*V-a*bv*(y(1)/H)*y(3)-muv*y(3);
dy(4)=a*bv*(y(2)/H)*y(3)-muv*y(4);
end
end
