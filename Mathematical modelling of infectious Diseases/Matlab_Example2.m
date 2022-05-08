function EulerSIR()
clear all; clc;
N=1000; mu=1/60;
beta1=20/100;beta2=15/10000;sigma=0.03;gamma=0.1;
options=odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-4]);
[T,Y]=ode45(@SIRmodel,[0 140],[N-1 1 0],options);
plot(T,Y(:,1),'r',T,Y(:,2),'b',T,Y(:,3),'k','Linewidth',2)
xlabel('time')
ylabel('population')
function dy=SIRmodel(t,y)
dy=zeros(3,1);
dy(1)=mu*N-beta1*beta2*y(2)-mu*y(1);
dy(2)=beta1*beta2*y(1)*y(2)-(mu+sigma+gamma)*y(2);
dy(3)=gamma*y(2)-mu*y(3);
end
end

