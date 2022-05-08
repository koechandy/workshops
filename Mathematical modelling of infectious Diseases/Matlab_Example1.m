function EulerSIR ()
clear all;clc;
a=0; b=40; nt=4000;
h=(b-a)/nt; t(1)=a;
N=1000; w2(1)=1;w3(1)=0;w1(1)=N-w2(1)-w3(1);
mu=1/60;beta1=20/100;beta2=15/100;sigma=0.03;gamma=0.1;
for i=1:nt
    t(i+1)=t(1)+i*h;
    w1(i+1)=w1(i)+h*fxnS(t(i),w1(i),w2(i),w3(i));
     w2(i+1)=w2(i)+h*fxnI(t(i),w1(i),w2(i),w3(i));
      w3(i+1)=w3(i)+h*fxnR(t(i),w1(i),w2(i),w3(i));
end
plot(t,w1,'b',t,w2,'r',t,w3,'g','linewidth',2)
xlabel('time');
ylabel('population');
function fS=fxnS(t,S,I,R)
fS=mu*N-beta1*beta2*S*I-mu*S;
end
function fI=fxnI(t,S,I,R)
    fI=beta1*beta2*S*I-(mu+sigma+gamma)*I;
end
function fR=fxnR(t,S,I,R)
    fR=gamma*I-mu*R;
end 
end
