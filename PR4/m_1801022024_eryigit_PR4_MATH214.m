clear
clear all
clc
%definition
a=0; b=0.6;
h1=0.05; h2=0.025;
n1=(b-a)/h1; n2=(b-a)/h2;
L=0.98; Vs=12; R=14.2;
w=zeros(12,1);
w(1,1)=0.01; t1= a:h1:b; t2= a:h2:b;
%euler method
for k=1: length(w)
    w(k+1)=w(k)+h1*((Vs-R*w(k))/L);
end
figure
plot(t1,w,'r-');
hold on
w=zeros(12,1);
w(1,1)=0.01;
%midpoint method
for k=1: length(w)
    w(k+1)=w(k)+h1*((h1*k*((Vs-R*w(k))/L))*(Vs-R*w(k))/L);
end
plot(t1,w,'k-');
hold on
w=zeros(12,1);
w(1,1)=0.01;
%modified euler method
for k=1: length(w)
    w(k+1)=w(k)+h1/2*(((Vs-R*w(k))/L)+(((h1*k+1)+h1)*(Vs-R*w(k))/L));
end
plot(t1,w,'g-');
hold on
w=zeros(12,1);
w(1,1)=0.01;
%runge-kutta method order four
for k=1: length(w)
    k1=h1*((Vs-R*w(k))/L);
    k2=h1*((h1*k*((Vs-R*w(k))/L))+(k1/2));
    k3=h1*((h1*k*((Vs-R*w(k))/L))+(k2/2));
    k4=h1*(h1*(k+1)*((Vs-R*w(k))/L)+k3);
    w(k+1)=w(k)+1/6*(k1+2*k2+2*k3+k4);
end
plot(t1,w,'b-');
hold on
xlabel('Time[ms]')
ylabel('Current[A]')
title('Current/Time graph delta(t)=0.05');
legend('Euler','Midpoint','Modified Euler','Runge Kutta');
w=zeros(24,1);
w(1,1)=0.01;

%euler method
for k=1: length(w)
    w(k+1)=w(k)+h2*((Vs-R*w(k))/L);
end
figure
plot(t2,w,'r-');
hold on
w=zeros(24,1);
w(1,1)=0.01;
%midpoint method
for k=1: length(w)
    w(k+1)=w(k)+h2*((h2*k*((Vs-R*w(k))/L))*(Vs-R*w(k))/L);
end
plot(t2,w,'k-');
hold on
w=zeros(24,1);
w(1,1)=0.01;
%modified euler method
for k=1: length(w)
    w(k+1)=w(k)+h2/2*(((Vs-R*w(k))/L)+(((h2*k+1)+h2)*(Vs-R*w(k))/L));
end
plot(t2,w,'g-');
hold on
w=zeros(24,1);
w(1,1)=0.01;
%runge-kutta method order four
for k=1: length(w)
    k1=h2*((Vs-R*w(k))/L);
    k2=h2*((h2*k*((Vs-R*w(k))/L))+(k1/2));
    k3=h2*((h2*k*((Vs-R*w(k))/L))+(k2/2));
    k4=h2*(h2*(k+1)*((Vs-R*w(k))/L)+k3);
    w(k+1)=w(k)+1/6*(k1+2*k2+2*k3+k4);
end
plot(t2,w,'b-');
hold on
xlabel('Time[ms]')
ylabel('Current[A]')
title('Current/Time graph delta(t)=0.025');
legend('Euler','Midpoint','Modified Euler','Runge Kutta');