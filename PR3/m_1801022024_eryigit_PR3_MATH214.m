clear all
clear 
clc

load pr3data.dat

t=pr3data(:,1)'; i=pr3data(:,2)'; v=pr3data(:,3)';
%Current graph
figure
plot(t,i);
xlabel('Time[ms]')
ylabel('Current[A]')
title('Current/Time graph');
%voltage graph
figure
plot(t,v);
xlabel('Time[ms]')
ylabel('Voltage[V]')
title('Voltage/Time graph');
%Derivative of Current
for k=2 : length(t)-1 
    h=pr3data(2,1);
    deriv_c(1,1)=(((-3*i(1,1))+4*i(1,2)-i(1,3))/(2*h)); 
    deriv_c(1,41)=(((3*i(1,41))-4*i(1,40)+i(1,39))/(2*h));
    deriv_c(1,k)=((i(1,(k+1))-i(1,k-1))/(2*h));
end
figure
plot(t,deriv_c)
xlabel('Time[ms]')
ylabel('Current[A]');
title('Numerical differentiation of voltage')
%Power graph
for k=1 :length(i) 
    p(:,k)=i(:,k)*v(:,k);
end
figure
plot(t,p)
xlabel('Time[ms]')
ylabel('Power[W]')
title('Power/Time graph');
%definition
a=0; %lower limit 
b=1; %upper limit
L=0.1;

%trapezoidal rule
n=length(t)-1; 
h=(b-a)/n;
cxc=0;
for k=2:40
    cxc=cxc+i(1,k)*v(1,k);
end
pptrapp=(h/2)*(i(1,1)*v(1,1)+(2*cxc)+i(1,41)*v(1,41));

%simpsons rule
n=80;
h=(b-a)/n;
cx1=0;
cx2=0;
for j=3:2:n/2-1
    cx1=cx1+i(1,j)*v(1,j);
end
for k=2:2:n/2
    cx2=cx2+i(1,k)*v(1,k);
end
ppsimpp=(2*cx1+4*cx2+i(1,1)*v(1,1)+v(1,41)*i(1,41))*h/3;

%midpoint rule
n=80+2;
h=(b-a)/n;
cvc=0;
for j=1:1:n/2
    cvc=cvc+i(1,j)*v(1,j);
end
ppmidpp=2*h*cvc;
%work
w=1/2*L*(i(1,41)^2);