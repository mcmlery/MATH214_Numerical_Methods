clear all 
clear
clc
load current1.dat
load current2.dat
load current3.dat
load current4.dat
%definition

c1=zeros(9,2);
c2=zeros(13,2);
c3=zeros(25,2);
c4=zeros(61,2);
R=14.2;
L=0.98;

for i=1 :9
for j=1 :2
c1(i,j)=current1(i,j);
end
end

for i=1 :13
for j=1 :2
c2(i,j)=current2(i,j);
end
end

for i=1 :25
for j=1 :2
c3(i,j)=current3(i,j);
end
end

for i=1 :61
for j=1 :2
c4(i,j)=current4(i,j);
end
end

%forward difference
for k=1 :8
    h=c1(2,1);
    f1(k,1)=(c1((k+1),2)-c1(k,2))/h;
    f1(9,1)=0;
end

for k=1 :12
    h=c2(2,1);
    f2(k,1)=(c2((k+1),2)-c2(k,2))/h;
    f2(13,1)=0;
end

for k=1 :24
    h=c3(2,1);
    f3(k,1)=(c3((k+1),2)-c3(k,2))/h;
    f3(25,1)=0;
end

for k=1 :60
    h=c4(2,1);
    f4(k,1)=(c4((k+1),2)-c4(k,2))/h;
    f4(61,1)=0;
end

%backward difference
for k=9:-1:2
    h=c1(2,1);
    b1(k,1)=(c1((k),2)-c1(k-1,2))/h;
    
end

for k=13:-1 :2
    h=c2(2,1);
    b2(k,1)=(c2((k),2)-c2(k-1,2))/h;
end

for k=25:-1 :2
    h=c3(2,1);
    b3(k,1)=(c3((k),2)-c3(k-1,2))/h;
end

for k=61:-1:2
    h=c4(2,1);
    b4(k,1)=(c4((k),2)-c4(k-1,2))/h;
end

%centered difference
for k=2 :8
    h=c1(2,1);
    cd1(1,1)=((-3*c1(1,2))+4*c1(1+1,2)-c1(1+2,2))/2*h; %endpoint formula
    cd1(k,1)=(c1((k+1),2)-c1(k-1,2))/2*h; %midpoint formula
    cd1(9,1)=0;
end

for k=2 :12
    h=c2(2,1);
    cd2(1,1)=((-3*c2(1,2))+4*c2(1+1,2)-c2(1+2,2))/2*h; %endpoint formula
    cd2(k,1)=(c2((k+1),2)-c2(k-1,2))/2*h; %midpoint formula
    cd2(13,1)=0;
end

for k=2 :24
    cd3(1,1)=((-3*c3(1,2))+4*c3(1+1,2)-c3(1+2,2))/2*h; %endpoint formula
    cd3(k,1)=(c3((k+1),2)-c3(k-1,2))/2*h; %midpoint formula
    cd3(25,1)=0;
end

for k=2 :60
    cd4(1,1)=((-3*c4(1,2))+4*c4(1+1,2)-c4(1+2,2))/2*h; %endpoint formula
    cd4(k,1)=(c4((k+1),2)-c4(k-1,2))/2*h; %midpoint formula
    cd4(61,1)=0;
end
%inductance value for forward difference
for e=1 :9
Ef1(e,1)=L*f1(e,1)+R*c1(e,2);
end
for e=1 :13
Ef2(e,1)=L*f2(e,1)+R*c2(e,2);
end
for e=1 :25
Ef3(e,1)=L*f3(e,1)+R*c3(e,2);
end
for e=1 :61
Ef4(e,1)=L*f4(e,1)+R*c4(e,2);
end
%inductance value for backward difference
for e=1 :9
Eb1(e,1)=L*b1(e,1)+R*c1(e,2);
end
for e=1 :13
Eb2(e,1)=L*b2(e,1)+R*c2(e,2);
end
for e=1 :25
Eb3(e,1)=L*b3(e,1)+R*c3(e,2);
end
for e=1 :61
Eb4(e,1)=L*b4(e,1)+R*c4(e,2);
end
%inductance value for centered difference
for e=1 :9
Ec1(e,1)=L*cd1(e,1)+R*c1(e,2);
end
for e=1 :13
Ec2(e,1)=L*cd2(e,1)+R*c2(e,2);
end
for e=1 :25
Ec3(e,1)=L*cd3(e,1)+R*c3(e,2);
end
for e=1 :61
Ec4(e,1)=L*cd4(e,1)+R*c4(e,2);
end

%plot// if you want to view any table use the ctrl+r key combination for comment line 


% plot(f1)
% hold on
% plot(f2)
% hold on
% plot(f3)
% hold on
% plot(f4), title('Derivative table for forward difference') ,legend('current_1', 'current_2', 'current_3', 'current_4');
% hold on
% 
% plot(b1)
% hold on
% plot(b2)
% hold on
% plot(b3)
% hold on
% plot(b4), title('Derivative table for backward difference') ,legend('current_1', 'current_2', 'current_3', 'current_4');
% hold on
% 
% plot(cd1)
% hold on
% plot(cd2)
% hold on
% plot(cd3)
% hold on
% plot(cd4), title('Derivative table for centered difference') ,legend('current_1', 'current_2', 'current_3', 'current_4');
% hold on
% 
% plot(Ef1)
% hold on
% plot(Ef2)
% hold on
% plot(Ef3)
% hold on
% plot(Ef4), title('inductance value table for forward difference') ,legend('current_1', 'current_2', 'current_3', 'current_4');
% hold on
% 
% plot(Eb1)
% hold on
% plot(Eb2)
% hold on
% plot(Eb3)
% hold on
% plot(Eb4), title('inductance value table for backward difference') ,legend('current_1', 'current_2', 'current_3', 'current_4');
% hold on
% 
% plot(Ec1)
% hold on
% plot(Ec2)
% hold on
% plot(Ec3)
% hold on
% plot(Ec4), title('inductance value table for centered difference') ,legend('current_1', 'current_2', 'current_3', 'current_4');
% hold on

