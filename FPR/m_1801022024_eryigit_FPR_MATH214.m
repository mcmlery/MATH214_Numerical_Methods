clear 
clc
%definition values
load fprdata.dat
v(:,1)=fprdata(:,1); a(:,1)=fprdata(:,2);
lnv=zeros(5,1);
lna=zeros(5,1);
L=0.98 ; Vs=2; R=14.2; n1=600/25; h1=25/1000; n2=600/2.5; h2=2.5/1000;
for i=1 : length(v)
    lnv(i,1)=log(v(i,1));
    lna(i,1)=log(a(i,1));
end
%least square
topa=0; topv=0; toplna=0; topv2=0; topvlna=0;
for i=1 : length(v)
    topa=a(i,1)+topa;
    topv=v(i,1)+topv;
    toplna=lna(i,1)+toplna;
    topv2=v(i,1)^2+topv2;
    topvlna=v(i,1)*lna(i,1)+topvlna;
    
end
kata=(length(v)*topvlna-topv*toplna)/(length(v)*topv2-topv^2);
katlnb=(topv2*toplna-topvlna*topv)/(length(v)*topv2-topv^2);
katb=exp(katlnb);

%functions
id = @ (ger)  katb*exp(kata*ger);
vd = @ (aki) (log(aki)-log(katb))/kata;

%euler method
w1=zeros(n1+1,1);
w2=zeros(n2+1,1);
tvd1=zeros(n1+1,1); tvd2=zeros(n2+1,1);
tvd1(1,1)=0; tvd2(1,1)=0;
tid1=zeros(n1+1,1); tid2=zeros(n2+1,1);
tid2(1,1)=0; tid2(1,1)=0;
w1(1,1)=h1*Vs/L;
w2(1,1)=h2*Vs/L;
san1=0 : h1 : 0.6;
san2=0 : h2 : 0.6;
for i=1 : n1
   w1(i+1,1)= w1(i,1)+h1*((Vs-vd(w1(i,1))-w1(i,1)*R)/L);
   tvd1(i+1,1)=vd(w1(i,1));
   tid1(i+1,1)=w1(i+1,1);
end
for i=1 : n2
   w2(i+1,1)= w2(i,1)+h2*((Vs-vd(w2(i,1))-w2(i,1)*R)/L);
   tvd2(i+1,1)=vd(w2(i,1));
   tid2(i+1,1)=w2(i+1,1);
end
say1=0 :h1 : 1.3;
say2=0 :h2 : 1.3;
%graphs
figure
plot(san1,tvd1)
hold on
plot(san2,tvd2)
xlabel('Time[ms]')
ylabel('Diode Voltage[V]')
title('Diode Voltage/Time graph');
legend('t=0.025 s','t=0.0025s')
figure
plot(san1,tid1)
hold on
plot(san2,tid2)
xlabel('Time[ms]')
ylabel('Current[A]')
title('Main Current/Time graph');
legend('t=0.025 s','t=0.0025 s')
figure 
plot (san1,Vs-tvd1-tid1*R)
hold on
plot (san2,Vs-tvd2-tid2*R)
xlabel('Time[ms]')
ylabel('Inductance Voltage[V]')
title('Inductance Voltage/Time graph');
legend('t=0.025 s','t=0.0025 s')
figure
plot(san1,tid1*R)
hold on
plot(san2,tid2*R)
xlabel('Time[ms]')
ylabel('Voltage[V]')
title('Resistor Voltage/Time graph ');
legend('t=0.025 s','t=0.0025 s')
figure
plot(say1,id(say1))
xlabel('Time[ms]')
ylabel('Voltage[V]')
title('Fitting graph');
hold on
plot(v,a,'o');
legend('Least square','fprdata')