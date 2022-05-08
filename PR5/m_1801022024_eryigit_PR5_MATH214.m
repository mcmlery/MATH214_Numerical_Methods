clc
clear 

load pr5data.dat;
xi=pr5data(:,1);  yi=pr5data(:,2);
%Linear least squares
m=length(xi);
txi=0;
tyi=0;
txi2=0;
txiyi=0;
plsa0=0;
plsa1=0;
for i=1 : m
    txi=xi(i,1)+txi;
    tyi=yi(i,1)+tyi;
    txi2=(xi(i,1)^2)+txi2;
    txiyi=(xi(i,1)*yi(i,1))+txiyi;
end
plsa0=(txi2*tyi-txiyi*txi)/(m*txi2-(txi^2));
plsa1=(m*txiyi-txi*tyi)/(m*txi2-(txi^2));

plsPx = @ (x) plsa1*x+plsa0;
plsy=zeros(m,1);
for i=1 : m
   plsy(i,1)=plsPx(xi(i,1));
end
result1=zeros(m,1);
for i=1 : m
    result1(i,1)=yi(i,1)-plsy(i,1);
end
plot(xi,plsy);
hold on
plot(xi,yi,'o');
xlabel('Distance[m]');
ylabel('Voltage[mV]');
title('Linear least squares');
legend('Aproximation value','Reel value');

%polinomial least squares degree 2
txin=0;
txi2n=0;
txi3n=0;
txi4n=0;
txiyi2n=zeros(3,1);
plsa2=zeros(3,1);
n=0;
n1=2;
for i=1 : m
    txin=xi(i,1)+txin;
    txi2n=(xi(i,1)^2)+txi2n;
    txi3n=(xi(i,1)^3)+txi3n;
    txi4n=(xi(i,1)^4)+txi4n;
end
for n=0 : 2
for i=1 : m
    txiyi2n(n+1,1)=((xi(i,1)^n)*yi(i,1))+txiyi2n(n+1,1);
end
end
                         %equation system definition
%     plsa2(1,1)*m+plsa2(2,1)*txin+plsa2(3,1)*txi2n==txiyi2n(1,1);
%     plsa2(1,1)*txin+plsa2(2,1)*txi2n+plsa2(3,1)*txi3n==txiyi2n(2,1);
%     plsa2(1,1)*txi2n+plsa2(2,1)*txi3n+plsa2(3,1)*txi4n==txiyi2n(3,1);

lsd2_1=[m txin txi2n; txin txi2n txi3n; txi2n txi3n txi4n];
lsd2_2=[plsa2(1,1);plsa2(2,1);plsa2(3,1)];
lsd2_3=[txiyi2n(1,1);txiyi2n(2,1);txiyi2n(3,1)];

invlsd2=inv(lsd2_1);
lsd2_2=invlsd2*lsd2_3;

plsPx2 = @ (x) lsd2_2(3,1)*x^2+lsd2_2(2,1)*x+lsd2_2(1,1);
plsy2=zeros(m,1);
for i=1 : m
   plsy2(i,1)=plsPx2(xi(i,1));
end
result2=zeros(m,1);
for i=1 : m
    result2(i,1)=yi(i,1)-plsy2(i,1);
end
figure
plot(xi,plsy2);
hold on
plot(xi,yi,'o');
xlabel('Distance[m]');
ylabel('Voltage[mV]');
title('Polinomial least squares n=2');
legend('Aproximation value','Reel value');

%polinomial least squares degree 3
txin=0;
txi2n=0;
txi3n=0;
txi4n=0;
txi5n=0;
txi6n=0;
txiyi3n=zeros(4,1);
plsa3=zeros(4,1);
n=0;
n2=3;
for i=1 : m
    txin=xi(i,1)+txin;
    txi2n=(xi(i,1)^2)+txi2n;
    txi3n=(xi(i,1)^3)+txi3n;
    txi4n=(xi(i,1)^4)+txi4n;
    txi5n=(xi(i,1)^5)+txi5n;
    txi6n=(xi(i,1)^6)+txi6n;
end
for n=0 : 3
for i=1 : m
        txiyi3n(n+1,1)=((xi(i,1)^n)*yi(i,1))+txiyi3n(n+1,1);
end
end
                                %equation system definition
%     plsa3(1,1)*m+plsa3(2,1)*txin+plsa3(3,1)*txi2n+plsa3(4,1)*txi3n==txiyi3n(1,1);
%     plsa3(1,1)*txin+plsa3(2,1)*txi2n+plsa3(3,1)*txi3n+plsa3(4,1)*txi4n==txiyi3n(2,1);
%     plsa3(1,1)*txi2n+plsa3(2,1)*txi3n+plsa3(3,1)*txi4n+plsa3(4,1)*txi5n==txiyi3n(3,1);
%     plsa3(1,1)*txi3n+plsa3(2,1)*txi4n+plsa3(3,1)*txi5n+plsa3(4,1)*txi6n==txiyi3n(4,1);

lsd3_1=[m txin txi2n txi3n; txin txi2n txi3n txi4n; txi2n txi3n txi4n txi5n ; txi3n txi4n txi5n txi6n];
lsd3_2=[plsa3(1,1);plsa3(2,1);plsa3(3,1);plsa3(4,1)];
lsd3_3=[txiyi3n(1,1);txiyi3n(2,1);txiyi3n(3,1);txiyi3n(4,1)];

invlsd3=inv(lsd3_1);
lsd3_2=invlsd3*lsd3_3;

plsPx3 = @ (x) lsd3_2(4,1)*x^3+lsd3_2(3,1)*x^2+lsd3_2(2,1)*x+lsd3_2(1,1);
plsy3=zeros(m,1);
for i=1 : m
   plsy3(i,1)=plsPx3(xi(i,1));
end
result3=zeros(m,1);
for i=1 : m
    result3(i,1)=yi(i,1)-plsy3(i,1);
end
figure
plot(xi,plsy3);
hold on
plot(xi,yi,'o');
xlabel('Distance[m]');
ylabel('Voltage[mV]');
title('Polinomial least squares n=3');
legend('Aproximation value','Reel value');