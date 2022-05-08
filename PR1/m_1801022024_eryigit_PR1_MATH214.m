clc
clear
clear all
%The Bisection Method
a1 = -3; %endpoint
b1 = 10; %endpoint
tol = 1.e-10; %tolerance
n0=50; %max iteration
n=1;
e0 = 1/36*pi*10^-9;%definition
E= @(x) (1/4*pi*e0)*[(13*(x+7))/(abs(x+7)^3)+(9*(x+4))/(abs(x+4)^3)+(6*(x-11))/(abs(x-11)^3)+(3*(x-14))/(abs(x-14)^3)]; %Function
while n<n0
    p=(a1+b1)/2; %bisection method formula
    if(E(p)==0 ||(b1-a1)/2 < tol) %iteration limit point
        fprintf('OUTPUT(%f) by Bisection Method \n',p); %printing the result to the screen
        break %break of the loop
    end
    root1(n)=p; %root array
    n=n+1; %increase iteration
    if (E(p)*E(b1))< 0
        a1=p; %assigment
    else
        b1=p; %assigment
    end
   if(n==50)
      fprintf('Method failed after %d iterations by Bisection Method \n',n0) %error massage
   end
end
plot(root1,'b--'); %drawing graphics
    xlabel('Number of iteration')
    ylabel('Root convergence')
    hold on

clear
clear all
%The Newton Method
tol = 1.e-10; %tolerance
n0=50; %max iteration
e0 = 1/36*pi*10^-9;%definition
n=1;
p0=2; %initial approximation
E= @(x) (1/4*pi*e0)*[(13*(x+7))/(abs(x+7)^3)+(9*(x+4))/(abs(x+4)^3)+(6*(x-11))/(abs(x-11)^3)+(3*(x-14))/(abs(x-14)^3)]; %Function
dE= @(x) (1/4*pi*e0)*[-(26/(abs(x+7)*(x+7)^2))-(18/(abs(x+4)*(x+4)^2))-(12/(abs(x-11)*(x-11)^2))-(6/(abs(x-14)*(x-14)^2))]; %derivative of the function
while(n<50)
p=p0-E(p0)/dE(p0); %newton method formula
if(abs(p-p0)<tol)%iteration limit point
    fprintf('OUTPUT(%f) by Newton Method\n',p); %printing the result to the screen
    break
end
root2(n)=p; %root array
n=n+1; %increase iteration
p0=p; %assigment
if(n==50)
      fprintf('Method failed after %d iterations by Newton Method\n',n0) %error massage
end

end
plot(root2,'r-'); %drawing graphics
    hold on
    
clear
clear all
%The Secant Method
tol = 1.e-10; %tolerance
n0=50; %max iteration
e0 = 1/36*pi*10^-9;%definition
n=2;
p0=-1; %initial approximation
p1=8; %initial approximation
E= @(x) (1/4*pi*e0)*[(13*(x+7))/(abs(x+7)^3)+(9*(x+4))/(abs(x+4)^3)+(6*(x-11))/(abs(x-11)^3)+(3*(x-14))/(abs(x-14)^3)]; %Function
q0=E(p0); %assigment
q1=E(p1); %assigment
while(n<50)
p=p1-q1*(p1-p0)/(q1-q0); %secant method formula
if(abs(p-p1)<tol) %iteration limit point
    fprintf('OUTPUT(%f) by Secant Method\n',p);
    break
end
n=n+1; %increase iteration
q0=q1; %assigment
p0=p1; %assigment
p1=p; %assigment
q1=E(p); %assigment
root3(n)=p; %root array
if(n==50)
      fprintf('Method failed after %d iterations by Secant Method',n0) %error massage
end

end
plot(root3,'k-.'); %drawing graphics
legend('Bisection','Newton','Secant')