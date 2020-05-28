
%Keep the file frekon.m in the sam folder

clear all; close all; clc;
nn=10;


for i=1:nn
    theta(i) = pi*(-1 + 2*(i-1)/(nn-1));
end

[T,Res]=lyapunov(nn,@kuramoto_le,@ode45,0,0.1,1000,theta);
plot(T,Res);