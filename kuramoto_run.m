
%Keep the file frekon.m in the sam folder

%clear all; close all; clc;
nn=10;


for i=1:nn
    theta(i) = pi*(-1 + 2*(i-1)./(nn-1));
end



[t,x]=ode45('kuramoto',[0 500],theta);

plot(t,x)