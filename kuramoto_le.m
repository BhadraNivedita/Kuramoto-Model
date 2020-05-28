function f=kuramoto_le(t,X);

%  Kuramoto equation 
%
%  theta_i_dot = omega_i + K/N * SUM_{j=1}^N sin(theta_j-theta_i)

% Values of parameters

nn=10; kk=0.5;

for i=1:nn
    w(i) = -1 + 2*(i-1)/(nn-1);
end

for i=1:nn
    for j=1:nn
    Y(i,j) = X(nn*j+i);
    end
end

f=zeros(nn^2,1);

% Flow equation

for i=1:nn
	f(i) = w(i);
	for j=1:nn
	f(i) = f(i) + (kk/nn)*sin(X(j)-X(i));
    end
end

% Linearized system

Jac = 0.0;

for i=1:nn
for j=1:nn
Jac(i,j) = cos(X(j)-X(i));
end
end

for i=1:nn
for j=1:nn
Jac(i,i) = Jac(i,i) - cos(X(j)-X(i));
end
end

Jac = (kk/nn)*Jac;
  
% Variational equation

f(nn+1:nn+nn^2)=Jac*Y;

% To run: "kuramoto_le_run.m"