% Dylan Jackson
% ME 3165-308
% 2/12/2018
% Thomas Algorithm
clear;clc;close all
% Domain
h=.25;
n=100;
dy=h/n;
U=zeros(n+1,5);
rho=1000;
L=.1;

% Boundary conditions
U(1,:)=0;
U(n+1,:)=.01;

% Pressure Gradient
dpdx=[-1e-3,-.5e-3,0,.5e-3,1e-3];

vis=8.9e-4;
b=ones(n,5)*-2;

d=ones(n,length(dpdx))*dy^2*(1/vis).*dpdx;
a=ones(n,5)*1;
c=ones(n,5)*1;
 for g=1:5
    
d(n,g)=d(n,g)-U(n+1,g);
for k = 2:n
m = a(k,g)/b(k-1,g); 
b(k,g) = b(k,g) - m*c(k-1,g);  
d(k,g) = d(k,g) - m*d(k-1,g);  
end 
% Backward substitution phase 
U(n,g) = d(n,g) / b(n,g); 
for k = n-1:-1:2 
U(k,g) = (d(k,g) - c(k,g)*U(k+1,g)) / b(k,g);  
end
 end
figure
z=linspace(0,h,n+1);
plot(U,z)
title('Velocity Profile in pipe')
xlabel('Velocity in m/s')
ylabel('Height in Pipe')
legend('dpdx=-1e-3','dpdx=-.5e-3','dpdx=0','dpdx=.5e-3','dpdx=1e-3')

Rey=rho*L/vis.*U;

