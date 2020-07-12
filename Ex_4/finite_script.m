N=50; h=1/N; c=1; k=1/10000; v=k/(h^2); B=1-2*v;
u=zeros(N-1,N+1); % defining a 2D array of zeros
for j=0:h:1 % initialize the array, u, with the desired input %function
%To plot one, comment out other two u(1,c) lines in the loop
% u(1,c)=(sin(2*pi*j)); % initializing u for the t=1 case 
% u(1,c)=abs(sin(2*pi*j));
% u(1,c)=cos((pi*j)/2); % to test for non-zero boundaries
%  u(1,c)=(j-0.5)^2; %additional testing case with polynomial term
u(1,c)=exp(j+1)-5;
c=c+1; % c increments as j increments from 0 to 1 in intervals of %0.02
end
%u=triangularPulse(0,1,x); %for tent function, comment out the above % for loop
 
for m=1:20 %m controls loop iterations, ie number of %lines plotted
    for j=2:N %runs N-1 times
% Implementing %central algorithm using the 2D array:
        u(m+1,j)=v*u(m,j-1)+B*u(m,j)+v*u(m,j+1); 
% TESTING FOR DIFFERENT BOUNDARY CONDITIONS: please comment out the % other 2 cases
% Zero Boundary Conditions:
%            u(m+1,1)=0; 
%            u(m+1,N+1)=0;
% Constant non-zero boundary conditions:
%         u(m+1,1)=0.5; 
%         u(m+1,N+1)=0.5;
% Time Varying Non-Zero boundary conditions:
%          u(m+1,1)=cos((pi*m*k)/4); 
%          u(m+1,N+1)=cos((pi*m*k)/4); 
    end
end
x=[0:h:1]; 
 for i=1:m
    hold on
plot(x,u(i,:))
title('Variation of heat distribution across a metal rod with respect to time');
ylabel ('Temperature/K'); 
xlabel ('Length of the rod/m'); 
 end
