%Exercise 3 - RLC Series Circuit 

R = 250; %resistance in Ohms 
L = 0.650; %inductance in henries 
C = 3*10^-6; %capacitance in farads
qi = 500*10^(-9); %initial value of the charge 
i0 = 0; %inital value of the current 
ti = 0; %t = 0; (initial time)

figure; 

vin = @(t) 5; %step input signal 
h = 0.00001; %step size 
tfinal = 0.05; 
N = round((tfinal-ti)/h); %to set the upper bound of the loop
Ya = zeros(1,N); %array for current with zeros 
Xb = zeros(1,N); %array for charge (state of the equation) with zeros 
Tc = zeros(1,N); %array for time with zeros 

Xb(1) = qi; %first element of the array equal to inital charge 
Ya(1) = i0; %first element of the array equal to inital current 
Tc(1) = ti; %first element of the array equal to initial time 

for i=1:N-1; 
    [Ya(i+1), Xb(i+1)] = RK4second(R,C,L,h,Xb(i), Ya(i), Tc(i),vin); 
    %calling the 3/8 runge kutta 4th order method 
    Tc(i+1)=Tc(i)+h; %increment the time 
end

subplot(3,3,1) 
plot(Tc,R*Ya); %Graph between the output voltage vs time
title('Vout with step signal input with 5 V amplitude');
ylabel ('Vout/V'); 
xlabel ('t/s'); 

%Impulsive signal with decay with Vin = 5V and tau = 3 (ms)^2
vin = @(t) 5 * exp(-t^2/(3*10^-6)); 

for i=1:N-1; 
    [Ya(i+1), Xb(i+1)] = RK4second(R,C,L,h,Xb(i), Ya(i), Tc(i),vin); 
    Tc(i+1)=Tc(i)+h; 
end

subplot(3,3,2); 
plot(Tc,R*Ya);  
title('Vout with a 5 V input impulsive signal with decay');
ylabel ('Vout/V'); 
xlabel ('t/s'); 

%square wave with amplitude Vin = 5V for different frequencies: 
%f = 5 Hz, 100 Hz, 500 Hz 

%f=5 Hz

vin = @(t) 5*square(10*pi*t);
h = 0.0001;
tfinal = 0.5; 
N = round((tfinal-ti)/h); 
Ya = zeros(1,N); 
Xb = zeros(1,N); 
Tc = zeros(1,N); 

Xb(1) = qi; 
Ya(1) = i0; 
Tc(1) = ti; 

for i=1:N-1 
    [Ya(i+1), Xb(i+1)] = RK4second(R,C,L,h,Xb(i), Ya(i), Tc(i),vin); 
    Tc(i+1)=Tc(i)+h; 
end

subplot(3,3,3); 
plot(Tc,R*Ya); 
title('Vout with square wave input with f = 5 Hz and amplitude 5 V');
ylabel ('Vout/V'); 
xlabel ('t/s'); 

Vin = vin(Tc); 
hold on; 
plot(Tc, Vin); 

%f=100 Hz

vin = @(t) 5*square(200*pi*t);
h = 0.00001;
tfinal = 0.05; 
N = round((tfinal-ti)/h); 
Ya = zeros(1,N); 
Xb = zeros(1,N); 
Tc = zeros(1,N); 

Xb(1) = qi; 
Ya(1) = i0; 
Tc(1) = ti; 

for i=1:N-1 
    [Ya(i+1), Xb(i+1)] = RK4second(R,C,L,h,Xb(i), Ya(i), Tc(i),vin); 
    Tc(i+1)=Tc(i)+h; 
end

subplot(3,3,4); 
plot(Tc,R*Ya); 
title('Vout with square wave input with f = 100 Hz and amplitude 5 V');
ylabel ('Vout/V'); 
xlabel ('t/s'); 

Vin = vin(Tc); 
hold on; 
plot(Tc, Vin); 

%f=500 Hz

vin = @(t) 5*square(1000*pi*t);
h = 0.00001;
tfinal = 0.01; 
N = round((tfinal-ti)/h); 
Ya = zeros(1,N); 
Xb = zeros(1,N); 
Tc = zeros(1,N); 

Xb(1) = qi; 
Ya(1) = i0; 
Tc(1) = ti; 

for i=1:N-1 
    [Ya(i+1), Xb(i+1)] = RK4second(R,C,L,h, Xb(i), Ya(i), Tc(i),vin); 
    Tc(i+1)=Tc(i)+h; 
end

subplot(3,3,5); 
plot(Tc,R*Ya); 
title('Vout with square wave input with f = 500 Hz and amplitude 5 V');
ylabel ('Vout/V'); 
xlabel ('t/s'); 

Vin = vin(Tc); 
hold on; 
plot(Tc, Vin); 

%sine wave with amplitude Vin = 5V for different frequencies: 
%f = 5 Hz, 100 Hz, 500 Hz 

%f=5 Hz


vin = @(t) 5 * sin(10*pi*t);
h = 0.0001;
tfinal = 0.5; 
N = round((tfinal-ti)/h); 
Ya = zeros(1,N); 
Xb = zeros(1,N); 
Tc = zeros(1,N); 

Xb(1) = qi; 
Ya(1) = i0; 
Tc(1) = ti; 

for i=1:N-1 
    [Ya(i+1), Xb(i+1)] = RK4second(R,C,L,h,Xb(i), Ya(i), Tc(i),vin); 
    Tc(i+1)=Tc(i)+h; 
end

subplot(3,3,6); 
plot(Tc,R*Ya); 
title('Vout with sine wave input with f = 5 Hz and amplitude 5 V');
ylabel ('Vout/V'); 
xlabel ('t/s'); 

Vin = vin(Tc); 
hold on; 
plot(Tc, Vin); 

%f=100 Hz

vin = @(t) 5 * sin(200*pi*t);
h = 0.00001;
tfinal = 0.05; 
N = round((tfinal-ti)/h); 
Ya = zeros(1,N); 
Xb = zeros(1,N); 
Tc = zeros(1,N); 

Xb(1) = qi; 
Ya(1) = i0; 
Tc(1) = ti; 

for i=1:N-1 
    [Ya(i+1), Xb(i+1)] = RK4second(R,C,L,h,Xb(i), Ya(i), Tc(i),vin); 
    Tc(i+1)=Tc(i)+h; 
end

subplot(3,3,7); 
plot(Tc,R*Ya); 
title('Vout with sine wave input with f = 100 Hz and amplitude 5 V');
ylabel ('Vout/V'); 
xlabel ('t/s'); 

Vin = vin(Tc); 
hold on; 
plot(Tc, Vin); 

%f=500 Hz

vin = @(t) 5*sin(1000*pi*t);
h = 0.00001;
tfinal = 0.01; 
N = round((tfinal-ti)/h); 
Ya = zeros(1,N); 
Xb = zeros(1,N); 
Tc = zeros(1,N); 

Xb(1) = qi; 
Ya(1) = i0; 
Tc(1) = ti; 

for i=1:N-1 
    [Ya(i+1), Xb(i+1)] = RK4second(R,C,L,h, Xb(i), Ya(i), Tc(i),vin); 
    Tc(i+1)=Tc(i)+h; 
end

subplot(3,3,8); 
plot(Tc,R*Ya); 
title('Vout with sine wave input with f = 500 Hz and amplitude 5 V');
ylabel ('Vout/V'); 
xlabel ('t/s'); 

Vin = vin(Tc); 
hold on; 
plot(Tc, Vin); 


