i0=0;                                                                %initial current value
R=0.5;                                                               %resistor value = 0.5ohm 
L=0.0015;                                                            %Inductor value = 1.5mH
tf =  5 * 10^-5;                                                     %stop time
h  = 0.000000001;                                                    %step size

Vin = @(t) 3.5*exp(-(t.^2)/(150e-12));                               %input voltage

func=@(t,i) (1/L)*(Vin(t)-R*i);                                      %function handle->i'=(1/L)*(Vin(t)-R*i)

[Heun_t,Heun_Vout] = Heun(func,i0,tf,h,R,L);                         %calling Heun function

figure(1)
plot(Heun_t,Heun_Vout);                                              %plotting output voltage vs time for Heun method
hold on
title('Vout against time for Heun method');
xlabel('Time / s');
ylabel('Vout / V');
legend('Vout', 'Vin');


[Midpoint_t,Midpoint_Vout] = Midpoint(func,i0,tf,h,R,L);             %calling Midpoint function

figure(2)
plot(Midpoint_t,Midpoint_Vout);                                      %plotting output voltage vs time for Midpoint method
hold on
title('Vout against time for Midpoint method');
xlabel('Time / s');
ylabel('Vout / V');
legend('Vout', 'Vin');


[Ralston_t,Ralston_Vout] = Ralston(func,i0,tf,h,R,L);                %calling Ralston function

figure(3)
plot(Ralston_t,Ralston_Vout);                                        %plotting output voltage vs time for Ralston method
hold on
title('Vout against time for Ralston method');                              
xlabel('Time / s');
ylabel('Vout / V');
legend('Vout', 'Vin');
