% Exercise-2 : Error Analysis

% This script uses the methods implemented in the previous exercise,
% (Heun, Midpoint and Ralston) and calculates the error function by finding 
% the exact solution of the ODE. The error function is finally used to
% carry the required error analysis.

% The following circuit components/conditions were initialised:
clc; 
clear;
R = 0.5;                % Value of resistor in ohms
L = 1.5 * 10^-3;        % Value of inductor in Henries
t_initial = 0;          % Initial time
i_initial = 0;          % Current at initial time
t_final = 0.0005;       % Final time
h = 0.000005;           % Time Step, maximum value calculated as per 
                        % Nyquist Sampling Rate, since T = 150*10^-6 s,
                        % maximum value of h is 75*10^-6 s.
N = round ((t_final - t_initial/h));

%The given input signal parameters are:

Vin_Amp = 6;                    % Amplitude of signal
T = 150 * 10^-6;                % Time period of cosine wave
f = 1/T;                        % Frequency of wave
w = 2*pi*f                      % Angular velocity
Vin=@(t) Vin_Amp*cos(w*t); % Input Signal
func=@(t,i) (1/L)*(Vin(t)-R*i); % To calculate (di/dt) at any instant

% Next step is to calculate output voltage, Vout, through the three
% proposed methods (Heun, Midpoint and Ralston):

[Heun_t, Heun_Vout]         = Heun(func, i_initial, t_final, h, R, L);
[Midpoint_t, Midpoint_Vout] = Midpoint (func, i_initial, t_final, h, R, L);
[Ralston_t, Ralston_Vout]   = Ralston (func, i_initial, t_final, h, R, L);

% The next step is to calculate the exact solution of the ODE. This Ode is
% solved using integration factor method. Refer to the report for the
% solution of the ODE.

exact = zeros(1,N);

% Time t will be the same for all three methods because h is equal.
exact = ((Vin_Amp*L*w)/((R^2)+(L*L*w*w)))*((L*w*cos(w*Heun_t)) - (R*sin(w*Heun_t)));


% Calculating the errors:
Heun_error      = exact - Heun_Vout;
Midpoint_error  = exact - Midpoint_Vout;
Ralston_error   = exact - Ralston_Vout;

% Plotting the three solutions (Exact solution, Numerical solution and Error):

% Heun

figure(1) ;

subplot (3 ,1 ,1) ;
plot(Heun_t, Heun_Vout, 'r');
title ([ 'Numerical solution by Heuns method, h = ' num2str(h)]);
xlabel('Time [s]');
ylabel('Vout [V]');

subplot (3 ,1 ,2) ;
plot(Heun_t, exact, 'b');
title ([ 'Exact solution by Heuns method, h = ' num2str(h)]);
xlabel('Time [s]');
ylabel('Vout [V]');

subplot (3 ,1 ,3) ;
plot(Heun_t, Heun_error, 'g');
title ([ 'Error in Heuns method, h = ' num2str(h)]);
xlabel('Time [s]');
ylabel('Error [V]');

% Midpoint

figure(2) ;

subplot (3 ,1 ,1) ;
plot(Midpoint_t, Midpoint_Vout, 'r');
title ([ 'Numerical solution by Midpoint, h = ' num2str(h)]);
xlabel('Time [s]');
ylabel('Vout [V]');

subplot (3 ,1 ,2) ;
plot(Midpoint_t, exact, 'b');
title ([ 'Exact solution by Midpoint, h = ' num2str(h)]);
xlabel('Time [s]');
ylabel('Vout [V]');

subplot (3 ,1 ,3) ;
plot(Midpoint_t, Midpoint_error, 'g');
title ([ 'Error in Midpoint, h = ' num2str(h)]);
xlabel('Time [s]');
ylabel('Error [V]');

% Ralston

figure(3) ;

subplot (3 ,1 ,1) ;
plot(Ralston_t, Ralston_Vout, 'r');
title ([ 'Numerical solution by Ralston, h = ' num2str(h)]);
xlabel('Time [s]');
ylabel('Vout [V]');

subplot (3 ,1 ,2) ;
plot(Ralston_t, exact, 'b');
title ([ 'Exact solution by Ralston, h = ' num2str(h)]);
xlabel('Time [s]');
ylabel('Vout [V]');

subplot (3 ,1 ,3) ;
plot(Ralston_t, Ralston_error, 'g');
title ([ 'Error in Ralston, h = ' num2str(h)]);
xlabel('Time [s]');
ylabel('Error [V]');


%%%%%%%

% Finally, we carry out the error analysis in this section. We keep on
% modifying the value of 'h' and the value of maximum error obtained for
% every method is recorded, as are the values of 'h'. An additional plot is
% used to show the order dependency of O(h^2).

K = 20;            % To plot for 20 different values of h

% Initialising new matrices for analysis:

Heun_err            = zeros(1,K);               % Arrays to store max error
Midpoint_err        = zeros(1,K);
Ralston_err         = zeros(1,K);
h_matrix            = logspace (-5, -4.1, K);   % Region where order 2 
                                                % properties are observed.
line_dependency     = zeros(1,K);               % Reference curve for order
                                                % 2 properties.

% Now we initialise the loop to calculate values for different values of
% 'h', keeping in mind the Nyquist Sampling Rate and effect on
% computing speed.

for j = 1:K
    N = round ((t_final - t_initial)/h_matrix(j));
    
    % Solving for the new value of 'h':
    
    [Heun_t, Heun_Vout]         = Heun(func, i_initial, t_final, h_matrix(j), R, L);
    [Midpoint_t, Midpoint_Vout] = Midpoint(func, i_initial, t_final, h_matrix(j), R, L);
    [Ralston_t, Ralston_Vout]   = Ralston(func, i_initial, t_final, h_matrix(j), R, L);
    
    % Solving for the new value of 'exact':
    
        
    exact = ((Vin_Amp*L*w)/((R^2)+(L^2*w^2)))*((L*w*cos(w*Heun_t)) - (R*sin(w*Heun_t)));
    
    % Solving for the new values of 'error':
    
    Heun_error      = abs(exact - Heun_Vout);
    Midpoint_error  = abs(exact - Midpoint_Vout);
    Ralston_error   = abs(exact - Ralston_Vout);
    
    % Solving for maximum errors for every method for every corresponding
    % value of 'h':
    
    Heun_err(j)         = max(Heun_error);
    Midpoint_err(j)     = max(Midpoint_error);
    Ralston_err(j)      = max(Ralston_error);
    line_dependency(j)  = h_matrix(j)^2;
    
end

figure(4);

% Heun

subplot(3, 1, 1);
loglog (h_matrix , Heun_err , '*r') ;
title('Max Error vs Step Size (Heuns Method)');
xlabel('Step Size (s)');
ylabel ('Error in Vout (V)') ;
hold on;
logx = log(h_matrix);
logy = log(Midpoint_err);
Const = polyfit(logx, logy, 1)
hold on
plot(h_matrix, exp(polyval(Const, logx))); 
legend('Numerical error','Exact error');

% Midpoint

subplot(3, 1, 2);
loglog (h_matrix , Midpoint_err , '*r') ;
title('Max Error vs Step Size (Midpoint Method)');
xlabel('Step Size (s)');
ylabel ('Error in Vout (V)') ;
hold on;
logx = log(h_matrix);
logy = log(Midpoint_err);
Const = polyfit(logx, logy, 1)
hold on
plot(h_matrix, exp(polyval(Const, logx)));
legend('Numerical error','Exact error');

% Ralston

subplot(3, 1, 3);
loglog (h_matrix , Ralston_err , '*r') ;
title('Max Error vs Step Size (Ralston Method)');
xlabel('Step Size (s)');
ylabel ('Error in Vout (V)') ;
hold on;
logx = log(h_matrix);
logy = log(Ralston_err);
Const = polyfit(logx, logy, 1)
hold on
plot(h_matrix, exp(polyval(Const, logx)));
legend('Numerical error','Exact error');

figure(5)

% Proof of second order

loglog (h_matrix , Heun_err);
hold on
loglog (h_matrix , Midpoint_err);
loglog (h_matrix , Ralston_err);
loglog(h_matrix, line_dependency);
legend( 'Heuns method', 'Midpoint method', 'Ralstons method', 'y=x^2');
title('Max Error vs Step Size');
xlabel('Step Size (s)');
ylabel ('Error in Vout (V)') ;





