function [t,Vout] = Ralston(func,i0,tf,h,R,L)
    %R resistor value 
    %L inductor value 
    %tf stop time 
    %i0 inital current value 
    %func function handle
    ti = 0;                                           %initial time
    N = round((tf-ti)/h);                             %N = number of iterations
    t(1)=ti;                                          %set t(1) as initial time value, ti
    i = zeros(1,N);                                   %initialise i array to zeros
    i(1) = i0;                                        %set i(1) as initial current value, i0
    Vout(1) = L*func(t(1),i(1));                      %initial value of output voltage calculated->Vout(1)=L*di(1)/dt(1);  
    
    for j=1:N-1                                       % loop for N iterations
        t(j+1) = t(j) + h;                            %Updating the next value of time with t(j) + h, where t(j) is the previous time value and h is the step size
        k1 = func(t(j),i(j));                         %gradient at t(j)
        k2 = func(t(j) + 2*h/3, i(j) + 2*h*k1/3);     %second gradient at t(j) + 2*h/3
        i(j+1) = i(j) + h*(k1 + 3*k2)/4;              %next value of i calculated from previous values of t,i 
        Vout(j+1) = L*func(t(j+1),i(j+1));            %Next value of output voltage calculated-> Vout(j+1) = L*di(j+1)/dt(j+1)
    end
end