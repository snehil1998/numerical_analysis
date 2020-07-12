function [Yi, Xi] = RK4second(R,C,L,h,X,Y,t,vin) %definining the function required

%generates Yi and Xi, given the Y and X
vin1 = vin(t);
coup=@(t,X,Y) Y; 
pre_coup=@(t,X,Y) (vin1 - (R*Y) - X/C)/L; %circuit equation derived above expressing dy/dt  


% setting the method of calculating Yi and Xi 

K_1 = feval(coup, t, X, Y); 
L_1 = feval(pre_coup, t, X, Y);

K_2 = feval(coup, t+h/3, X+h*K_1/3, Y+h*L_1/3); 
L_2 = feval(pre_coup, t+h/3, X+h*K_1/3, Y+h*L_1/3); 

K_3 = feval(coup, t+2*h/3, X+h*(-K_1/3 + K_2), Y+h*(-L_1/3 + L_2)); 

L_3 = feval(pre_coup, t+2*h/3, X+h*(-K_1/3 + K_2), Y+h*(-L_1/3 + L_2)); 

K_4 = feval(coup, t+h, X+h*(K_1 - K_2 + K_3), Y+h*(L_1 - L_2 + L_3)); 

L_4 = feval(pre_coup, t+h, X+h*(K_1 - K_2 + K_3), Y+h*(L_1 - L_2 + L_3)); 

X_new = (K_1 + 3*K_2 + 3*K_3 + K_4)/8; 
Y_new = (L_1 + 3*L_2 + 3*L_3 + L_4)/8; 

Yi = Y + h*Y_new; %Yi formula 
Xi = X + h*X_new; %Xi formula 

end