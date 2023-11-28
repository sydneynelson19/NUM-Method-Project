# NUM-Method-Project
clear, clc, close all
%%%%%%%%%% Driver Script
% Time span
tspan = linspace(0, 10, 1000);
[t, y] = rk4(@odefun, tspan, y0);




%%%%%%%%%% Functions

% ODEFUN function
function dydt = odefun(t, y)
beta = 1;
mu_L_inv = 6; 
mu_I_inv = 10; 
e = 0.001;
A_p = 5000;
P_i = 930.27249;
S_i = P_i/A_p;
L_i = 0.01*S_i;
I_i = 0;
R_i = mu_I_inv*I_i;
B_i = 1;

S = y(1);
I = y(2);
L = y(3);
R = y(4);
P = y(5);
P_b = y(6);

% System of ODEs
dSdt = -beta*S*I + (1/A_p)*(1/mu_L_inv)*P;
dLdt = beta*S*I-(1/mu_L_inv)*L + e;
dIdt = (1/mu_L_inv) - (1/mu_I_inv)*I;
dRdt = (1/mu_I_inv)*I;
dPdt = (1/mu_L_inv)*P_b + (1/mu_L_inv)*P;
dP_bdt = (0.1724*P_b - 0.0000212*(P_b^2))*T_E*A_p;

dydt = [dSdt; dIdt; dLdt; dRdt; dPdt; dP_bdt];
    end
% SLIRP FUNCTION
function slopes = SLIRP

end

% RUNGE KUTTA 4TH ORDER FUNCTION
function [t, y] = rk4(odefun, tspan, y0)
 t0 = tspan(1);
tf = tspan(end);
h = tspan(2) - tspan(1);
num_steps = length(tspan);

t = zeros(num_steps, 1);
y = zeros(num_steps, length(y0));

% initial Conditions
t(1) = t0;
y(1, :) = y0;

% Runge Kutta
    for i = 2:num_steps
        t(i) = t(i-1) + h;
        k1 = h*odefun(t(i-1), y(i-1, :));
        k2 = h*odefun(t(i-1) + 0.5*h, y(i-1, :) + 0.5*k1);
        k3 = h*odefun(t(i-1) + 0.5*h, y(i-1, :)+0.5*k2);
        k4 = h*odefun(t(i-1) + h, y(i-1, :) + k3);
        
        y(i, :) = y(i-1, :)+(1/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
end
