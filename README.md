# NUM-Method-Project
clear, clc, close all
%%%%%%%%%% Driver Script
% Initial Parameters
load EnvironmentalForcing.mat
h = 0.01;
beta = 1;
mu_L = 6; 
mu_I = 10; 
e = 0.001;
A_p = 5000;
P_i = 930.27249;
S_i = P_i/A_p;
L_i = 0.01*S_i;
I_i = 0;
R_i = mu_I*I_i;
B_i = 1;
y0 = [S_i,L_i,I_i,R_i,B_i,0];
T_E = -0.35968 + 0.10789*T - 0.00214*(T.^2);
% Time span
tspan = 0:60;
[t, y] = rk4(@odefun, tspan, y0);



%%%%%%%%%% Functions

% ODEFUN function
function dydt = odefun(t, y)
S = y(1);
I = y(2);
L = y(3);
R = y(4);
P = y(5);
P_b = y(6);

% System of ODEs
dSdt = -beta*S*I + dPdt*(1/A_p);
dLdt = beta*S*I-(1/mu_L)*L + e;
dIdt = (1/mu_L) - (1/mu_I)*I;
dRdt = (1/mu_I)*I;
dP_bdt = (0.1724*P_b - 0.0000212*(P_b^2))*T_E;
dPdt = dP_bdt + (1.33*t_d)*T_E;
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%NEW AND IMPROVED%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Driver Script
load EnvironmentalForcing.mat
beta_max = 1;
mu_L_min = 6; 
mu_I = 10; 
e = 0.001;
A_p = 5000;
P_i = 930.27249;
S_i = P_i/A_p;
L_i = 0.01*S_i;
I_i = 0;
R_i = mu_I*I_i;
B_i = 1;
y = [S_i,L_i,I_i,R_i,B_i,0];
y0 = y;
t0 = 0;
t_d = 1:length(T);

T_B = zeros(size(T)); 
for i = 1:length(T)
    if T(i) <= 0
        T_B(i) = 0;
    elseif T(i) > 0 && T(i) < 35
        T_B(i) =  (0.00241.*(T(i).^2.06737)).*(35 - T(i)).^(0.72859);
    elseif T(i) >= 35
        T_B(i) = 0;
    end
end

beta = beta_max*T_B;

mu_L = 0;
while mu_L < mu_L_min
    for i = 1:length(t_d)
        tinit = 0;
mu_L(i) = sum(tinit:T_B(i));
    end
end


slopes = @(t, y) SLIRP(t, y, beta, mu_L, mu_I, e, T, t_d, A_p);
[t, y] = rk4(slopes, tspan, y0);


% Plotting

plot(t, y(1,:),'k-')
hold on
plot(t, y(2, :), 'b-')
hold on
plot(t, y(3, :), 'g-')
hold on
plot(t, y(4, :), 'r-')
hold on
plot(t, y(5, :), 'p-')
hold on
plot(t, y(6, :), 'm-')
hold off

%% Functions


function slopes = SLIRP(t, y, beta, mu_L, mu_I, e, T, t_d, A_p)

S = y(1);
I = y(3);
L = y(2);
R = y(4);
P = y(5);
P_b = y(6);

T_E = -0.35968 + 0.10789.*T - 0.00214.*(T.^2);
dP_ldt = (1.33.*t_d).*T_E;
dLdt = beta.*S.*I-(1./mu_L).*L + e;
dIdt = (1./mu_L) - (1./mu_I).*I;
dRdt = (1./mu_I).*I;
dP_bdt = (0.1724.*P_b - 0.0000212.*(P_b.^2)).*T_E;
dPdt = dP_bdt + dP_ldt;
dSdt = -beta.*S.*I + dPdt.*(1/A_p);

slopes = [dSdt, dLdt, dIdt, dRdt, dP_bdt, dPdt];
end


function [t,y] = rk4(odefun,tspan,y0)
%UNTITLED Summary of this function goes here
% Detailed explanation goes here
num_steps = length(tspan);
t = zeros(1, num_steps);
y = zeros(length(y0), num_steps);
y(:, 1) = y0;

% Loop through time steps
for i = 1:num_steps-1
    h = tspan(i+1) - tspan(i);
    k1 = odefun(tspan(i), y(:, i));
    k2 = odefun(tspan(i) + h/2, y(:, i) + h*k1/2);
    k3 = odefun(tspan(i) + h/2, y(:, i) + h*k2/2);
    k4 = odefun(tspan(i) + h, y(:, i) + k3*h);

    y(:, i+1) = y(:, i) + (k1(i) + 2*k2(i) + 2*k3(i) + k4(i))/6;
    t(i+1) = tspan(i+1);
end
end




% Sydney Nelson
% u1399500
% Lab 9
clear, clc, close all


load EnvironmentalForcing.mat

B_max = 1;
mu_L_min = 6; 
mu_I = 10; 
e    = 0.001; 
A_P    = 5000; 
P_I  = 1.33*30*(-0.35968+0.10789*15-0.00214*15*15)*30; 
S_I  = P_I/A_P; 
L_I  = S_I*0.01;
I_I  = 0; 
R_I  = I_I*mu_I; 
B_I  = 1; 

tinit=1;
T_B = zeros(size(T)); 
mu_L = zeros(size(T));

for i=1:length(T)

    T_B(i) = solve(T(i));
    mu_L(i) = sum(T_B(tinit:i));

    while(mu_L(i)>mu_L_min)

        tinit=tinit+1;
        mu_L(i) = sum(T_B(tinit:i));

    end
end
mu_L = 1./mu_L; 

p{1} = B_max; 
p{2} = mu_L; 
p{3} = mu_I; 
p{4} = e;    
p{5} = T;    
p{6} = tspan;
p{7} = A_P;    

odefun = @(t,y) slirp_function(t,y,p);

y0(1) = B_I; 
y0(2) = P_I; 
y0(3) = S_I; 
y0(4) = L_I;
y0(5) = I_I; 
y0(6) = R_I; 

[tspan,y] = rk4(odefun,tspan,y0');

B = y(:,1);
P = y(:,2);
S = y(:,3);
L = y(:,4);
I = y(:,5);
R = y(:,6);

hold on
figure;
plot(tspan,P/A_P,'-g',tspan,B/A_P,'--m',tspan,S,'-.k',tspan,L,'--c',tspan,I,':b',tspan,R,'-.r','LineWidth',2);
title('slirp model');
xlabel('time in days');
ylabel('population fraction')
legend('total pop','susceptible pop','berry pop','latent pop ','infected pop','removed pop','Location','northwest')
xlim([0,61]);
hold off

%% Changing B and mul (parts e&f)

figure
Bvals = [1 0.5 1.5];
for k = 1:3
    B_max = Bvals(k);

% beta_max = 1;
mu_L_min = 6; 
mu_I = 10; 
e    = 0.001; 
A_P    = 5000; 
P_I  = 1.33*30*(-0.35968+0.10789*15-0.00214*15*15)*30; 
S_I  = P_I/A_P; 
L_I  = S_I*0.01;
I_I  = 0; 
R_I  = I_I*mu_I; 
B_I  = 1; 

tinit=1;
T_B   = zeros(size(T)); 
mu_L = zeros(size(T));
for i=1:length(T)
    T_B(i) = solve(T(i));
    mu_L(i) = sum(T_B(tinit:i));
    while(mu_L(i)>mu_L_min)
        tinit=tinit+1;
        mu_L(i) = sum(T_B(tinit:i));
    end
end
mu_L = 1./mu_L; 

p{1} = B_max; 
p{2} = mu_L; 
p{3} = mu_I; 
p{4} = e;    
p{5} = T;    
p{6} = tspan;
p{7} = A_P;

odefun = @(t,y) slirp_function(t,y,p);

y0(1) = B_I; 
y0(2) = P_I; 
y0(3) = S_I; 
y0(4) = L_I;
y0(5) = I_I; 
y0(6) = R_I; 

[tspan,y] = rk4(odefun,tspan,y0');

B = y(:,1);
P = y(:,2);
S = y(:,3);
L = y(:,4);
I = y(:,5);
R = y(:,6);

subplot(2, 3, k);
plot(tspan,P/A_P,'-g',tspan,B/A_P,'--m',tspan,S,'-.k',tspan,L,'--c',tspan,I,':b',tspan,R,'-.r','LineWidth',2);
title('b = ', B_max);
xlabel('time in days');
ylabel('population fraction')
legend('total pop','susceptible pop','berry pop','latent pop ','infected pop','removed pop','Location','northwest')
xlim([0,61]);
hold on

end


mul_vals = [6 3 9];
for k = 1:3
    mu_L_min = mul_vals(k);

B_max = 1;
% mul_min = 6; 
mu_I = 10; 
e = 0.001; 
A_P= 5000; 
P_I = 1.33*30*(-0.35968+0.10789*15-0.00214*15*15)*30; 
S_I = P_I/A_P; 
L_I = S_I*0.01;
I_I = 0; 
R_I = I_I*mu_I; 
B_I = 1; 

tinit=1;
T_B= zeros(size(T)); 
mu_L= zeros(size(T));
for i=1:length(T)
    T_B(i) = sall_effect(T(i));
    mu_L(i) = sum(T_B(tinit:i));
    while(mu_L(i)>mu_L_min)
        tinit=tinit+1;
        mu_L(i) = sum(T_B(tinit:i));
    end
end
mu_L = 1./mu_L; 

p{1} = B_max; 
p{2} = mu_L; 
p{3} = mu_I; 
p{4} = e;    
p{5} = T;    
p{6} = tspan;
p{7} = A_P;

odefun = @(t,y) slirp_function(t,y,p);

y0(1) = B_I; 
y0(2) = P_I; 
y0(3) = S_I; 
y0(4) = L_I;
y0(5) = I_I; 
y0(6) = R_I; 

[tspan,y] = rk4(odefun,tspan,y0');

B = y(:,1);
P = y(:,2);
S = y(:,3);
L = y(:,4);
I = y(:,5);
R = y(:,6);

subplot(2, 3, k + 3);
plot(tspan,P/A_P,'-g',tspan,B/A_P,'--m',tspan,S,'-.k',tspan,L,'--c',tspan,I,':b',tspan,R,'-.r','LineWidth',2);
title('mul = ', mu_L_min);
xlabel('time in days');
ylabel('population fraction')
legend('total pop','susceptible pop','berry pop','latent pop ','infected pop','removed pop','Location','northwest')
xlim([0,61]);
hold on

end

function T_B = solve(T)
    if (T<=0)
        T_B = 0;
    elseif (T<35)
        T_B = 0.000241*T^2.06737*(35 - T)^0.72859;
    else
        T_B = 0;
    end
end

function [dydt] = slirp_function(t,y,p)
    B_max = p{1};
    mu_L     = p{2};
    mu_I     = 1/p{3}; 
    e        = p{4};
    T        = p{5};
    day      = p{6};
    A        = p{7};

    B = y(1);
    P = y(2);
    S = y(3);
    L = y(4);
    I = y(5);
    R = y(6);

    if(ceil(t)==floor(t))
        T_new = T(t);
        day_new = day(t);
        mul_new = mu_L(t);
    else 
        t = floor(t);
        T_new = 0.5*(T(t)+T(t+1));
        day_new = 0.5*(day(t)+day(t+1));
        mul_new = 0.5*(mu_L(t)+mu_L(t+1));
    end
    beta = B_max*solve(T_new); 

    TE = -0.35968 +0.10789*T_new-0.00214*T_new^2;
    
    dP_Bdt = (0.1724*B-0.0000212*B*B)*TE;
    dP_dt = (1.33*(day_new+30))*TE+dP_Bdt; 
    dSdt = -beta*S*I+dP_dt/A; 
    dLdt = beta*S*I-mul_new*L+e;  
    dIdt = mul_new*L-mu_I*I;      
    dRdt = mu_I*I;          

    dydt = [dSdt; dIdt; dLdt; dRdt; dP_dt; dP_Bdt];
    dydt = dydt';

end

function [t, y] = rk4(odefun, tspan, y0)
    num_steps = length(tspan);
    num_states = length(y0);
    t = tspan;
    y = zeros(num_steps, num_states);
    y(1,:) = y0(:);
for i = 1:(num_steps-1)
         h = tspan(i+1) - tspan(i);
         k1 = odefun((i), y(i,:));
         k2 = odefun((i)+0.5, y(i,:) + 0.5*k1*h);
         k3 = odefun((i)+0.5, y(i,:) + 0.5*k2*h);
         k4 = odefun((i+1), y(i,:) + k3*h);
         slope = (k1+ 2*k2 + 2*k3 + k4) / 6;
         y(i+1,:) = y(i,:) + h*slope;
         t(i+1) = tspan(i+1);
end 
end 
