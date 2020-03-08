%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AEM 349-001 HW 1

clear all
close all

% Givens
dydt = @(t,y) t.*y-3.*y;
exact = @(t) exp(((t.^2)/2)-3.*t);
y0 = 1;
t0 = 0;
tf = 2;
dt = [ 0.5 0.2 0.1];
% For sach step size, Use all methods
for i = 1:length(dt)
    h = dt(i);
    [t_out_euler,y_out_euler] = euler(dydt,y0,t0,tf,h);
    [t_out_heun, y_out_heun]  = heun(dydt,y0,t0,tf,h);
    [t_out_rk2,  y_out_rk2]   = rk2(dydt,y0,t0,tf,h);
    [t_out_rk4,  y_out_rk4]   = rk4(dydt,y0,t0,tf,h);
    fprintf("\n Step Size = %f\n",h)
    % For each step, Determine Errors
    for j = 1:length(t_out_euler) 
        act = exact(t_out_euler(j));
        err_eu = abs((act-y_out_euler(j))/act)*100;
        err_heun = abs((act-y_out_heun(j))/act)*100;
        err_rk2 = abs((act-y_out_rk2(j))/act)*100;
        err_rk4 = abs((act-y_out_rk4(j))/act)*100;
        % Print out Table
        fprintf('\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
            t_out_euler(j),act,y_out_euler(j),err_eu,y_out_heun(j), ...
            err_heun,y_out_rk2(j),err_rk2,y_out_rk4(j),err_rk4);
    end
    
end
hold on
plot(t_out_euler,y_out_euler);
plot(t_out_heun,y_out_heun);
plot(t_out_rk2,y_out_rk2);
plot(t_out_rk4,y_out_rk4);
plot(t_out_euler,act);
ylabel('y(t)')
xlabel('t')
legend('Euler','Heun','RK2','RK4','Actual')
box on

%%%%%%%%%%%%%%%%%%%%%%%
% Euler's explicit Method
function [t_out_euler, y_out_euler] = euler(dydt,y0,t0,tf,h)

    N = (tf-t0)/h;
    tn = t0:h:tf; % Steps
    % Set Initial Values
    y_out_euler(1) = y0;
    t_out_euler = t0;
    for i=1:length(tn)-1
        % Euler Steps
        t_out_euler(i+1) = t_out_euler(i) + h;
        y_out_euler(i+1) = y_out_euler(i) + ...
            h*double(dydt(t_out_euler(i),y_out_euler(i)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%
% Heun's Method
function [t_out_heun, y_out_heun] = heun(dydt,y0,t0,tf,h)
    % Set Intial Values
    t_out_heun(1) = t0;
    y_out_heun(1) = y0;
    tn = t0:h:tf; % Steps
    for i = 1:length(tn)-1;
        % Heun's Method Steps
        t_out_heun(i+1) = t_out_heun(i) + h;
        m1 = double(dydt(t_out_heun(i),y_out_heun(i)));
        m2 = double(dydt((t_out_heun(i)+h),(y_out_heun(i)+h*m1)));
        y_out_heun(i+1) = y_out_heun(i) + double(h*((m1+m2)/2));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runge Kutta 2nd Order
function [t_out_rk2,y_out_rk2] = rk2(dydt,y0,t0,tf,h)
    % Set Initial Values
    t_out_rk2(1) = t0;
    y_out_rk2(1) = y0;
    N = (tf-t0)/h;
    t=t0;
    y=y0;
    for i =1:N
        % rk2 Steps
        k1 = h*dydt(t,y);
        k2 = h*dydt(t+0.5*h,y+0.5*k1);
        dy=k2;
        t=t+h;
        y=y+dy;
        t_out_rk2(i+1) = t;
        y_out_rk2(i+1) = y;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runge Kutta 4th Order
function [t_out_rk4,y_out_rk4] = rk4(dydt,y0,t0,tf,h)
    % Set Initial Values
    t_out_rk4(1) = t0;
    y_out_rk4(1) = y0;
    N = (tf-t0)/h;
    t=t0;
    y=y0;
    for i = 1:N
        % rk4 Steps
        k1 = h*double(dydt(t,y));
        k2 = h*double(dydt((t+h/2),(y+k1/2)));
        k3 = h*double(dydt((t+h/2),(y+k2/2)));
        k4 = h*double(dydt((t+h),(y+k3)));
        dy = (1/6)*(k1+2*k2+2*k3+k4);
        t = t+h;
        y = y+dy;
        t_out_rk4(i+1) = t;
        y_out_rk4(i+1) = y;
    end
end