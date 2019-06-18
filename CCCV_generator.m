 %% Script for generating CCCV Charging Profiles for a Given Battery Chemistry via a PI Loop implementation
% By: Zach Gima

% % % Current | Positive <=> Discharge, Negative <=> Charge
clearvars
clc
close all

fs = 25;

%% Initialize Parameters
% General Parameters, NCR18650PF (Bosch NCA cell)
run param/params_NCA
dt = 0.5; %seconds
p.delta_t = dt;
p = set_discretization(p);

V0 = 2.8; % SOC = 0%: --> 2.5V; SOC = 100%: --> 4.2V
T_amb = 25; 
V_max = 4; %V

capacity = 2900; %mAh
C_rate = 0.5; %[-]
I_chrgmax = -capacity*C_rate/1000; %A
I_chrgcutoff = -0.01; %A i.e. 5mA cutoff

% PI Charge Controller Settings
Kp = 1.8; %0.75 %2 for just P-controller
Ki = 0.08*dt;
Kd = 0;
T_i = 10; % reset time [seconds]
Kt = dt/T_i; % reset gain: large Kt = short reset time; small = long; typically should be larger than Ki

% error function
% Because charging current is defined as negative in the model, if we use
% the typical PI error law (ref - output), then we end up needing to
% implement a reverse-action PI controller (i.e. decrease in error -->
% increase in control). By changing the sign on the voltage error, we can
% more easily implement the PI control here as direct.
% err = @(V) V_max - V;
err = @(V) V - V_max;

%% Simulate Model & Determine Charging Profile
% Initialize 
k = 1;
V = zeros(2,1);
I = zeros(2,1);
V(1) = V0;
I(1) = 0;
Int = 0; % integral term
Cur = I(1); %I(2);
% err_old = 0;
% err_history = zeros(2,1);
% err_history(1) = 0;

% Determine Current Trajectory until Cutoff Current reached
while k < 3 || abs(I(k)) > abs(I_chrgcutoff) 
    %% Simulate Model for initial states and voltage output
    % Initialize states for each time step
    [p,x0_nom] = init_spmet_nocasadi(p,V0,T_amb,0);  % generate initial conditions for spmet & load inputs, parameters

    % Simulate Voltage from Model
    [~,V(k+1)] = ode_spmet(x0_nom,Cur,p);
    V0 = V(k+1); %update voltage used to initialize states at next time step
    
    %% PI w/ Anti-Windup Routine
    % Adopted from Aström, Feedback Systems
    
    % Compute updated voltage error
    err_V = err(V(k+1));

    %%% Input saturation reset form
    % Proportional & Derivative terms
    Prop = Kp*err_V;
    Deriv = 0;

    %Compute Current
    % Positive error between Vmax and V should correspond to a negative
    % (i.e. charging) current
    % I_ctrl = controller input to actuator; I(k) = actuator output
    I_ctrl = Prop + Int + Deriv;
    I(k+1) = max(I_ctrl,I_chrgmax); % bound the max current; here use max instead of min because charge current is (-)
    
    % Update Integral Term:
    % controller-actuator error: 0 when no saturation
    err_sat = I(k+1) - I_ctrl;  
    Int = Int + Ki*err_V + Kt*err_sat;
    
    % Debug variables
%     err_V_track(k,1) = err_V; %debug
%     err_sat_track(k,1) = err_sat; %debug
%     Int_track(k,1) = Int; %debug
%     Prop_track(k,1) = Prop; %debug
    
    %%% Discrete velocity form
%     delta_I = -(Kp*(1+dt/T_i)*err_new - Ki*err_old);
%     I(k+1) = I(k) + delta_I;
%     I(k+1) = max(I(k+1),I_chrgmax); % bound the max current; here use max instead of min because charge current is (-)
%     err_old = err_new;
    
    %% Debugging stuff
%     if mod(k,1000) == 0
%         fprintf('Current is %1.6f \n',Cur);
%         fprintf('Voltage is %1.6f \n',V(k+1));
%         
%         % Debugging plots
%         %%% errors
%         figure('Position', [100 100 900 700])
%         subplot(3,1,1)
%         plot(Ki*err_V_track)
%         title('Proportional Error')
%         set(gca,'Fontsize',fs)
% 
%         subplot(3,1,2)
%         plot(Kt*err_sat_track)
%         title('Saturation error')
%         set(gca,'Fontsize',fs)
% 
%         subplot(3,1,3)
%         plot(Int_track)
%         title('Total Integral Error')
%         set(gca,'Fontsize',fs)
% 
%         %%% Voltage and current
%         figure('Position', [100 100 900 700])
%         subplot(2,1,1);
%         plot(I)
%         xlabel('Time (s)')
%         ylabel('Current (A)')
%         set(gca,'Fontsize',fs)
% 
%         subplot(2,1,2);
%         plot(V)
%         xlabel('Time (s)')
%         ylabel('Voltage (V)')
%         set(gca,'Fontsize',fs)
%         
%         close all
%     end
    
%     if mod(k,1200) == 0
%         figure
%         hold on
%         plot(Int_track)
%         plot(Prop_track,'b')
% %         title('Integral & Proportional Error')
%         legend('Integral','Proportional');
%         
%         figure
%         plot(I)
%         figure
%         plot(V)
%         
%         close all
%     end
    
    %% Update k and Cur
    Cur = I(k+1);
    k = k + 1;
end
time = (dt:dt:k*dt)';

figure('Position', [100 100 900 700])
subplot(2,1,1);
plot(time,I)
% xlabel('Time (s)')
ylabel('Current (A)')
set(gca,'Fontsize',fs)

subplot(2,1,2);
plot(time,V)
xlabel('Time (s)')
ylabel('Voltage (V)')
set(gca,'Fontsize',fs)

% Debugging plots
%%% errors
% figure('Position', [100 100 900 700])
% subplot(3,1,1)
% plot(Ki*err_V_track)
% title('Proportional Error')
% set(gca,'Fontsize',fs)
% 
% subplot(3,1,2)
% plot(Kt*err_sat_track)
% title('Saturation error')
% set(gca,'Fontsize',fs)
% 
% subplot(3,1,3)
% plot(Int_track)
% title('Total Integral Error')
% set(gca,'Fontsize',fs)
