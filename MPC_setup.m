
clear; clc;
%% MPC setup

%% Function Overview: 

%   MPC_setup -     All the parameter settings for the MPC simulation, 
%                   the calculation of the Error/Dose Index and a simple 
%                   plot.  

%   MPC_loop -      initiates the real time loop of MPC, includes the plant 
%                   block of the system through the ode_solver and has been 
%                   left to include the estimator variables, but this has 
%                   not been added yet, but a Kalman filter can be added. 

%   MPC_algorithm - sets up the cost function for all of the predicted
%                   future steps.

%   X0 -            Initial conditions for the states of a NSCLC cell.

%   Xnp0 -          Initial conditions for the states of a Wild Type cell.

%   XT -            The total consentration of each molecule (not 
%                   inculding pEGFR or pIGFR that activate the system).

%   Xol  -          The Open Loop response of either the Wild type cells 
%                   or the NSCLC cells (depending on initial conditions).

%   State_Space -   the linearisation of the system, add how the linear
%                   model is produced from the state conditions. The NSCLC 
%                   is kept as an example, please replace this.

%   ode_solver -    acts as the Plant block in the MPC, it is currently 
%                   the non-linear model of the NSCLC system. Add the 
%                   system to be acted apon under function 'f'.

%   menten -        The 6 Menten parameters for the 3 michaleas menten 
%                   terms describing how the inputs interact with the
%                   model.

% The controller is set to follow the Wild type reference and only has 
% input constraints. Currently there is no Estimator added as the 
% controller is set to perfectly measure all of the internal states of the 
% system. When this is not the case, please use the vector x_esti as the 
% output of the kalman filter estimation of the states from the observable 
% states within MPC_loop.


%% parameters

Par.time.Ts = 1;                      % sample time of the controller(h)
Par.time.dt = 0.01;                   % Time step for the ODE solver(h)
Par.time.Tend =2200;                  % final time (h)

Par.sim.x_hat = X0;                   % the state of the system at time j
Par.sim.x0 = X0;                      % initial state
Par.sim.N = 10;                       % steps in the prediction horizon
Par.sim.n_inputs = 3;                 % the number of inputs

Par.Init_input.INPUT_act = [0 0 0]';  % initial inputs

Par.ctrl.weights.alpha = 0;           % cost on the internal state errors
Par.ctrl.weights.beta = 0;            % cost on the outputs
Par.ctrl.weights.gamma = [1 1e5 1e9]; % cost on the input
Par.ctrl.weights.theta = 0;           % rate of change of the input
Par.ctrl.weights.eta = 1;             % cost on the intergral of the error 
                                      % signal of the outputs
Par.ctrl.Lb = [0 0 0];                % lower bound for the control input
Par.ctrl.Ub = [1 1 1];                % upper bound for the control input
Par.ctrl.options = optimoptions...    % options for the quadratic solver 
    ('quadprog','display','none');    % which reduce the runtime of the 
                                      % programe, but will show no warning
                                      % messages.
                                      
%% Signals

xol = Xol(Par.time,Par.sim.x0);       % Free response of the NSCLC Cell
xnp = Xol(Par.time,Xnp0)';            % Wild type response
xT = XT;                              % Total concentration of molecules

ref.Xref = xnp;                       % set Wild Type as the reference. 

%% MPC function

[STATES, Par, endtime] = MPC_loop(Par, ref);


%% Controller Performance 

% Error Index - integral squared error of the output

[~,~,C,~,~,~,~,~] = State_Space(Par.sim.x0,Par.Init_input.INPUT_act,Par.time.Ts);
const = ((Par.time.Ts^2))/4;
err1=STATES.err;
err2=[err1(:,2:end) zeros(length(C(1,:)),1)];
err = err1+err2;
err = C'*C*err(:,1:end-1);
EI = const*sum(sum(err.^2)); 
EI = round(EI,3);

% Dose Index - Integral of the inputs 
INPUT_act = STATES.input;
DI_1 = Par.time.Ts*(sum(INPUT_act(1,:)) - INPUT_act(1,end)/2) ;
DI_1 = round(DI_1,0);
DI_2 = Par.time.Ts*(sum(INPUT_act(2,:)) - INPUT_act(2,end)/2) ;
DI_2 = round(DI_2,0);
DI_3 = Par.time.Ts*(sum(INPUT_act(3,:)) - INPUT_act(3,end)/2) ;
DI_3 = round(DI_3,2);

%% plots 


x_real = STATES.x_real; % Plant states
INPUT_act = STATES.input; % input profiles
t = linspace(0,length(x_real(1,:))*Par.time.dt,length(x_real(1,:))); %time step of output states
tnp = linspace(0,Par.time.Tend,Par.time.Tend/Par.time.Ts); % time steps of the inputs
t=t(1:end-1);
x_real=x_real(:,1:end-1);
 
 % normalisation of reference and outputs by the total molecule
 % consentration
 
 x_real([8 12],:) = x_real([8 12],:)./xT([7 13]);
 xnp(:,[8 12]) = xnp(:,[8 12])./xT([7 13])';
 xol([8 12],:) = xol([8 12],:)./xT([7 13]);

figure
n = 18; % font size
m = 2.5;% line thickness
   
   
subplot(6,3,[1 4 7 10 13 16])
plot(t(1:end),xnp(1:length(x_real(1,1:end)),8),'LineWidth',m)
hold on
plot(t(1:end),xol(8,1:length(x_real(1,1:end))),'LineWidth',m)
plot(t(1:end),x_real(8,1:end),'k','LineWidth',m)
axis([0 125 0 0.7])
xticks([0 25 50 75 100 125])
xticklabels({'0', '25', '50', '75','100', '125'})
%xlabel('Time(min)')
ylabel(' pERK / Total ERK')
title('y1 - ERK')
set(findall(gcf,'-property','FontSize'),'FontSize',n)

subplot(6,3,[2 5 8 11 14 17])
plot(t,xnp(1:length(x_real(1,:)),12),'LineWidth',m)
hold on
plot(t(1:end),xol(12,1:length(x_real(1,1:end))),'LineWidth',m)
plot(t,x_real(12,:),'k','LineWidth',m)
 axis([0 2200 0 0.7])
%xlabel('Time(min)')
xticks([0 440 880 1320 1760 2200])
xticklabels({'0', '440','880', '1320', '1760', '2200'})
ylabel(' pAkt / Total Akt')
legend('Wild Type', 'Free', 'MIMO')
legend('boxoff')
legend('Location','best')
title('y2 - Akt')
txt = 'Time(min)';
text(-2200,-0.07,txt);
text(740,-0.075,txt);
set(findall(gcf,'-property','FontSize'),'FontSize',n)
hold on

subplot(6,3,[3 6])
plot(tnp,INPUT_act(1,:),'k','LineWidth',m)
axis([0 2200 0 1.100])
xticks([])
title('I_1')
txt = ['DI_1 = ' num2str(DI_1)];
text(1250,0.95,txt);
ylabel('\mu M')
set(findall(gcf,'-property','FontSize'),'FontSize',n)
set(get(gca,'title'),'Position',[1100 1.075 1])
hold on



subplot(6,3,[9 12])
plot(tnp,INPUT_act(2,:),'k','LineWidth',m)
axis([0 2200 0 1.100])
ylabel('\mu M')
xticks([])
title('I_2')
txt = ['DI_2 = ' num2str(DI_2)];
text(1250,0.95,txt);
set(findall(gcf,'-property','FontSize'),'FontSize',n)
set(get(gca,'title'),'Position',[1100 1.075 1])
hold on


subplot(6,3,[15 18])
plot(tnp,INPUT_act(3,:),'k','LineWidth',m)
axis([0 2200 0 1.100])
ylabel('\mu M')
title('I_3')
txt = ['DI_3 = ' num2str(DI_3)];
text(1250,0.95,txt);
xticks([0 440 880 1320 1760 2200])
xticklabels({'0', '440','880', '1320', '1760', '2200'})
txt = 'Time(min)';
text(740,-0.38,txt);
hold on
set(findall(gcf,'-property','FontSize'),'FontSize',n)
set(get(gca,'title'),'Position',[1100 1.075 1])
% sgtitle('T_s = 1 min') 

annotation('textbox',...
    [0.13 0.775 0.99 0.15],...
    'String',{['EI = ' num2str(EI)]},...
    'FontSize',n,...
    'FontName','Arial',...
    'LineStyle','none',...
    'LineWidth',2,...
    'Color','k');

