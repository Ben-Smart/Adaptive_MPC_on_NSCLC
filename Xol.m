
function xol = Xol(time,x0)
%% Open loop simulation - no inputs
% the initial conditions (x0) used dictates wheather the Wild type or NSCLC
% cells are being simulated

Ts = time.Ts; %controller time step
dt = time.dt; %ode solver step
Tend = time.Tend; %end time for sinulation


I =  [0 0 0]; % zero vector for the inputs
z = [x0' I];  % states for the ode solver
 tspan = linspace(0,Tend,Tend*Ts/dt); %time span to solve for
[ ~ , z ] = ode_solver(z, tspan); % solving the open loop simulation
xol = z(:,1:length(x0))'; % open loop states

end