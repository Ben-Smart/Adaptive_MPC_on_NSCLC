function [u_opt, cost_val, X_err]= MPC_algorithm(x0, ref, J, sim, ctrl, u_opt, Ts,dt)

%% Sim parameters

N = sim.N; % length of prediction horizon
n_inputs = sim.n_inputs; % the number of inputs

%% cost

j = J*Ts/dt; %stepping through the reference matrix


% if the prediction horizon is longer than the reference signal the 
% controller will aim to keep the state at the last state in the reference 
% signal 
if length(ref(1,:))<(j+N+1)
    ref = [ref ref(:,end).*ones(length(ref(:,1)),N+1)]; 
end

%% Cost function
[a,cost_val, X_err] = costfunction(N,Ts,x0, ref(:,(j+1):(j+N+1)),u_opt,ctrl);
%[a, cost_val, X_err, x_prediction_error] = costfunction(N,Ts,x0, ref(:,(j+1):Ts/dt:(j+N*Ts/dt+1)),u_opt,ctrl);

u_opt=a(1:n_inputs,:); % Setting the first step of the predicted inputs to the actual inputs

%% Drug holidays when the weights are changed online

% if j < 600/0.01
%     
%  u_opt([2 3]) = [0 0];
% ctrl.weights.gamma = [1, 1e100, 1e100];
% 
% end
% 
% 
% 
% if j > 600/0.01
%     
% u_opt([1 3]) = [0 0];
% ctrl.weights.gamma = [1e100, 1e5, 1e100];
% 
% end


end

function [a,cost, X_tilde] = costfunction(N,Ts,x_e,Xref,I,ctrl)
%% COST FUNCTION

X = x_e-Xref(:,1); % error state so equilibrium is at the origin
x0 = X0; % Initial conditions to normalise the states by total 
         % consentration of each molecule

x_norm = [-log(0);x0(2);x0(2);x0(4);x0(4);x0(17);x0(7);x0(7);-log(0);...
          x0(10);x0(10);x0(13);x0(13);x0(14);x0(15);x0(15);x0(17);x0(18)...
          ;x0(19);x0(20);x0(20)]; % each state normalised by the total
                                  % consentration of the molecule for the
                                  % error calculations

X_tilde = X./x_norm; % normalised states for the error calculations

[A,B,c,~,~,~,~,~] = State_Space(x_e,I,Ts); % updated linear model of the 
                                           % system at the current state

%% weights for the cost function                                           
alpha = ctrl.weights.alpha;
beta = ctrl.weights.beta;
gamma = ctrl.weights.gamma;
theta = ctrl.weights.theta;
eta = ctrl.weights.eta;


%% matrices for the traditional cost function

    q = alpha*eye(length(A)) + beta*(c'*c); % Proportional weight on the states
    r = repmat(gamma,1,N); % Weights accosiated with each input


M=zeros(length(A)*(N+1),length(A));                 %Past states affect on future states
C=zeros(length(B(:,1))*(N+1),length(B(1,:))*(N));   %Past inputs affect on future states
Q=zeros(length(A)*(N+1),length(A)*(N+1));           %weights of future states

M(1:length(A),1:length(A)) = eye(length(A));
Q(1:length(A),1:length(A)) = q;

for i = 1:N
    M(length(A)*(i)+1:length(A)*(i+1),1:length(A)) =...
        A*M(length(A)*(i-1)+1:length(A)*(i),1:length(A));
    
    C(length(B(:,1))*(i)+1:length(B(:,1))*(i+1),:) = ...
        A*C(length(B(:,1))*(i-1)+1:length(B(:,1))*(i),:);
    
    C(length(B(:,1))*(i)+1:length(B(:,1))*(i+1),...
        length(B(1,:))*(i-1)+1:length(B(1,:))*(i)) = B;
    
    Q(length(A)*(i)+1:length(A)*(i+1),length(A)*(i)+1:length(A)*(i+1)) = q;
end

R = diag(r);
H = C'*Q*C + R;
F = C'*Q*M;

%% Cost on the rate of change of the inputs 

%%Quadratic term

% main diagonal elements
D = 2*eye(N*length(B(1,:))); 
D(end-length(B(1,:))+1:end,end-length(B(1,:))+1:end) = eye(length(B(1,:)));

%off diagonal elements
off_diag = -ones(1,(N-1)*length(B(1,:)));
D = D + diag(off_diag,length(B(1,:)));
D = D + diag(off_diag,-length(B(1,:)));

D = theta*D; 

%%proportional term
d = zeros(1,N*length(B(1,:)));
d(1:length(B(1,:))) = -2*I(1:length(B(1,:)))'*eye(length(B(1,:)));
d=theta*d;


%% Cost on the integral of the output states


% Quadratic term
P_tilde = zeros(length(A)*(N+1));
for i = 1:N
P_tilde(1 + (i-1)*length(A):(i)*length(A),1 + (i-1)*length(A):(i)*length(A)) = 2*(c')*c;
P_tilde(1 + (i)*length(A):(i+1)*length(A),1 + (i-1)*length(A):(i)*length(A)) = c'*c;
P_tilde(1 + (i-1)*length(A):(i)*length(A),1 + (i)*length(A):(i+1)*length(A)) = c'*c;
end
P_tilde(end-length(A)+1:end,end-length(A)+1:end) = c'*c;

% Proportional term
p_tilde = [zeros(1,length(A(1,:))) (2*(c')*c*X)' zeros(1,length(A(1,:))*(N-1))];



%substitute in relation between error and input
P = C'*(eta*P_tilde)*C;
p = eta*(p_tilde*C + 2*X'*M'*P_tilde*C);


   %% Constraints

 A0 = [];       %There are no state constraints added here, but A0 and b0 
 b0 = [];       %can be edited to add inequaliy constraints on the states

%inequality constraints of inputs, Upper and lower bounds
lb = repmat(ctrl.Lb,1,N); %lower bounds
ub = repmat(ctrl.Ub,1,N); %upper bounds


%last input to guide solver
x_minus = I.*ones(1,N);     % to reduce the runtime the previous input is
                            % told to the solver

%% Solving the cost function

f = (X'*F' + d + p)'; % terms proportional to the inputs
sqr = H + D + P;      % terms with a quadratic relationship to the input
sqr=(sqr+sqr')/2; %added to ensure that no approximation errors result in a none symetic quadratic term

options = ctrl.options;
 [a,cost] = quadprog(sqr,f,A0,b0,[],[],lb,ub,x_minus,options); %finding the inputs that cause the minimum

end