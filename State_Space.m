               %% State_Space 
% the linearisation of the system - how the linear model is produced
% from the state conditions. The NSCLC system is kept as an example.

function [Ad,Bd,Cd,Dd,A,B,C,D] = State_Space(x,Iss,Ts)
%% Parameters in the model

p = zeros(22,3);


% Parameter:   id =  gamma_EGFR, name = gamma_EGFR
	p(3,1)=0.02;

% Reaction: id = SOS_conformational_activation, name = SOS conformational activation	% Local Parameter:   id =  k_SOS_E, name = k_SOS_E
	p(1,1)=694.731;
	% Local Parameter:   id =  n_SOS, name = n_SOS
	p(1,2)=1.0;
	% Local Parameter:   id =  KM_SOS_E, name = KM_SOS_E
	p(1,3)=6086070.0;


% Reaction: id = kRas_Activation, name = kRas Activation	% Local Parameter:   id =  k_Ras_SOS, name = k_Ras_SOS
	p(2,1)=32.344;
	% Local Parameter:   id =  n_Ras_SOS, name = n_Ras_SOS
	p(2,2)=1.0;
	% Local Parameter:   id =  KM_Ras_SOS, name = KM_Ras_SOS
	p(2,3)=35954.3;


% Reaction: id = EGFR_degradation, name = EGFR degradation


% Reaction: id = ERK_activationBy_Mek, name = ERK activation by Mek	% Local Parameter:   id =  k_ERK_MekActive, name = k_ERK_MekActive
	p(4,1)=9.85367;
    p(4,2) = 1;
	% Local Parameter:   id =  KM_ERK_MekActive, name = KM_ERK_MekActive
	p(4,3)=1007340.0;


% Reaction: id = SOS_deactivationBy_P90, name = SOS deactivation by P90	% Local Parameter:   id =  k_D_SOS_P90Rsk, name = k_D_SOS_P90Rsk
	p(5,1)=161197.0;
	% Local Parameter:   id =  n_D_SOS, name = n_D_SOS
	p(5,2)=1.0;
	% Local Parameter:   id =  KM_D_SOS_P90Rsk, name = KM_D_SOS_P90Rsk
	p(5,3)=896896.0;

% Reaction: id = SOS_activationBy_IGF, name = SOS activation by IGF	% Local Parameter:   id =  k_A_SOS_I, name = k_A_SOS_I
	p(6,1)=500.0;
	% Local Parameter:   id =  n_A_SOS_I, name = n_A_SOS_I
    p(6,2)=1.0;
	% Local Parameter:   id =  KM_A_SOS_I, name = KM_A_SOS_I
	p(6,3)=100000.0;

% Reaction: id = PI3KCA_activationBy_IGF1R, name = PI3KCA activation by IGF1R	% Local Parameter:   id =  k_PI3K_IGF1R, name = k_PI3K_IGF1R
	p(7,1)=10.6737;
	% Local Parameter:   id =  n_PI3K_I, name = n_PI3K_I
	p(7,2)=1.0;
	% Local Parameter:   id =  KM_PI3K_IGF1R, name = KM_PI3K_IGF1R
	p(7,3)=184912.0;

% Reaction: id = PI3KCA_activationBy_EGF, name = PI3KCA activation by EGF	% Local Parameter:   id =  k_PI3K_EGF1R, name = k_PI3K_EGF1R
	p(8,1)=10.6737;
	% Local Parameter:   id =  n_PI3K_E, name = n_PI3K_E
    p(8,2)=1.0;
	% Local Parameter:   id =  KM_PI3K_EGF1R, name = KM_PI3K_EGF1R
	p(8,3)=184912.0;

% Reaction: id = Akt_activationBy_PI3KCA, name = Akt activation by PI3KCA	% Local Parameter:   id =  k_AKT_PI3K, name = k_AKT_PI3K
	p(9,1)=0.0566279;
	% Local Parameter:   id =  n_AKT_PI3K, name = n_AKT_PI3K
	p(9,2)=1.0;
	% Local Parameter:   id =  KM_AKT_PI3K, name = KM_AKT_PI3K
	p(9,3)=653951.0;

% Reaction: id = Akt_deactivation, name = Akt deactivation	% Local Parameter:   id =  kd_AKT, name = kd_AKT
	p(10,1)=0.005;

% Reaction: id = ERK_deactivationBy_PP2A, name = ERK deactivation by PP2A	% Local Parameter:   id =  k_ERKactive_PP2A, name = k_ERKactive_PP2A
	p(11,1)=8.8912;
	% Local Parameter:   id =  n_ERKactive_PP2A, name = n_ERKactive_PP2A
	p(11,2)=1.0;
	% Local Parameter:   id =  KM_ERKactive_PP2A, name = KM_ERKactive_PP2A
	p(11,3)=3496490.0;

% Reaction: id = PI3KCA_activationBy_kRas, name = PI3KCA activation by kRas	% Local Parameter:   id =  k_PI3K_Ras, name = k_PI3K_Ras
	p(12,1)=0.0771067;
	% Local Parameter:   id =  n_PI3K_Ras, name = n_PI3K_Ras
	p(12,2)=1.0;
	% Local Parameter:   id =  KM_PI3K_Ras, name = KM_PI3K_Ras
	p(12,3)=272056.0;

% Reaction: id = Raf_activationBy_kRas, name = Raf activation by kRas	% Local Parameter:   id =  k_Raf_RasActive, name = k_Raf_RasActive
	p(13,1)=0.884096;
	% Local Parameter:   id =  n_Raf_RasActive, name = n_Raf_RasActive
	p(13,2)=1.0;
	% Local Parameter:   id =  KM_Raf_RasActive, name = KM_Raf_RasActive
	p(13,3)=62464.6;

% Reaction: id = Mek_activationBy_Raf, name = Mek activation by Raf	% Local Parameter:   id =  k_Mek_PP2A, name = k_Mek_PP2A
	p(14,1)=185.759;
	% Local Parameter:   id =  n_Mek_PP2A, name = n_Mek_PP2A
	p(14,2)=1.0;
	% Local Parameter:   id =  KM_MekPP2A, name = KM_MekPP2A
	p(14,3)=4768350.0;

% Reaction: id = Raf_deactivationBy_Akt, name = Raf deactivation by Akt	% Local Parameter:   id =  k_Raf_AKT, name = k_Raf_AKT
	p(15,1)=15.1212;
	% Local Parameter:   id =  n_Raf_AKT, name = n_Raf_AKT
	p(15,2)=1.0;
	% Local Parameter:   id =  KM_Raf_AKT, name = KM_Raf_AKT
	p(15,3)=119355.0;

% Reaction: id = Ras_deactivation, name = Ras deactivation by RasGab	% Local Parameter:   id =  k_RasActiveRasGap, name = k_RasActiveRasGap
	p(16,1)=1509.36;
	% Local Parameter:   id =  n_RasActiveRasGap, name = n_RasActiveRasGap
	p(16,2)=1.0;
	% Local Parameter:   id =  KM_RasActiveRasGap, name = KM_RasActiveRasGap
	p(16,3)=1432410.0;

% Reaction: id = Mek_deactivation, name = Mek deactivation by PP2A	% Local Parameter:   id =  k_MekActivePP2A, name = k_MekActivePP2A
	p(17,1)=2.83243;
	% Local Parameter:   id =  n_MekActivePP2A, name = n_MekActivePP2A
	p(17,2)=1.0;
	% Local Parameter:   id =  KM_MekActivePP2A, name = KM_MekActivePP2A
	p(17,3)=518753.0;

% Reaction: id = IGFR_active_degradation, name = IGFR active degradation
    % Parameter:   id =  gamma_IGFR, name = gamma_IGFR
	p(18,1)=0.02;

% Reaction: id = PI3KCA_deactivation, name = PI3KCA deactivation
% Parameter:   id =  kd_PI3K_a, name = kd_PI3K_a
	p(19,1)=0.005;


% Reaction: id = Raf_deactivation, name = Raf deactivation by RafPP	% Local Parameter:   id =  k_RasActive_RafPP, name = k_RasActive_RafPP
	p(20,1)=0.126329;
	% Local Parameter:   id =  n_RasActive_RafPP, name = n_RasActive_RafPP
	p(20,2)=1.0;
	% Local Parameter:   id =  KM_RasActive_RafPP, name = KM_RasActive_RafPP
	p(20,3)=1061.71;

	
% Reaction: id = P90_activationBy_ERK, name = P90 activation by ERK
% Parameter:   id =  k_P90Rsk_ERKActive, name = k_P90Rsk_ERKActive
	p(21,1)=0.0213697;
% Parameter:   id =  KM_P90Rsk_ERKActive, name = KM_P90Rsk_ERKActive
	p(21,3)=763523.0;
    p(21,2)=1;

% Reaction: id = P90_deactivation, name = P90 deactivation	% Local Parameter:   id =  kd_P90Rsk, name = kd_P90Rsk
	p(22,1)=0.005;
   
    %% %% %% State space formulation %% %% %%
    
%% reactions needed for the inputs to the modeled    
[kon_I1,km_I1,kon_I2,km_I2, kon_I3,km_I3] = menten;

%% A matrix
A = zeros(21,21);

A(1,1) = -p(3,1)*x(1);


[A(2,1),A1] = hill(p(1,:),x(1),x(2),1);
[A(2,21),A(2,3)] = hill(p(5,:),x(21),x(3),0);
[A(2,9),A2] = hill(p(6,:),x(9),x(2),1);

A(2,2) = A1+A2;



A(3,1) = -A(2,1);
A(3,2) = -A(2,2);
A(3,9) = -A(2,9);
A(3,21) = -A(2,21);
A(3,3) = -A(2,3);


[A(4,5),A(4,4)] = hill(p(13,:),x(5),x(4),1);
[A(4,12),A1] = hill(p(15,:),x(12),x(16),0);
[A(4,19),A2] = hill(p(20,:),x(19),x(16),0);

A(4,16) = A1 +A2;


[A(5,3),A(5,15)] = hill(p(2,:),x(3),x(15),0);
[A(5,18),A(5,5)] = hill(p(16,:),x(18),x(5),1);


[A(6,16),A(6,17)] = hill(p(14,:),x(16),x(17),0);
[A(6,14),A(6,6)] = hill(p(17,:),x(14),x(6),1);


[A(7,6),A(7,7)] = hill(p(4,:),x(6),x(7),1);
[A(7,14),A(7,8)] = hill(p(11,:),x(14),x(8),0);


[A(8,6),A(8,7)] = hill(p(4,:),x(6),x(7),0);
[A(8,14),A(8,8)] = hill(p(11,:),x(14),x(8),1);


A(9,9) = -p(18,1)*x(9);


[A(10,9),A1] = hill(p(7,:),x(9),x(10),1);
[A(10,1),A2] = hill(p(8,:),x(1),x(10),1);
A(10,1) = 2*x(1)*A(10,1);
A1 = x(1)*A1;
[A(10,5),A3] = hill(p(12,:),x(5),x(10),1);
A(10,11) = p(19,1)*x(11);
A(10,10) = A1 + A2 + A3;


[A(11,9),A1] = hill(p(7,:),x(9),x(10),0);
[A(11,1),A2] = hill(p(8,:),x(1),x(10),0);
A(11,1) = 2*x(1)*A(11,1);
A1 = x(1)*A1;
[A(11,5),A3] = hill(p(12,:),x(5),x(10),0);
A(11,11) = -p(19,1)*x(11);
A(11,10) = A1 + A2 + A3;


[A(12,11),A(12,13)] = hill(p(9,:),x(11),x(13),0);
A(12,12) = -p(10,1)*x(10);


[A(13,11),A(13,13)] = hill(p(9,:),x(11),x(13),1);
A(13,12) = p(10,1)*x(10);


[A(15,3),A(15,15)] = hill(p(2,:),x(3),x(15),1);
[A(15,18),A(15,5)] = hill(p(16,:),x(18),x(5),0);


[A(16,5),A(16,4)] = hill(p(13,:),x(5),x(4),0);
[A(16,12),A1] = hill(p(15,:),x(12),x(16),1);
[A(16,19),A2] = hill(p(20,:),x(19),x(16),1);
A(16,16) = A1 + A2;


[A(17,16),A(17,17)] = hill(p(14,:),x(16),x(17),1);
[A(17,14),A(17,6)] = hill(p(17,:),x(14),x(6),0);


[A(20,8),A(20,20)] = hill(p(21,:),x(8),x(20),1);
A(20,21) = p(22,1)*x(21);

[A(21,8),A(21,20)] = hill(p(21,:),x(8),x(20),0);
A(21,21) = -p(22,1)*x(21);

%% input interactions with the A matrix

A(10,10) = A(10,10) - kon_I1*Iss(1)/(km_I1+Iss(1));
A(11,11) = A(11,11) - kon_I1*Iss(1)/(km_I1+Iss(1));

A(12,12) = A(12,12) - kon_I2*Iss(2)/(km_I2+Iss(2));
A(13,13) = A(13,13) - kon_I2*Iss(2)/(km_I2+Iss(2));

A(17,17) = A(17,17) - kon_I3*Iss(3)/(km_I3+Iss(3));
A(6,6)   = A(6,6)   - kon_I3*Iss(3)/(km_I3+Iss(3));


%% B matrix

B = zeros(21,3);

    
   % I1

    %I1 - PI3K
    U21=kon_I1*km_I1*x(10)/(km_I1+Iss(1))^2;
    
    %I2 - pPI3K-active
    U22=kon_I1*km_I1*x(11)/(km_I1+Iss(1))^2;
    
    
    % I2
    
    %I2 - Akt_active
    U31=kon_I2*km_I2*x(12)/(km_I2+Iss(2))^2;
    
    %I2 - Akt
    U32=kon_I2*km_I2*x(13)/(km_I2+Iss(2))^2;
    
    
    
    % I3
    
    %I3 - MEK
    U51=kon_I3*km_I3*x(17)/(km_I3+Iss(3))^2;
    %I3 - MEK
    U52=kon_I3*km_I3*x(17)/(km_I3+Iss(3))^2;

 
%I1 as the input

 B(10,1) = -U21;
 B(11,1) = -U22;


%I2 as the input

 B(12,2) = -U31; 
 B(13,2) = -U32; 

%I3 as the input

B(17,3) = -U51; 
B(6,3)  = -U52;


%% C matrix

C = zeros(2,21);

C(1,8) = 1; 

C(2,12) = 1; 

%% D matrix

D=[];

%% create a state space of the system

sys = ss(A,B,C,D);

%% descreatised with Zero Order Hold


sysd1 = c2d(sys,Ts,'zoh');

Ad = sysd1.A;

Bd = sysd1.B;

Cd = sysd1.C;

Dd = sysd1.D;

end

%% derivative of the hills terms for the state space

function [Rx,Ry] = hill(p,x,y,m) 

k = p(1);
n = p(2);
KM = p(3);

Rx = (k*y^n)/(KM^n + y^n);
Ry = (k*n*x*(y^(n-1))*KM^n)/(KM^n + y^n)^2;

if m ==1
    Rx = -Rx;
    Ry = -Ry;
end

end