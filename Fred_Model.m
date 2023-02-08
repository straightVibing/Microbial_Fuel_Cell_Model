%% Microbial Fuel Cell Model

clc
clear  
close all 

%% Assumptions/Key Details
% Model based on acetate 
% Anode reaction: (CH2O)2 + 2H20 > 2CO2 + 8H+ + 8e-
% Cathode reaction: O2 + 4e- + 2H20 > 4OH-


%% Timestep definition
tmax = 1000;
d_t=1;
t = 0:d_t:tmax;
%% Parameter definition

Va = 5; % Volume of anodic compartment (m3)
Qa = 4; % Volumetric flowrate into the anode (m3 s-1)
Am = 3; % Membrane cross section (m^2)
CacIN = 2; % Initial concentration of acetate (mol m-3)
Cco2IN = 3; % Initial concentration of dissolved CO2 (mol m-3)
ChIN = 4; % Initial concentration of H+ (mol m-3)
CxIN = 5; % Initial concentration of bacteria (mol m-3)
r1 = 1;
fx = 10; % Reciprical of washout fraction, Yac (bacterial yield) and Kdec (decay constant) (dimensionless)
Yac = 0.05; % Bacterial yield
Kdec = 8.33E-4; % decay constant (h-1)

%% Matrix creation
% Uses static allocation to reduce compute time
Cac=zeros(1,length(t));
Cco2 = zeros(1,length(t));
Ch = zeros(1,length(t));
Cx = zeros(1,length(t));


%% Initial Value Assignment
Cac(1) = CacIN;
Cco2(1) = Cco2IN;
Ch(1) = ChIN;
Cx(1) = CxIN;



%% Equations

for i=1:(length(t)-1)

% Reaction Rate in anode

% Mass balances in anode
Cac(i+1) = Cac(i) + d_t*(Qa*(CacIN - Cac(i)) - Am*r1)/Va; % Acetate mass balance

Cco2(i+1) = Cco2(i) + d_t*(Qa*(Cco2IN - Cco2(i)) + 2*Am*r1)/Va; % Dissolved CO2 mass balance

Ch(i+1) = Ch(i) + d_t*(Qa*(ChIN - Ch(i)) + 8*Am*r1)/Va; % H+ ions mass balance

Cx(i+1) = Cx(i) + d_t*(Qa*(CxIN-Cx(i))/fx + Am*Yac*r1 - Va*Kdec*Cx(i))/Va; % Bacteria mass balance


% Reaction Rate in cathode

% Mass balance in cathode

end


%% Matrix creation and value assignment

%% Initial Values

%% Solving ODEs

%% Calculations

%% Plotting
