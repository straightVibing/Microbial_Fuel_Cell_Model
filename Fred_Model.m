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

% Constants

F = 96485.4; % Faraday's constant (Coulombs mol-1))
R = 8.3144; % Gas Constant (J mol-1 K-1))
P = 1; % Pressure (atm)
TREF = 303; % Reference temperature (K)
            % From "A 1D mathematical model..."
Co2REF = (0.21*P)/(R*TREF); % Reference concentration of oxygen
avo =  6.022E23; % Avogadro's constant

% Operational Parameters

T = 303; % Operational temperature (K)
           % Will be varying this later

Va = 5; % Volume of anodic compartment (m3)
Qa = 4; % Volumetric flowrate into the anode (m3 s-1)
Am = 3; % Membrane cross section (m^2)
Vc = 5; % Volume of cathodic compartment (m3)
Qc = 4; % Volumetric flowrate into the cathode (m3 s-1)
CacIN = 2; % Initial concentration of acetate (mol m-3)
Cco2IN = 3; % Initial concentration of dissolved CO2 (mol m-3)
ChIN = 4; % Initial concentration of H+ (mol m-3)
CxIN = 5; % Initial concentration of bacteria (mol m-3)
fx = 10; % Reciprical of washout fraction, Yac (bacterial yield) and Kdec (decay constant) (dimensionless)
Yac = 0.05; % Bacterial yield
Kdec = 8.33E-4; % decay constant (h-1)
Co2IN = 6; % Initial concentration of O2 (mol m-3)
CohIN = 7; % Initial concentration of OH- (mol m-3)
CmIN = 8; % Initial concentration of M+ cations (mol m-3)


% Reaction rate
k01 = 0.207; % Rate constant of anode reaction at standard conditions (mol m-2 h-1)
k02 = 3.288E-5; % Rate constant of cathode reaction at standard conditions

alpha = 0.051; % Charge transfer coefficient in the anode
beta = 0.663; % Charge transfer coefficient in the cathode
Kac = 0.592; % Half velocity rate constant for acetate (mol m-3)
Ko2 = 0.004; % Half velocity rate constant for dissolved oxygen (mol m-3)


% Current density
CapA = 4E2; % Capacitance of anode (F m-2)
CapC = 5E2; % Capacitance of cathode (F m-2)

%alphaC = 0.44; % cathodic transfer coefficient from "A 1D mathematical model" (dimensionless)
%etaC = ; % cathode overpotential from "A 1D mathematical model" (V)
%io2REF = 4.222E-2*exp((73200/R)*(1/353 - 1/T)) ; % Exchange current density of oxygen from "A 1D mathematical model" (A m-2)
%% Matrix creation
% Uses static allocation to reduce compute time

% Anode mass balance values
Cac=zeros(1,length(t)); % concentration of acetate (mol m-3)
Cco2 = zeros(1,length(t)); % concentration of dissolved CO2 (mol m-3)
Ch = zeros(1,length(t)); % concentration of H+ (mol m-3)
Cx = zeros(1,length(t)); % concentration of bacteria (mol m-3)


% Cathode mass balance values 
Co2 = zeros(1,length(t));
Coh = zeros(1,length(t));
Cm = zeros(1,length(t));

% Current density
Nm = zeros(1,length(t)); % Superficial flux of cations
icell = zeros(1,length(t)); % cell current density


% Reaction rates
% Anode reaction rate
r1 =zeros(1,length(t));

% Cathode reaction rate
r2 = zeros(1,length(t));

% Overpotentials 
etaA = zeros(1,length(t)); % Anode overpotential 
etaC = zeros(1,length(t)); % Cathode overpotential


%% Initial Value Assignment

% Mass balance concentration values
Cac(1) = CacIN;
Cco2(1) = Cco2IN;
Ch(1) = ChIN;
Cx(1) = CxIN;
Co2(1) = Co2IN;
Coh(1) = CohIN;
Cm(1) = CmIN;

% Current density 
%icell(1) = io2REF*Co2(1)/Co2REF*exp((alphaC*etaC*F)/(R*T));

%% Equations

for i=1:(length(t)-1)

% Reaction Rate in anode
r1 = 1;

% Mass balances in anode

Cac(i+1) = Cac(i) + d_t*(Qa*(CacIN - Cac(i)) - Am*r1)/Va; % Acetate mass balance

Cco2(i+1) = Cco2(i) + d_t*(Qa*(Cco2IN - Cco2(i)) + 2*Am*r1)/Va; % Dissolved CO2 mass balance

Ch(i+1) = Ch(i) + d_t*(Qa*(ChIN - Ch(i)) + 8*Am*r1)/Va; % H+ ions mass balance

Cx(i+1) = Cx(i) + d_t*(Qa*(CxIN-Cx(i))/fx + Am*Yac*r1 - Va*Kdec*Cx(i))/Va; % Bacteria mass balance


% Reaction Rate in cathode
r2 = 1;

% Flux of cations through the membrane



% Mass balance in cathode

Co2(i+1) = Co2(i) + d_t*(Qc*(Co2IN - Co2(i)) + Am*r2)/Vc;

Coh(i+1) = Coh(i) + d_t*(Qc*(CohIN - Coh(i)) - 4*Am*r2)/Vc;

Cm(i+1) = Cm(i) + d_t*(Qc*(CmIN - Cm(i)) + Am*Nm)/Vc;

end


%% Matrix creation and value assignment

%% Initial Values

%% Solving ODEs

%% Calculations

%% Plotting
