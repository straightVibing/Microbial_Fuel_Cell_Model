%% Microbial Fuel Cell Model

clc
clear  
close all 

%% Assumptions/Key Details
% Model based on acetate 
% Anode reaction: (CH2O)2 + 2H20 > 2CO2 + 8H+ + 8e-
% Cathode reaction: O2 + 4e- + 2H20 > 4OH-


%% Timestep definition
tmax = 40;
d_t=0.1;
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

Va = 5.5E-5; % Volume of anodic compartment (m3)
Qa = 2.25E-5; % Volumetric flowrate into the anode (m3 h-1)
Am = 5E-4; % Membrane cross section (m^2)
Vc = 5.5E-5; % Volume of cathodic compartment (m3)
Qc = 1.11E-3; % Volumetric flowrate into the cathode (m3 h-1)
CacIN = 1.56; % Initial concentration of acetate (mol m-3)
Cco2IN = 0.0; % Initial concentration of dissolved CO2 (mol m-3)
ChIN = 0.0; % Initial concentration of H+ (mol m-3)
CxIN = 0.0; % Initial concentration of bacteria (mol m-3)
fx = 10; % Reciprical of washout fraction, Yac (bacterial yield) and Kdec (decay constant) (dimensionless)
Yac = 0.05; % Bacterial yield
Kdec = 8.33E-4; % decay constant (h-1)
Co2IN = 0.3125; % Initial concentration of O2 (mol m-3)
CohIN = 0.0; % Initial concentration of OH- (mol m-3)
CmIN = 0.0; % Initial concentration of M+ cations (mol m-3)


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

% Cac(1) = CacIN;
% Cco2(1) = 1;
% Ch(1) = 1;
% Cx(1) = 1;
% Co2(1) = Co2IN;
% Coh(1) = 1;
% Cm(1) = 1;

% Overpotentials - Need to find paper values for these 

% etaA(1) = 0.5; % Need actual!
etaA(1) = R*T/(alpha*F)*log((Qa+Va*Kdec*fx)/(k01*Yac*Am*fx)*((Kac)/(CacIN) +1)); % r1 is initially 0
                % Taken from page 7 of Zheng and gives -0.251395962275225
                % Based on Figure 4 (d) I believe this is correct 

etaC(1) = etaA(1); % Figure 4 (d) looks like its the same as etaA(1)

% Current density 
%icell(1) = io2REF*Co2(1)/Co2REF*exp((alphaC*etaC*F)/(R*T));

% Current density
% Nm(1) = Am*Cm(1);
% icell(1) = Nm(1)*F/3600;

icell(1) = F*Vc*Cm(1);
Nm(1) = 3600*icell(1)/F;

% Reaction rates
% Anode reaction rate
r1(1) = k01*exp((alpha*F)/(R*T)*etaA(1))*(Cac(1)/(Kac+Cac(1)))*Cx(1);

%r(1) = 450*icell(1)/F;

% Cathode reaction rate
r2(1) = -k02*Co2(1)/(Ko2+Co2(1))*exp((beta-1)*F/(R*T)*etaC(1));

%r2(1) = -900*icell(1)/F;

%% Equations

for i=1:(length(t)-1)

% Cell density and cation flux    

% Anode overpotential
etaA(i+1) = etaA(i) + d_t*(3600*icell(i)-8*F*r1(i)*1/CapA);

% Reaction Rate in anode
r1(i+1) = k01*exp((alpha*F)/(R*T)*etaA(i))*(Cac(i)/(Kac+Cac(i)))*Cx(i);
%r1(i+1) = 450*icell(i)/F;
% Mass balances in anode

Cac(i+1) = Cac(i) + d_t*(Qa*(CacIN - Cac(i)) - Am*r1(i))/Va; % Acetate mass balance

Cco2(i+1) = Cco2(i) + d_t*(Qa*(Cco2IN - Cco2(i)) + 2*Am*r1(i))/Va; % Dissolved CO2 mass balance

Ch(i+1) = Ch(i) + d_t*(Qa*(ChIN - Ch(i)) + 8*Am*r1(i))/Va; % H+ ions mass balance

Cx(i+1) = Cx(i) + d_t*(Qa*(CxIN-Cx(i))/fx + Am*Yac*r1(i) - Va*Kdec*Cx(i))/Va; % Bacteria mass balance


% Cathode overpotential

etaC(i+1) = etaC(i) + d_t*(-3600*icell(i)-4*F*r2(i)*1/CapC); 

% Reaction Rate in cathode
r2(i+1) = -k02*Co2(i)/(Ko2+Co2(i))*exp((beta-1)*F/(R*T)*etaC(i));
%r2(i+1) = -900*icell(i)/F;



% Mass balance in cathode

Co2(i+1) = Co2(i) + d_t*(Qc*(Co2IN - Co2(i)) + Am*r2(i))/Vc;

Coh(i+1) = Coh(i) + d_t*(Qc*(CohIN - Coh(i)) - 4*Am*r2(i))/Vc;

Cm(i+1) = Cm(i) + d_t*(Qc*(CmIN - Cm(i)) + Am*Nm(i))/Vc; % When Cm(i) = CmIN nothing happens

% Current density
% Nm(i+1) = Am*Cm(i+1);
% icell(i+1) = Nm(i+1)*F/3600;

icell(i+1) = F*Vc*Cm(i);
Nm(i+1) = 3600*icell(i)/F;

end

%% Plotting
figure(1)
% Top two plots
tiledlayout(2,2)

nexttile
plot(t,r1,'LineWidth',1,'Displayname','r1')
hold on
plot(t,r2,'LineWidth',1,'Displayname','r2')
hold off
xlabel('Time (h)','FontWeight','bold')
ylabel('Reaction rate (mol m^{-2} h^{-1})','FontWeight','bold')
title('Reaction rates')
legend

nexttile
plot(t,Cac,'k','LineWidth',1)
xlabel('Time (h)','FontWeight','bold')
ylabel('Acetate Concentration (mol m^{-3})','FontWeight','bold')
title('Acetate concentration')

nexttile([1 2])
plot(t,Cac,'LineWidth',1,'Displayname','Acetate conc')
hold on
plot(t,Cco2,'LineWidth',1,'Displayname','CO2 conc')
hold on 
plot(t,Co2,'LineWidth',1,'Displayname','O2 conc')
hold on 
plot(t,Cx,'LineWidth',1,'Displayname','Bacteria conc')
hold on
plot(t,Cm,'LineWidth',1,'Displayname','Cation conc')
xlabel('Time (h)','FontWeight','bold')
ylabel('Concentration (mol m^{-3})','FontWeight','bold')
title('All concentrations')
legend


figure(2)
plot(icell,Cac,'k','LineWidth',1)
xlabel('Current density (A m^{-2})','FontWeight','bold')
ylabel('Acetate Concentration (mol m^{-3})','FontWeight','bold')
title('Acetate concentration')
