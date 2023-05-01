%% Microbial Fuel Cell Model

% Steady state model
% Represents the end of Phase 1: Creation of model 
% DO NOT MODIFY
clc
clear  
close all 

%% Assumptions/Key Details
% Model based on acetate 
% Anode reaction: (CH2O)2 + 2H20 > 2CO2 + 8H+ + 8e-
% Cathode reaction: O2 + 4e- + 2H20 > 4OH-
% Anode and cathode overpotentials remain fixed


%% Timestep definition
tmax = 600;
d_t=0.01;
t = 0:d_t:tmax;
%% Parameter definition

% Constants

F = 96485.4; % Faraday's constant (Coulombs mol-1))
R = 8.3144; % Gas Constant (J mol-1 K-1))
P = 1; % Pressure (atm)

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


% Surface area of electrode page 2 of Zheng in the experimental section
% 2 Strips of graphite felt 4.5cm by 1cm by 0.5 cm
Ac = 9E-4;

% Cell architecture
U0 = 0.77; % Cell open circuit potential (V)
dm = 1.778E-4; % Thickness of membrane (m)
dcell = 2.2E-2; % Distance between anode and cathode in cell (m)
km = 17; % Electrical conductivity of membrane (ohm-1 m-1)
kaq = 5; % Electrical conductivity of the aqueous solution (ohm-1 m-1)
         % Will need to specify membrane and solution compositions for
         % report

% Reaction rate
% These are temperature dependent
% Will need to account for that when varying temperature
k01 = 0.207; % Rate constant of anode reaction at standard conditions (mol m-2 h-1)
k02 = 3.288E-5; % Rate constant of cathode reaction at standard conditions


% These are also temperature dependent 
alpha = 0.051; % Charge transfer coefficient in the anode
beta = 0.663; % Charge transfer coefficient in the cathode
Kac = 0.592; % Half velocity rate constant for acetate (mol m-3)
Ko2 = 0.004; % Half velocity rate constant for dissolved oxygen (mol m-3)


% Capacitances
CapA = 4E2; % Capacitance of anode (F m-2)
CapC = 5E2; % Capacitance of cathode (F m-2)

%% Matrix creation
% Uses static allocation to reduce compute time compared to dynamic 

% Anode mass balance values
Cac=zeros(1,length(t)); % concentration of acetate (mol m-3)
Cco2 = zeros(1,length(t)); % concentration of dissolved CO2 (mol m-3)
Ch = zeros(1,length(t)); % concentration of H+ (mol m-3)
Cx = zeros(1,length(t)); % concentration of bacteria (mol m-3)


% Cathode mass balance values 
Co2 = zeros(1,length(t)); % concentration of oxygen (mol m-3)
Coh = zeros(1,length(t)); % concentration of OH- ions (mol m-3)
Cm = zeros(1,length(t));  % concentration of cations (mol m-3)

% Reaction rates
% Anode reaction rate
r1 = zeros(1,length(t)); % (mol m-2 h-1)

% Cathode reaction rate
r2 = zeros(1,length(t)); % (mol m-2 h-1)

% Overpotentials
etaA = zeros(1,length(t)); % (V)
etaC = zeros(1,length(t)); % (V)

% Cell voltage
Ucell = zeros(1,length(t)); % (V)



% loop for different current densities
% Data to save - final values inmatrix as I'm interested in the steady
% state values once the system has stabilised 
% Reaction rates
% Concentrations
% Overpotentials

icellMAX = 11.5;

icellM = zeros(1,length(icellMAX));
r1M = zeros(1,length(icellMAX));
r2M = zeros(1,length(icellMAX));
CacM = zeros(1,length(icellMAX));
Cco2M = zeros(1,length(icellMAX));
ChM = zeros(1,length(icellMAX));
CmM = zeros(1,length(icellMAX));
CohM = zeros(1,length(icellMAX));
Co2M =zeros(1,length(icellMAX));
CxM = zeros(1,length(icellMAX));
etaAM = zeros(1,length(icellMAX));
etaCM = zeros(1,length(icellMAX));
UcellM = zeros(1,length(icellMAX));

inc = 0; % Incremental value to let me do smaller increments of the current density

for icell = 0.1:0.1:icellMAX

    Nm = 3600*icell/F;
    
    inc = inc +1;

    icellM(inc) = icell;
    %% Initial Value Assignment

    % Mass balance concentration values
    % As defined by Zheng et al
    Cac(1) = CacIN;
    Cco2(1) = Cco2IN;
    Ch(1) = ChIN;
    Cx(1) = CxIN;
    Co2(1) = Co2IN;
    Coh(1) = CohIN;
    Cm(1) = CmIN;
    

   
    % Anode steady state reaction rate  
    r1(1) = 3600*icell/(8*F);


    % Anode overpotential at steady state
    etaA(1) = R*T/(alpha*F)*log((Qa+Va*Kdec*fx)/(k01*Yac*Am*fx)*((Kac)/(CacIN -r1(1)*(Am/Qa)) +1)); % (V)
                    % Taken from page 7 of Zheng and gives -0.251395962275225
                    % Based on Figure 4 (d) I believe this is correct 
    

    
    % Cathode steady state reaction rate
    
    r2(1) = -3600*icell/(4*F);


    % Cathode overpotential at steady state
    etaC(1) = R*T/(F*(beta-1))*log(r2(1)/-k02*(Ko2+Co2(1))/(Co2(1)));

    
    % Cell Voltage
    
    Ucell(1) = U0 - etaA(1) + etaC(1) -(dm/km + dcell/kaq)*icell;


    %% Equations
    
    for i=1:(length(t)-1)
    
    % Mass balances in anode
    
    Cac(i+1) = Cac(i) + d_t*(Qa*(CacIN - Cac(i)) - Am*r1(i))/Va; % Acetate mass balance
    

    
    
    Cco2(i+1) = Cco2(i) + d_t*(Qa*(Cco2IN - Cco2(i)) + 2*Am*r1(i))/Va; % Dissolved CO2 mass balance
    
    Ch(i+1) = Ch(i) + d_t*(Qa*(ChIN - Ch(i)) + 8*Am*r1(i))/Va; % H+ ions mass balance
    
    Cx(i+1) = Cx(i) + d_t*(Qa*(CxIN-Cx(i))/fx + Am*Yac*r1(i) - Va*Kdec*Cx(i))/Va; % Bacteria mass balance
    
    % Anode overpotential
    
    etaA(i+1) = etaA(i) + d_t*(3600*icell - 8*F*r1(i))/CapA; % Change in anode overpotential
    
    % Reaction Rate in anode
    %r1(i+1) = k01*exp((alpha*F)/(R*T)*etaA(i))*(Cac(i)/(Kac+Cac(i)))*Cx(i); % Anode reaction rate
    
    r1(i+1) = 3600*icell/(8*F);
    
    % Mass balance in cathode
    
    Co2(i+1) = Co2(i) + d_t*(Qc*(Co2IN - Co2(i)) + Am*r2(i))/Vc;
    
    Coh(i+1) = Coh(i) + d_t*(Qc*(CohIN - Coh(i)) - 4*Am*r2(i))/Vc;
    
    Cm(i+1) = Cm(i) + d_t/Vc*(Qc*(CmIN - Cm(i)) + Am*Nm); % When Cm(i) = CmIN nothing happens
    
    % Cathode overpotential
    
    %etaC(i+1) = R*T/F*(beta -1)*log(r2(i)/-k02*(Ko2+Co2(i))/(Co2(i)));

    etaC(i+1) = etaC(i) + d_t*(-3600*icell - 4*F*r2(i))/CapC; % Change in Cathode overpotential
    
    
    % Reaction Rate in cathode
    
    %r2(i+1) = -k02*Co2(i)/(Ko2+Co2(i))*exp((beta-1)*F/(R*T)*etaC(i));
    r2(i+1) = -3600*icell/(4*F);
    
    % Cell voltage
    
    Ucell(i+1) = U0 - etaA(i) + etaC(i) -(dm/km + dcell/kaq)*icell;
    
    end

    r1M(inc) = r1(end);
    r2M(inc) = r2(end);
    CacM(inc) = Cac(end);
    Cco2M(inc) = Cco2(end);
    ChM(inc) = Ch(end);
    CmM(inc) = Cm(end);
    CohM(inc) = Coh(end);
    Co2M(inc) =Co2(end);
    CxM(inc) = Cx(end);
    etaAM(inc) = etaA(end);
    etaCM(inc) = etaC(end);
    UcellM(inc) = Ucell(end);


       
end

%% Calculations for plots and results

% Calculate pH and pOH
phAnode = -log10(ChM/1E3); % Convert to moles per litre
phCathode = 14 - -log10(CohM/1E3);

% Calculate power density
powerDensityM = icellM.*UcellM; % W m-2
                                % This works because the current density is
                                % A m-2 multiplied by voltage = W m-2


%% Plotting
MS = 1.75; % 'MarkerSize' value for plots 


figure(1)
tiledlayout(2,2)
nexttile
plot(icellM,r1M,'LineWidth',1,'Displayname','r1','Marker','o','MarkerSize',MS)
hold on
plot(icellM,r2M,'LineWidth',1,'Displayname','r2','Marker','o','MarkerSize',MS)
hold off
title("Reaction Rates")
grid
grid minor
ylim([-0.12 0.08])
yticks(-0.12:0.04:0.08)
ylabel('Reaction rate (mol m^{-2} h^{-1})','FontWeight','bold')
xlabel('Cell Current Density (A m^{-2})','FontWeight','bold')
legend

nexttile
plot(icellM,CacM,'LineWidth',1,'Displayname','Acetate','Marker','o','MarkerSize',MS)
hold on
plot(icellM,Cco2M,'LineWidth',1,'Displayname','CO2','Marker','o','MarkerSize',MS)
hold on
% plot(icellM,ChM,'LineWidth',1,'Displayname','H^{+} ions','Marker','o','MarkerSize',MS)
% hold on
plot(icellM,CmM,'LineWidth',1,'Displayname','Cations','Marker','o','MarkerSize',MS)
hold on
% plot(icellM,CohM,'LineWidth',1,'Displayname','OH^{-}','Marker','o','MarkerSize',MS)
% hold on
plot(icellM,Co2M,'LineWidth',1,'Displayname','O2','Marker','o','MarkerSize',MS)
hold on
plot(icellM,CxM,'LineWidth',1,'Displayname','Bacteria','Marker','o','MarkerSize',MS)
hold off
grid 
grid minor
legend
title('All concentrations')
ylabel('Concentration (mol m^{-3})','FontWeight','bold')
xlabel('Cell Current Density (A m^{-2})','FontWeight','bold')


nexttile
plot(icellM,phAnode,'LineWidth',1,'Displayname','Anode','Marker','o','MarkerSize',MS)
hold on
plot(icellM,phCathode,'LineWidth',1,'Displayname','Cathode','Marker','o','MarkerSize',MS)
hold off
title("Cell pH")
grid
grid minor
ylabel('pH','FontWeight','bold')
xlabel('Cell Current Density (A m^{-2})','FontWeight','bold')
legend

nexttile
plot(icellM,etaAM,'LineWidth',1,'Displayname','Anode','Marker','o','MarkerSize',MS)
hold on
plot(icellM,etaCM,'LineWidth',1,'Displayname','Cathode','Marker','o','MarkerSize',MS)
hold off
title("Overpotentials")
grid
grid minor
ylabel('Overpotentials (V)','FontWeight','bold')
xlabel('Cell Current Density (A m^{-2})','FontWeight','bold')
legend

%matlab2tikz('myfile.tex');

% figure(2)
% plot(icellM,UcellM,'LineWidth',1,'Displayname','Voltage','Marker','o','MarkerSize',MS)
% title("Cell Voltage")
% grid
% grid minor
% ylabel('Cell Voltage (V)','FontWeight','bold')
% xlabel('Cell Current Density (A m^{-2})','FontWeight','bold')
% legend
% 
% 
% figure(3)
% plot(icellM,powerDensityM,'LineWidth',1,'Displayname','Cell Power Density','Marker','o','MarkerSize',MS)
% title("Power Density")
% grid
% grid minor
% ylabel('Power Density (W m^{-2})','FontWeight','bold')
% xlabel('Cell Current Density (A m^{-2})','FontWeight','bold')
% legend

figure(4)
yyaxis left
plot(icellM,UcellM,'LineWidth',1,'Displayname','Cell Voltage','Marker','o','MarkerSize',MS)
ylabel('Cell Voltage (V)','FontWeight','bold')
yyaxis right
plot(icellM,powerDensityM,'LineWidth',1,'Displayname','Cell Power Density','Marker','o','MarkerSize',MS)
ylabel('Power Density (W m^{-2})','FontWeight','bold')
title("Polarisation and Power Curve")
grid
grid minor
xlabel('Cell Current Density (A m^{-2})','FontWeight','bold')
legend

% figure(1)
% % Top two plots
% tiledlayout(2,2)
% 
% nexttile
% plot(t,r1,'LineWidth',1,'Displayname','r1')
% hold on
% plot(t,r2,'LineWidth',1,'Displayname','r2')
% hold off
% xlabel('Time (h)','FontWeight','bold')
% ylabel('Reaction rate (mol m^{-2} h^{-1})','FontWeight','bold')
% grid
% title('Reaction rates')
% legend
% 
% nexttile
% plot(t,Cac,'k','LineWidth',1)
% xlabel('Time (h)','FontWeight','bold')
% ylabel('Acetate Concentration (mol m^{-3})','FontWeight','bold')
% ylim([0 1.6])
% grid
% title('Acetate concentration')
% 
% nexttile([1 2])
% plot(t,Cac,'LineWidth',1,'Displayname','Acetate conc')
% hold on
% plot(t,Cco2,'LineWidth',1,'Displayname','CO2 conc')
% hold on 
% plot(t,Co2,'LineWidth',1,'Displayname','O2 conc')
% hold on 
% plot(t,Cx,'LineWidth',1,'Displayname','Bacteria conc')
% hold on
% plot(t,Cm,'LineWidth',1,'Displayname','Cation conc')
% xlabel('Time (h)','FontWeight','bold')
% ylabel('Concentration (mol m^{-3})','FontWeight','bold')
% ylim([-0.1 1.6])
% grid
% title('All concentrations')
% legend
% 
% 
% % figure(2)
% % % Top two plots
% % tiledlayout(2,2)
% % 
% % nexttile([1 2])
% % plot(icell,Cac,'k','LineWidth',1)
% % xlabel('Current density (A m^{-2})','FontWeight','bold')
% % ylabel('Acetate Concentration (mol m^{-3})','FontWeight','bold')
% % title('Acetate concentration')
% % 
% % nexttile([1 2])
% % plot(icell,Ucell,'k','LineWidth',1)
% % xlabel('Current density (A m^{-2})','FontWeight','bold')
% % ylabel('Cell Voltage (V)','FontWeight','bold')
% % title('Voltage')

