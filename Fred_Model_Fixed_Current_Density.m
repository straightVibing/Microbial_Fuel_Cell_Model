%% Microbial Fuel Cell Model

% Steady state model

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

% Fixed current density

%icell = 10; % Cell current density (A m-2)
%Nm = 3600*icell/F; % Superficial flux of cations through the membrane (mol m-2 h-1)
                   % Should this be fixed as well?


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


% Current density
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

icellMAX = 12;

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

for icell = 0:1:icellMAX

    Nm = 3600*icell/F;

    icellM(icell+1) = icell;
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
    
    % Mass balance concentration values
    % As defined by me to get over issues with the mass balance
    % If Cac(1) = CacIn then dCac/dt = 0
    % A result of incorrect or poorly ordered equations
    % If I can't fix them then I should find concentrations for specific times
    % and not t = 0
    
    
    % Cac(1) = CacIN;
    % Cco2(1) = 1;
    % Ch(1) = 1;
    % Cx(1) = 1;
    % Co2(1) = Co2IN;
    % Coh(1) = 1;
    % Cm(1) = 1;
    
    % Overpotentials
    
    
    % Reaction rates
    % From Zheng et al
    % Anode reaction rate
    %r1(1) = k01*exp((alpha*F)/(R*T)*etaA(1))*(Cac(1)/(Kac+Cac(1)))*Cx(1);
    
    r1(1) = 3600*icell/(8*F);
    
    etaA(1) = R*T/(alpha*F)*log((Qa+Va*Kdec*fx)/(k01*Yac*Am*fx)*((Kac)/(CacIN -r1(1)*(Am/Qa)) +1)); % (V)
                    % r1 is initially 0 @ t = 0
                    % Taken from page 7 of Zheng and gives -0.251395962275225
                    % Based on Figure 4 (d) I believe this is correct 
    
    etaC(1) = etaA(1); % Figure 4 (d) looks like its the same as etaA(1)
                 % A bit of theory as to why would be good for diss
                 % justification
    
    
    % Cathode reaction rate
    %r2(1) = -k02*Co2(1)/(Ko2+Co2(1))*exp((beta-1)*F/(R*T)*etaC(1));
    
    r2(1) = -3600*icell/(4*F);
    
    % Cell Voltage
    
    Ucell(1) = U0 - etaA(1) + etaC(1) -(dm/km + dcell/kaq)*icell;


    %% Equations
    
    for i=1:(length(t)-1)
    
    % Mass balances in anode
    
    Cac(i+1) = Cac(i) + d_t*(Qa*(CacIN - Cac(i)) - Am*r1(i))/Va; % Acetate mass balance
    
    % Test for weirdness
    
    if isnan(r1(i))
        disp(i)
        break
    end
    
    
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
    
    etaC(i+1) = etaC(i) + d_t*(-3600*icell - 4*F*r2(i))/CapC; % Change in Cathode overpotential
    
    
    % Reaction Rate in cathode
    
    %r2(i+1) = -k02*Co2(i)/(Ko2+Co2(i))*exp((beta-1)*F/(R*T)*etaC(i));
    r2(i+1) = -3600*icell/(4*F);
    
    % Cell voltage
    
    Ucell(i+1) = U0 - etaA(i) + etaC(i) -(dm/km + dcell/kaq)*icell;
    
    end

    r1M(icell+1) = r1(end);
    r2M(icell+1) = r2(end);
    CacM(icell+1) = Cac(end);
    Cco2M(icell+1) = Cco2(end);
    ChM(icell+1) = Ch(end);
    CmM(icell+1) = Cm(end);
    CohM(icell+1) = Coh(end);
    Co2M(icell+1) =Co2(end);
    CxM(icell+1) = Cx(end);
    etaAM(icell+1) = etaA(end);
    etaCM(icell+1) = etaC(end);


       
end




figure(1)
tiledlayout(2,2)
nexttile
plot(icellM,r1M,'LineWidth',1,'Displayname','r1','Marker','o')
hold on
plot(icellM,r2M,'LineWidth',1,'Displayname','r2','Marker','o')
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
plot(icellM,CacM,'LineWidth',1,'Displayname','Acetate','Marker','o')
hold on
plot(icellM,Cco2M,'LineWidth',1,'Displayname','CO2','Marker','o')
hold on
% plot(icellM,ChM,'LineWidth',1,'Displayname','H^{+} ions','Marker','o')
% hold on
plot(icellM,CmM,'LineWidth',1,'Displayname','Cations','Marker','o')
hold on
% plot(icellM,CohM,'LineWidth',1,'Displayname','OH^{-}','Marker','o')
% hold on
plot(icellM,Co2M,'LineWidth',1,'Displayname','O2','Marker','o')
hold on
plot(icellM,CxM,'LineWidth',1,'Displayname','Bacteria','Marker','o')
hold off
grid 
grid minor
legend
title('All concentrations')
ylabel('Concentration (mol m^{-3})','FontWeight','bold')
xlabel('Cell Current Density (A m^{-2})','FontWeight','bold')

% Calculate pH and pOH
phAnode = -log10(ChM/1E3); % Convert to moles per litre
phCathode = -log10(CohM/1E3);

nexttile
plot(icellM,phAnode,'LineWidth',1,'Displayname','Anode','Marker','o')
hold on
plot(icellM,phCathode,'LineWidth',1,'Displayname','Cathode','Marker','o')
hold off
title("Cell pH")
grid
grid minor
ylabel('pH','FontWeight','bold')
xlabel('Cell Current Density (A m^{-2})','FontWeight','bold')
legend

nexttile
plot(icellM,etaAM,'LineWidth',1,'Displayname','Anode','Marker','o')
hold on
plot(icellM,etaCM,'LineWidth',1,'Displayname','Cathode','Marker','o')
hold off
title("Overpotentials")
grid
grid minor
ylabel('Overpotentials (V)','FontWeight','bold')
xlabel('Cell Current Density (A m^{-2})','FontWeight','bold')
legend


%% Plotting
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

