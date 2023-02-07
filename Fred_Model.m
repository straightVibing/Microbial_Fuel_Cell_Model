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

%% Equations

% Reaction Rate in anode

% Mass balance in anode

% Reaction Rate in cathode

% Mass balance in cathode
f1 = @(cs) rmax*(cs/(ks+cs)) ;

%% Matrix creation and value assignment

%% Initial Values

%% Solving ODEs

%% Calculations

%% Plotting