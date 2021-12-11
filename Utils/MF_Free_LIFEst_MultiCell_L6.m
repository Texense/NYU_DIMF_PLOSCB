% MF_Free_LIFEst_MultiCell_L6.m

% This file is part of NYU_DIMF_PLOSCB, the code
% accompanying the paper "A data-informed mean-field
% approach to mapping of cortical parameter landscapes," to
% appear in PLoS Computational Biology
% (https://arxiv.org/abs/2110.12286).
% 
% Copyright (C) 2021 by Zhuo-Cheng Xiao <zx555@nyu.edu>
% 
% This program is free software; you can redistribute it
% and/or modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation;
% either version 2 of the License, or (at your option) any
% later version.
% 
% This program is distributed in the hope that it will be
% useful, but WITHOUT ANY WARRANTY; without even the implied
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
% PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public
% License along with this program; if not, write to the Free
% Software Foundation, Inc., 51 Franklin Street, Fifth
% Floor, Boston, MA 02110-1301 USA.

%% single cell simulation: to collect mean V (and firing rates)
% Input: NeuronType           'e' or 'i'
%        InitialStates        V and all Gs
%        NeuronNum            Number of NeuronSimulated. Generally we don't
%        want this number too large.

% Output: meanV, fr           Collected from all neurons, and time
%         EndStates           V and all Gs, for the next iteration

function [meanV,fr,EndStates] = MF_Free_LIFEst_MultiCell_L6(NeuronType, f_EnI, ...
                                         N_EE,N_EI,N_IE,N_II,...
                                         S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                         lambda_E,S_Elgn,rE_amb,S_amb,...
                                         lambda_I,S_Ilgn,rI_amb,...
                                         S_EL6,S_IL6,rE_L6,rI_L6,...%4 L6 added
                                         tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
                                         rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
                                         gL_E,gL_I,Ve,Vi,InitialStates,NeuronNum)
%% attribute parameters for different types of neurons
if strcmpi(NeuronType,'e')
    N_E = N_EE; N_I = N_EI;
    S_E = S_EE; S_I = S_EI;
    p_Fail = p_EEFail;
    lambda = lambda_E; S_lgn = S_Elgn; r_amb = rE_amb;
    S_L6 = S_EL6; r_L6 = rE_L6;
    gL = gL_E;
    rho_ampa = rhoE_ampa; rho_nmda = rhoE_nmda;
elseif strcmpi(NeuronType,'i')
    N_E = N_IE; N_I = N_II;
    S_E = S_IE; S_I = S_II;
    p_Fail = 0;
    lambda = lambda_I; S_lgn = S_Ilgn; r_amb = rI_amb;
    S_L6 = S_IL6; r_L6 = rI_L6;
    gL = gL_I;
    rho_ampa = rhoI_ampa; rho_nmda = rhoI_nmda;
else 
    disp('***Unrecognized Neuron Type')    
end
rE = f_EnI(1)/1000; rI = f_EnI(2)/1000; % f_EnI in s^-1, but here we use ms^-1
%% Evolve single neurons
T = 1*1e3; % in ms
dt = 0.1; t = 0:dt:T;
SampleProp =0.9; % last half time for meanV

v = zeros(length(t),NeuronNum); % time by rows and neuronInd by columns 
G_gaba_D = zeros(length(t),NeuronNum); G_gaba_R = zeros(length(t),NeuronNum);
G_ampa_D = zeros(length(t),NeuronNum); G_ampa_R = zeros(length(t),NeuronNum);
G_nmda_D = zeros(length(t),NeuronNum); G_nmda_R = zeros(length(t),NeuronNum);
spike = [];
G_I = zeros(length(t),NeuronNum); G_E = zeros(length(t),NeuronNum);

v(1,:) = InitialStates.v;
G_gaba_D(1,:)  = InitialStates.G_gaba_D; G_gaba_R(1,:)  = InitialStates.G_gaba_R;
G_ampa_D(1,:)  = InitialStates.G_ampa_D; G_ampa_R(1,:)  = InitialStates.G_ampa_R;
G_nmda_D(1,:)  = InitialStates.G_nmda_D; G_nmda_R(1,:)  = InitialStates.G_nmda_R;
%rng(100)
% input determination: Assume all Poisson
%rng(100)
p_lgn = 1-exp(-dt*lambda);            Sp_lgn = double(rand(size(v))<=p_lgn);
%rng(101)
p_amb = 1-exp(-dt*r_amb);             Sp_amb = double(rand(size(v))<=p_amb);
%rng(102)
p_EV1 = dt*rE*full(N_E)*(1-p_Fail);   Sp_EV1 = double(rand(size(v))<=p_EV1); % preserving expectation correct
%rng(103)
p_IV1 = dt*rI*full(N_I);              Sp_IV1 = double(rand(size(v))<=p_IV1); % preserving expectation correct
p_L6  = dt*r_L6;                      Sp_L6  = double(rand(size(v))<=p_L6);

%% NOTE!! I am cheating a little here: RefTimer always restart from 0.
RefTimer = zeros(size(InitialStates.v));
for tInd = 1:length(t)-1
     % Taking all previous states, compute new v as is (defer ref decisions later)
     G_I(tInd,:) = 1/(tau_gaba_D-tau_gaba_R) * (G_gaba_D(tInd,:) - G_gaba_R(tInd,:)); % S_EI is included in amplitude of GE_gaba
     G_E(tInd,:) = 1/(tau_ampa_D-tau_ampa_R) * (G_ampa_D(tInd,:) - G_ampa_R(tInd,:)) ...
           + 1/(tau_nmda_D-tau_nmda_R) * (G_nmda_D(tInd,:) - G_nmda_R(tInd,:)); % 
     vv = v(tInd,:) + dt*(-gL*v(tInd,:) - G_E(tInd,:) .*(v(tInd,:)-Ve) - G_I(tInd,:).*(v(tInd,:)-Vi)); %Note: This ensures [nan v] gives [nan vv]
     spike = [spike,tInd*dt*ones(1,sum(vv>=1,'all'))];
     vv(vv>=1) = nan; 
     
     %Now, refrectory neurons get out due to exponetial distributed time    
     RefTimer(isnan(v(tInd,:))) = RefTimer(isnan(v(tInd,:)))+dt;%RefTimer goes up for those in ref period
     vv(RefTimer>=tau_ref) = 0;%if timer reach tau_ref, kick v out of refrectory
     RefTimer(RefTimer>=tau_ref) = 0; %Restart        
     v(tInd+1,:) = vv;     

     % conductances
     G_gaba_R(tInd+1,:) = (G_gaba_R(tInd,:) + S_I*Sp_IV1(tInd,:)) * exp(-dt/tau_gaba_R); 
     G_gaba_D(tInd+1,:) = (G_gaba_D(tInd,:) + S_I*Sp_IV1(tInd,:)) * exp(-dt/tau_gaba_D); 
     G_ampa_R(tInd+1,:) = (G_ampa_R(tInd,:) + S_lgn * Sp_lgn(tInd,:) + S_amb * Sp_amb(tInd,:) + S_E * Sp_EV1(tInd,:) * rho_ampa + S_L6*Sp_L6(tInd,:)*rho_ampa) * exp(-dt/tau_ampa_R);
     G_ampa_D(tInd+1,:) = (G_ampa_D(tInd,:) + S_lgn * Sp_lgn(tInd,:) + S_amb * Sp_amb(tInd,:) + S_E * Sp_EV1(tInd,:) * rho_ampa + S_L6*Sp_L6(tInd,:)*rho_ampa) * exp(-dt/tau_ampa_D);
     G_nmda_R(tInd+1,:) = (G_nmda_R(tInd,:) + S_E * Sp_EV1(tInd,:) * rho_nmda + S_L6*Sp_L6(tInd,:)*rho_nmda) * exp(-dt/tau_nmda_R);
     G_nmda_D(tInd+1,:) = (G_nmda_D(tInd,:) + S_E * Sp_EV1(tInd,:) * rho_nmda + S_L6*Sp_L6(tInd,:)*rho_nmda) * exp(-dt/tau_nmda_D);
end

meanV = nanmean(v(floor(end*(1-SampleProp)):end,:),'all');
fr = length(spike)/(T/1e3)/NeuronNum;
EndStates = struct('v',v(end,:),'G_gaba_D',G_gaba_D(end,:),'G_gaba_R',G_gaba_R(end,:),...
                                'G_ampa_D',G_ampa_D(end,:),'G_ampa_R',G_ampa_R(end,:),...
                                'G_nmda_D',G_nmda_D(end,:),'G_nmda_R',G_nmda_R(end,:));
end