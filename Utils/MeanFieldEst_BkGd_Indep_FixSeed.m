% MeanFieldEst_BkGd_Indep_FixSeed.m

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

%% mean-field stuff. We use parameters andto estimate firing rates
%% NOTE! mean Vs are now autonomus from single cell simulations
% Input: C_EE,C_EI,C_IE,C_II connectivity matrices
%        S_EE,S_EI,S_IE,S_II synaptic strength
%        p_EEFail            Synaptic failure prob
%        lambda_E lambda_I   LGN input
%        rE_amb rI_amb       Ambient input
%        S_Elgn S_Ilgn S_amb Synaptic strength of drives
%        gL_E,gL_I           Leaky time constants
%        Ve,Vi               Reversak potentials
%        mVE,mVI             Mean Vs, collected from simulation
% Output:f_EnI               Estimation of firing rates, E and I; A sequences
%        meanVs              mean V of E and I
function [f_EnIOut,meanVs,loop] = MeanFieldEst_BkGd_Indep_FixSeed(C_EE,C_EI,C_IE,C_II,... %4
                                   S_EE,S_EI,S_IE,S_II,p_EEFail,... %5
                                   lambda_E,S_Elgn,rE_amb,S_amb,... %4
                                   lambda_I,S_Ilgn,rI_amb,... %3
                                   tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,... %7
                                   rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,... %4
                                   gL_E,gL_I,Ve,Vi,... %4
                                   N_HC,n_E_HC,n_I_HC,varargin) %3+x
                                   %This line: we are taking spatial center 
% First define fr and mV
f_EnIOut = [];

% Get connectivity for the center HC
E_sideInd = floor(1*n_E_HC+1):2*n_E_HC;
[E_Ind_X,E_Ind_Y] = meshgrid(E_sideInd,E_sideInd);
E_Ind = (reshape(E_Ind_X,size(E_Ind_X,1)*size(E_Ind_X,2),1)-1)*n_E_HC*N_HC + reshape(E_Ind_Y,size(E_Ind_X,1)*size(E_Ind_X,2),1);

I_sideInd = floor(1*n_I_HC+1):2*n_I_HC;
[I_Ind_X,I_Ind_Y] = meshgrid(I_sideInd,I_sideInd);
I_Ind = (reshape(I_Ind_X,size(I_Ind_X,1)*size(I_Ind_X,2),1)-1)*n_I_HC*N_HC + reshape(I_Ind_Y,size(I_Ind_X,1)*size(I_Ind_X,2),1);

% Take averaged number of input neurons 
% NOTE! Only picking up the middle part
N_EE = mean(sum(C_EE(E_Ind,:),2)); N_EI = mean(sum(C_EI(E_Ind,:),2)); 
N_IE = mean(sum(C_IE(I_Ind,:),2)); N_II = mean(sum(C_II(I_Ind,:),2));
% initialize with a reasonable guess     
mVE = 0.57; mVI = 0.67; 
%mVE =  0; mVI = 0; 
meanVs = [mVE;mVI];
mVEpre  = -1; mVIpre = -1;
f_EnI0 = [0;0]; f_EnIpre = [-10;-10]; % start with an impossible value

loop = 0;
while( norm([mVEpre;mVIpre] - [mVE;mVI]) > 0.01 || norm(f_EnIpre - f_EnI0)>0.2)
mVEpre = mVE; mVIpre = mVI;
f_EnIpre = f_EnI0;
f_EnI0 = MeanFieldEst_BkGd(N_EE,N_EI,N_IE,N_II,...
                                   S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                   lambda_E,S_Elgn,rE_amb,S_amb,...
                                   lambda_I,S_Ilgn,rI_amb,...
                                   gL_E,gL_I,Ve,Vi,mVE,mVI);

f_EnI0 = max([f_EnI0,[0;0]],[],2);
%f_EnI = zeros(size(f_EnI0)); % from simulation
% simulate one neuron                              
[mVE,~] = MEanFieldEst_SingleCell('e', f_EnI0, ...
                                         N_EE,N_EI,N_IE,N_II,...
                                         S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                         lambda_E,S_Elgn,rE_amb,S_amb,...
                                         lambda_I,S_Ilgn,rI_amb,...
                                         tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
                                         rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
                                         gL_E,gL_I,Ve,Vi);
[mVI,~] = MEanFieldEst_SingleCell('i', f_EnI0, ...
                                         N_EE,N_EI,N_IE,N_II,...
                                         S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                         lambda_E,S_Elgn,rE_amb,S_amb,...
                                         lambda_I,S_Ilgn,rI_amb,...
                                         tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
                                         rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
                                         gL_E,gL_I,Ve,Vi);
loop = loop+1;

f_EnIOut = [f_EnIOut,f_EnI0];
meanVs = [meanVs, [mVE;mVI]];
end
if nargin > 34 % Specify: I don't need the trajectory but the final ones
    f_EnIOut = f_EnIOut(:,end);
    meanVs = meanVs(:,end);
end

end



%% mean-field stuff. We use parameters and mean Vs to estimate firing rates
% Input: C_EE,C_EI,C_IE,C_II connectivity matrices
%        N_EE,N_EI,N_IE,N_II number of presynaptic neurons
%        S_EE,S_EI,S_IE,S_II synaptic strength
%        p_EEFail            Synaptic failure prob
%        lambda_E lambda_I   LGN input
%        rE_amb rI_amb       Ambient input
%        S_Elgn S_Ilgn S_amb Synaptic strength of drives
%        gL_E,gL_I           Leaky time constants
%        Ve,Vi               Reversak potentials
%        mVE,mVI             Mean Vs, collected from simulation
% Output:f_EnI               Estimation of firing rates, E and I

function f_EnI = MeanFieldEst_BkGd(N_EE,N_EI,N_IE,N_II,...
                                   S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                   lambda_E,S_Elgn,rE_amb,S_amb,...
                                   lambda_I,S_Ilgn,rI_amb,...
                                   gL_E,gL_I,Ve,Vi,mVE,mVI)
                               
% E_sideInd = floor(1*n_E_HC+1):2*n_E_HC;
% [E_Ind_X,E_Ind_Y] = meshgrid(E_sideInd,E_sideInd);
% E_Ind = (reshape(E_Ind_X,size(E_Ind_X,1)*size(E_Ind_X,2),1)-1)*n_E_HC*N_HC + reshape(E_Ind_Y,size(E_Ind_X,1)*size(E_Ind_X,2),1);
% 
% I_sideInd = floor(1*n_I_HC+1):2*n_I_HC;
% [I_Ind_X,I_Ind_Y] = meshgrid(I_sideInd,I_sideInd);
% I_Ind = (reshape(I_Ind_X,size(I_Ind_X,1)*size(I_Ind_X,2),1)-1)*n_I_HC*N_HC + reshape(I_Ind_Y,size(I_Ind_X,1)*size(I_Ind_X,2),1);
% 
% % Take averaged number of input neurons 
% % NOTE! Only picking up the middle part
% N_EE = mean(sum(C_EE(E_Ind,:),2)); N_EI = mean(sum(C_EI(E_Ind,:),2)); 
% N_IE = mean(sum(C_IE(I_Ind,:),2)); N_II = mean(sum(C_II(I_Ind,:),2));

% Matrix of firing rates
MatEI = [N_EE*S_EE*(Ve-mVE)*(1-p_EEFail), N_EI*S_EI*(Vi-mVE);
         N_IE*S_IE*(Ve-mVI),              N_II*S_II*(Vi-mVI)];

% External Drive and leaky current. Note that simulaton is by ms, so need *1000
VecInput = [(Ve-mVE) * (lambda_E*S_Elgn + rE_amb*S_amb);
            (Ve-mVI) * (lambda_I*S_Ilgn + rI_amb*S_amb)]*1000;
Leak = -[gL_E*mVE;
        gL_I*mVI]*1000; % leak has opposite signs of mean Vs
    
f_EnI = -(MatEI-eye(2))^-1*(Leak+VecInput);

end
%% single cell simulation: to collect mean V (and firing rates)


function [meanV,fr] = MEanFieldEst_SingleCell(NeuronType, f_EnI, ...
                                         N_EE,N_EI,N_IE,N_II,...
                                         S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                         lambda_E,S_Elgn,rE_amb,S_amb,...
                                         lambda_I,S_Ilgn,rI_amb,...
                                         tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
                                         rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
                                         gL_E,gL_I,Ve,Vi)
%% attribute parameters for different types of neurons
if strcmpi(NeuronType,'e')
    N_E = N_EE; N_I = N_EI;
    S_E = S_EE; S_I = S_EI;
    p_Fail = p_EEFail;
    lambda = lambda_E; S_lgn = S_Elgn; r_amb = rE_amb;
    gL = gL_E;
    rho_ampa = rhoE_ampa; rho_nmda = rhoE_nmda;
elseif strcmpi(NeuronType,'i')
    N_E = N_IE; N_I = N_II;
    S_E = S_IE; S_I = S_II;
    p_Fail = 0;
    lambda = lambda_I; S_lgn = S_Ilgn; r_amb = rI_amb;
    gL = gL_I;
    rho_ampa = rhoI_ampa; rho_nmda = rhoI_nmda;
else 
    disp('***Unrecognized Neuron Type')    
end
rE = f_EnI(1)/1000; rI = f_EnI(2)/1000; % f_EnI in s^-1, but here we use ms^-1
%% Evolve single neurons
T = 10000; % in ms
dt = 0.1; t = 0:dt:T;
SampleProp = 1/2; % last half time for meanV

v = zeros(size(t)); 
G_gaba_D = zeros(size(t)); G_gaba_R = zeros(size(t));
G_ampa_D = zeros(size(t)); G_ampa_R = zeros(size(t));
G_nmda_D = zeros(size(t)); G_nmda_R = zeros(size(t));
spike = [];

%rng(100)
% input determination: Assume all Poisson
rng(200)
p_lgn = 1-exp(-dt*lambda);            Sp_lgn = double(rand(size(t))<=p_lgn);
rng(201)
p_amb = 1-exp(-dt*r_amb);             Sp_amb = double(rand(size(t))<=p_amb);
rng(202)
p_EV1 = dt*rE*full(N_E)*(1-p_Fail); Sp_EV1 = double(rand(size(t))<=p_EV1);
rng(203)
p_IV1 = dt*rI*full(N_I);            Sp_IV1 = double(rand(size(t))<=p_IV1);
RefTimer = 0;
for tInd = 1:length(t)-1
    % Firstly, refrectory neurons get out due to exponetial distributed time
     if isnan(v(tInd))
        RefTimer = RefTimer+dt; %RefTimer goes up
        if RefTimer>=tau_ref    %if timer reach tau_ref, kick v out of refrectory
            v(tInd+1) = 0;
            RefTimer = 0;
        else
            v(tInd+1) = nan;
        end
     else
         G_I = 1/(tau_gaba_D-tau_gaba_R) * (G_gaba_D(tInd) - G_gaba_R(tInd)); % S_EI is included in amplitude of GE_gaba
         G_E = 1/(tau_ampa_D-tau_ampa_R) * (G_ampa_D(tInd) - G_ampa_R(tInd)) ...
             + 1/(tau_nmda_D-tau_nmda_R) * (G_nmda_D(tInd) - G_nmda_R(tInd)); % 
         vv = v(tInd) + dt*(-gL*v(tInd) - G_E.*(v(tInd)-Ve) - G_I.*(v(tInd)-Vi)); 
         if vv >= 1
             v(tInd+1) = nan;
             spike = [spike,t(tInd)];
         else
             v(tInd+1) = vv;
         end
     end

     % conductances
     G_gaba_R(tInd+1) = (G_gaba_R(tInd) + S_I*Sp_IV1(tInd)) * exp(-dt/tau_gaba_R); 
     G_gaba_D(tInd+1) = (G_gaba_D(tInd) + S_I*Sp_IV1(tInd)) * exp(-dt/tau_gaba_D); 
     G_ampa_R(tInd+1) = (G_ampa_R(tInd) + S_lgn * Sp_lgn(tInd) + S_amb * Sp_amb(tInd) + S_E * Sp_EV1(tInd) * rho_ampa) * exp(-dt/tau_ampa_R);
     G_ampa_D(tInd+1) = (G_ampa_D(tInd) + S_lgn * Sp_lgn(tInd) + S_amb * Sp_amb(tInd) + S_E * Sp_EV1(tInd) * rho_ampa) * exp(-dt/tau_ampa_D);
     G_nmda_R(tInd+1) = (G_nmda_R(tInd) + S_E * Sp_EV1(tInd) * rho_nmda) * exp(-dt/tau_nmda_R);
     G_nmda_D(tInd+1) = (G_nmda_D(tInd) + S_E * Sp_EV1(tInd) * rho_nmda) * exp(-dt/tau_nmda_D);
end
meanV = nanmean(v(floor(end*SampleProp):end));
fr = length(find(spike>SampleProp*T))/(T*SampleProp/1000);
end