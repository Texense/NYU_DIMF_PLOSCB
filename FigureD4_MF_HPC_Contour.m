% FigureD4_MF_HPC_Contour.m

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

%% Script for HPC:
% All the same as Fig_234 but only:
% Figure D4:
%           We start from the ref point and cancel other relations
% There are 7 pars: SEE SEI SIE SII
%                   SELGN SILGN FIL6 1-7
% ParInd1,2: Choose 2 from 7, others stay ref point
function [] = FigureD4_MF_HPC_Contour(ParInd1, ParInd2)
%% A Rough Estimation Contour for S_EI and S_IE
% first, setup connctivity map
CurrentFolder = pwd;
%FigurePath = [CurrentFolder '/Figures'];
addpath(CurrentFolder)
addpath([CurrentFolder '/Utils'])
addpath([CurrentFolder '/Data'])
addpath([CurrentFolder '/HPCData'])

%% Start ParPool.
cluster = gcp('nocreate');
if isempty(cluster)
parpool("local",[4,48]);
end

%% If now, we start computation below...
% Frst setup network
N_HC = 3;
% Number of E and I neurons
n_E_HC = 54; n_I_HC = 31; % per side of HC
N_E = n_E_HC^2 * N_HC^2; % neuron numbers In all
N_I = n_I_HC^2 * N_HC^2;
% Grid sizes of E and I neurons;
Size_HC = 0.500; % in mm;
Size_E = Size_HC/n_E_HC; Size_I = Size_HC/n_I_HC;
% Projection: SD of distances
SD_E = 0.2/sqrt(2); SD_I = 0.125/sqrt(2);
Dist_LB = 0.36; % ignore the connection probability of dist>0.3mmFreeParUse
% Peak probability of projection
Peak_EE = 0.15; Peak_I = 0.6;

% spatial indexes of E and I neuronsCopy_of_
[NnE.X,NnE.Y] = V1Field_Generation(N_HC,1:N_E,'e');
[NnI.X,NnI.Y] = V1Field_Generation(N_HC,1:N_I,'i');

% determine connections between E, I sparse metrices containing 0 or 1.
% Row_i Column_j means neuron j projects to neuron i add periodic boundary
C_EE = ConnectionMat(N_E,NnE,Size_E,...currentVec(2)^2
    N_E,NnE,Size_E,...
    Peak_EE,SD_E,Dist_LB,1);

C_EI = ConnectionMat(N_E,NnE,Size_E,...
    N_I,NnI,Size_I,...
    Peak_I,SD_I,Dist_LB,0);
% of rE_L6
C_IE = ConnectionMat(N_I,NnI,Size_I,...
    N_E,NnE,Size_E,...
    Peak_I,SD_E,Dist_LB,0);

C_II = ConnectionMat(N_I,NnI,Size_I,...
    N_I,NnI,Size_I,...
    Peak_I,SD_I,Dist_LB,1);

%% Variables and Parameters
% Fixed Pars
p_EEFail = 0.2; 

tau_ampa_R = 0.5; tau_ampa_D = 3;
tau_nmda_R = 2;   tau_nmda_D = 80;
tau_gaba_R = 0.5; tau_gaba_D = 5;
tau_ref = 2; % time unit is ms
%dt = 0.2;
gL_E = 1/20;  Ve = 14/3;  rhoE_ampa = 0.8; rhoE_nmda = 0.2; %S_Elgn = 2*S_EE;
gL_I = 1/15;  Vi = -2/3;  rhoI_ampa = 0.67;rhoI_nmda = 0.33;%S_Ilgn = 0.084;

% of rE_L6
lambda_E = 0.08; % ~16 LGN spike can excite a E neurons. 0.25 spike/ms makes 64 ms for such period.
lambda_I = 0.08;
%rE_amb = 0.72; rI_amb = 0.36;
rE_amb = 0.50; rI_amb = 0.50; S_amb = 0.01;

rE_L6 = 0.25; % rI_L6 to be determined
% Free pars: Fix or Range:in order of SEE SEI SIE SII
%                                     SELGN SILGN FIL6 1-7
GridNum = 160; %160
FreeParFix = {0.024, 0.024*1.51, 0.120*0.147, 0.120,...2.4
              0.048, 0.096,      0.75};
SEERan = [0.88, 1.12];
SEIRan = [0.88, 1.12];  
SIERan = [0.88, 1.12]; 
SIIRan = [0.88, 1.12];
SElgnRan = [0.8 1.2]; 
SIlgnRan = [0.5 1.5];  
rIL6Ran  = [0.4 1.6]; 
FreeParRan = {linspace(FreeParFix{1}*SEERan(1),  FreeParFix{1}*SEERan(2),  GridNum), ...% SEE
              linspace(FreeParFix{2}*SEIRan(1),  FreeParFix{2}*SEIRan(2),  GridNum), ...% SEI
              linspace(FreeParFix{3}*SIERan(1),  FreeParFix{3}*SIERan(2),  GridNum), ...% SIE
              linspace(FreeParFix{4}*SIIRan(1),  FreeParFix{4}*SIIRan(2),  GridNum), ...% SII
              linspace(FreeParFix{5}*SElgnRan(1),FreeParFix{5}*SElgnRan(2),GridNum), ...% SElgn
              linspace(FreeParFix{6}*SIlgnRan(1),FreeParFix{6}*SIlgnRan(2),GridNum), ...% SElgn
              linspace(FreeParFix{7}*rIL6Ran(1), FreeParFix{7}*rIL6Ran(2), GridNum)} ;  %FIL6 
FreeParUse = FreeParFix; 
FreeParUse([ParInd1, ParInd2]) = FreeParRan([ParInd1, ParInd2]);
% SQL6 Put later?
% S_EL6 = 1/3*S_EE;
% S_IL6test = 1/3 * S_IEtest;

%% MF estimation:
%SBound = 3.3; % multipliers of S_EE
Fr_NoFix = zeros(2,GridNum,GridNum );
mV_NoFix = zeros(2,GridNum,GridNum );
Fr_NoFixVar = zeros(2,GridNum,GridNum );
mV_NoFixVar = zeros(2,GridNum,GridNum );
Fr_NoFixTraj = cell(GridNum,GridNum );
mV_NoFixTraj = cell(GridNum,GridNum );

loopCount = zeros(GridNum,GridNum ); % count the number of loops
ConvIndi = logical(loopCount); % converged or not
FailIndi = zeros(size(loopCount));

SampleNum = 50;
StopNum = 500;
h = 1;
SimuT = 20*1e3;


tic
parfor RanInd1 = 1:GridNum   
    for RanInd2 = 1:GridNum
        % Use a new function to read the Pars.....
        [S_EE,S_EI,S_IE,S_II,...
         S_Elgn, S_Ilgn, rI_L6,...
         S_EL6,S_IL6] = AssignParVals(FreeParUse, ParInd1, ParInd2, RanInd1, RanInd2)
        
        tic
        [Fr_NoFixTraj{RanInd1,RanInd2 },mV_NoFixTraj{RanInd1,RanInd2 },...
            loopCount(RanInd1,RanInd2 ),   ConvIndi(RanInd1,RanInd2 ), ...
            FailIndi(RanInd1,RanInd2 )]...
            = MeanFieldEst_BkGd_Indep_StepSize_ref_testL6(C_EE,C_EI,C_IE,C_II,...
            S_EE,S_EI,S_IE,S_II,p_EEFail,...
            lambda_E,S_Elgn,rE_amb,S_amb,...
            lambda_I,S_Ilgn,rI_amb,...
            S_EL6,S_IL6,rE_L6,rI_L6,...
            tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
            rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
            gL_E,gL_I,Ve,Vi,...
            N_HC,n_E_HC,n_I_HC,'End',SampleNum,StopNum,h,SimuT);
        
        Fr_NoFix(:,RanInd1,RanInd2 ) = mean(Fr_NoFixTraj{RanInd1,RanInd2 }(:,end-SampleNum+1:end),2);
        mV_NoFix(:,RanInd1,RanInd2 ) = mean(mV_NoFixTraj{RanInd1,RanInd2 }(:,end-SampleNum:end),2);
        Fr_NoFixVar(:,RanInd1,RanInd2 ) = var(Fr_NoFixTraj{RanInd1,RanInd2 }(:,end-SampleNum:end),0,2);
        mV_NoFixVar(:,RanInd1,RanInd2 ) = var(mV_NoFixTraj{RanInd1,RanInd2 }(:,end-SampleNum:end),0,2);

        
        toc
    end
end
toc

%% Get the Mean of the EarlyTerminate Loops

% save data
ContourData_7D = ws2struct();
% add important info to the end of filename
ParName = {'S_EE','S_EI','S_IE','S_II',...
           'S_Elgn', 'S_Ilgn', 'rI_L6'};
CommentString = sprintf('_%s%s',ParName{ParInd1}, ParName{ParInd2});
save([pwd '/HPCData/FigD4N' CommentString '.mat'],'ContourData_7D')
end

% A function assining free parameters
function [S_EE,S_EI,S_IE,S_II,...
          S_Elgn, S_Ilgn, rI_L6,...
          S_EL6,S_IL6] = AssignParVals(FreeParUse, ParInd1, ParInd2, RanInd1, RanInd2)
      FreeParNow = zeros(size(FreeParUse));
      for ParInd = 1:length(FreeParUse)
          if ParInd~= ParInd1 && ParInd~= ParInd2
              FreeParNow(ParInd) = FreeParUse{ParInd};
          elseif ParInd == ParInd1
              FreeParNow(ParInd1) = FreeParUse{ParInd1}(RanInd1);
          elseif ParInd == ParInd2
              FreeParNow(ParInd2) = FreeParUse{ParInd2}(RanInd2);
          end
      end
S_EE = FreeParNow(1); S_EI = FreeParNow(2);
S_IE = FreeParNow(3); S_II = FreeParNow(4);
S_Elgn= FreeParNow(5); S_Ilgn= FreeParNow(6); 
rI_L6 = FreeParNow(7);
S_EL6 = S_EE/3; S_IL6 = S_IE/3;
end