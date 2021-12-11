% Figure234_MF_7DHPC_Contour.m

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
% Run each panel for selected pars
% S_EEtest = [0.018 0.021 0.024 0.027 0.030]; 
% S_IItest = [0.08  0.12  0.16  0.20];
% S_Elgntest = [1.5 2 2.5 3.0]*S_EE;
% S_Ilgn_Mtp = [1.5 2 2.5 3]; % of S_Elgn
% rI_L6_Mtp  = [1.5 3 4.5 6]; % of rE_L6
function [] = Figure234_MF_7DHPC_Contour(S_ElgnInd,S_IlgnInd,rI_L6Ind,...
                                       S_EEInd,S_IIInd)
%% A Rough Estimation Contour for S_EI and S_IE
% first, setup connctivity map
CurrentFolder = pwd;
%FigurePath = [CurrentFolder '/Figures'];
addpath(CurrentFolder)
addpath([CurrentFolder '/Utils'])
addpath([CurrentFolder '/Data'])
addpath([CurrentFolder '/HPCData'])

%% See if data file exists. If true return.
S_EEtest = [0.018 0.021 0.024 0.027 0.030]; 
S_IItest = [0.08  0.12  0.16  0.20];
S_EE = S_EEtest(S_EEInd);
S_II = S_IItest(S_IIInd);

CommentString = sprintf('_S_EE=%.3f_S_II=%.2f_S_ElgnInd=%d_S_IlgnInd%d_rI_L6Ind%d',S_EE,S_II,S_ElgnInd,S_IlgnInd,rI_L6Ind);
SearchFile = [pwd '/HPCData/FigContourL' CommentString '.mat'];
if isfile(SearchFile)
    disp('Data File exists. Return Now')
    return
end

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
Dist_LB = 0.36; % ignore the connection probability of dist>0.3mm
% Peak probability of projection
Peak_EE = 0.15; Peak_I = 0.6;

% spatial indexes of E and I neurons
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

C_IE = ConnectionMat(N_I,NnI,Size_I,...
    N_E,NnE,Size_E,...
    Peak_I,SD_E,Dist_LB,0);

C_II = ConnectionMat(N_I,NnI,Size_I,...
    N_I,NnI,Size_I,...
    Peak_I,SD_I,Dist_LB,1);

%% Variables and Parameters
%RefTimeE = zeros(N_E,1); VE = 0.5*rand(N_E,1)-0.5; SpE = sparse(N_E,1);
%GE_ampa_R = zeros(N_E,1); GE_nmda_R = zeros(N_E,1); GE_gaba_R =
%zeros(N_E,1); GE_ampa_D = zeros(N_E,1); GE_nmda_D = zeros(N_E,1);
%GE_gaba_D = zeros(N_E,1); RefTimeI = zeros(N_I,1); VI =
%1.5*rand(N_I,1)-0.5; SpI = sparse(N_I,1); GI_ampa_R = zeros(N_I,1);
%GI_nmda_R = zeros(N_I,1); GI_gaba_R = zeros(N_I,1); GI_ampa_D =
%zeros(N_I,1); GI_nmda_D = zeros(N_I,1); GI_gaba_D = zeros(N_I,1);
%load('Initials.mat')
%parameters
% S_EEtest = [0.018 0.021 0.024 0.027 0.030]; 
% S_IItest = [0.08  0.12  0.16  0.20];
% S_EE = S_EEtest(S_EEInd);
% S_II = S_IItest(S_IIInd);
p_EEFail = 0.2; S_amb = 0.01;

tau_ampa_R = 0.5; tau_ampa_D = 3;
tau_nmda_R = 2; tau_nmda_D = 80;
tau_gaba_R = 0.5; tau_gaba_D = 5;
tau_ref = 2; % time unit is ms
%dt = 0.2;
gL_E = 1/20;  Ve = 14/3;  rhoE_ampa = 0.8; rhoE_nmda = 0.2; %S_Elgn = 2*S_EE;
gL_I = 1/15;  Vi = -2/3; %S_Ilgn = 0.084;
rhoI_ampa = 0.67;rhoI_nmda = 0.33;

lambda_E = 0.08; % ~16 LGN spike can excite a E neurons. 0.25 spike/ms makes 64 ms for such period.
lambda_I = 0.08;
%rE_amb = 0.72; rI_amb = 0.36;
rE_amb = 0.50; rI_amb = 0.50;

%L6 input S_IEOneTime = 0.20*S_II;
S_EL6 = 1/3*S_EE; % S_IL6 = 1/3*S_IEOneTime; Now S_IL6 is porp to S_IE
rE_L6 = 0.25; % rI_L6 to be determined

% Replace S_EI by testing values
GridNum1 = 160; %160
GridNum2 = 160; %160
S_EI_Mtp = [0.9, 2.4]; % of S_EE
S_IE_Mtp = [0.1, 0.25]; % of S_II
S_EItest = linspace(S_EI_Mtp(1),S_EI_Mtp(2),GridNum1)*S_EE;
S_IEtest = linspace(S_IE_Mtp(1),S_IE_Mtp(2),GridNum2)*S_II;%*S_EE; I only specify a vecter length here

S_IL6test = 1/3 * S_IEtest;
% First determine S_Elgn
S_Elgntest = [1.5 2 2.5 3.0]*S_EE;
S_Elgn = S_Elgntest(S_ElgnInd);

% Panel: Two proportions
S_Ilgn_Mtp = [1.5 2 2.5 3]; % of S_Elgn
rI_L6_Mtp  = [1.5 3 4.5 6]; % of rE_L6
S_Ilgntest = S_Ilgn_Mtp * S_Elgn;
rI_L6test  = rI_L6_Mtp * rE_L6;

S_Ilgn = S_Ilgntest(S_IlgnInd);
rI_L6 = rI_L6test(rI_L6Ind);

% Add lines boundaries
LineL1 = polyfit([0.1  0.2 ],[0.8 0.8],1); % S_IEMtp first, second S_EIMtp. Those numbers are multipliers of S_II and S_EE
LineL2 = polyfit([0.06 0.28],[0.8 0.4],1);
LineU1 = polyfit([0.1  0.3 ],[2.5 2.5  ],1); % LineU1 = polyfit([0.1  0.3 ],[2.5 0.8],1);

% creat a 10-hr parallel
% cluster = gcp('nocreate');
% if isempty(cluster)
%     cluster = parpool([4 64]);
% %    cluster.IdleTimeout = 1200;
% end

%% MF estimation:
%SBound = 3.3; % multipliers of S_EE
Fr_NoFix = zeros(2,length(S_EItest),length(S_IEtest) );
mV_NoFix = zeros(2,length(S_EItest),length(S_IEtest) );
Fr_NoFixVar = zeros(2,length(S_EItest),length(S_IEtest) );
mV_NoFixVar = zeros(2,length(S_EItest),length(S_IEtest) );
Fr_NoFixTraj = cell(length(S_EItest),length(S_IEtest) );
mV_NoFixTraj = cell(length(S_EItest),length(S_IEtest) );

loopCount = zeros(length(S_EItest),length(S_IEtest) ); % count the number of loops
ConvIndi = logical(loopCount); % converged or not
FailIndi = zeros(size(loopCount));

SampleNum = 50;
StopNum = 300;
h = 1;
SimuT = 20*1e3;
aa = floor(length(S_IEtest)); % Matlab always fail to directly see this as a whole number!!!
%bb = floor(length(S_EItest));

tic

% S_IlgnInd = ceil(PanelInd/S_IIInd);
% rI_L6Ind  = mod(PanelInd,S_IIInd);
% rI_L6Ind(rI_L6Ind==0) = S_IIInd;

parfor S_EIInd = 1:length(S_EItest)
    S_EI = S_EItest(S_EIInd);
    
    
    for S_IEInd = 1:aa
        S_IE  = S_IEtest(S_IEInd);
        S_IL6 = S_IL6test(S_IEInd);
        %cut some redundant regime out of our interest
        if (S_EI/S_EE<=S_IE/S_II*LineL1(1)+LineL1(2) || S_EI/S_EE<=S_IE/S_II*LineL2(1)+LineL2(2) )
            disp(['S_IE = ' num2str(S_IE/S_II,'%.3f') '*S_II, S_EI = ' num2str(S_EI/S_EE,'%.3f') '*S_EE; Fr may be too high, break...'])
            continue
        end
        if (S_EI/S_EE>=S_IE/S_II*LineU1(1)+LineU1(2) )
            disp(['S_IE = ' num2str(S_IE/S_II,'%.3f') '*S_II, S_EI = ' num2str(S_EI/S_EE,'%.3f') '*S_EE; Fr may be too low, break...'])
            continue
        end
        
        tic
        [Fr_NoFixTraj{S_EIInd,S_IEInd },mV_NoFixTraj{S_EIInd,S_IEInd },...
            loopCount(S_EIInd,S_IEInd ),   ConvIndi(S_EIInd,S_IEInd ), ...
            FailIndi(S_EIInd,S_IEInd )]...
            = MeanFieldEst_BkGd_Indep_StepSize_ref_testL6(C_EE,C_EI,C_IE,C_II,...
            S_EE,S_EI,S_IE,S_II,p_EEFail,...
            lambda_E,S_Elgn,rE_amb,S_amb,...
            lambda_I,S_Ilgn,rI_amb,...
            S_EL6,S_IL6,rE_L6,rI_L6,...
            tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
            rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
            gL_E,gL_I,Ve,Vi,...
            N_HC,n_E_HC,n_I_HC,'End',SampleNum,StopNum,h,SimuT);
        
        toc
    end
end
toc

% for PanelInd = 1:PanelNum1*PanelNum2
%     S_IlgnInd = ceil(PanelInd/PanelNum2);
%     rI_L6Ind  = mod(PanelInd,PanelNum2);
%     rI_L6Ind(rI_L6Ind==0) = PanelNum2;
    
%     S_Ilgn = S_Ilgntest(S_IlgnInd);
%     rI_L6 = rI_L6test(rI_L6Ind);
    for S_EIInd = 1:length(S_EItest)
        S_EI = S_EItest(S_EIInd);
        
        for S_IEInd = 1:length(S_IEtest)
            S_IE = S_IEtest(S_IEInd);
            % cut some redundant regime out of our interest
            if (S_EI/S_EE<=S_IE/S_II*LineL1(1)+LineL1(2) || S_EI/S_EE<=S_IE/S_II*LineL2(1)+LineL2(2) )
                %disp(['S_IE = ' num2str(S_IE/S_II,'%.2f') '*S_II, S_EI = '
                %num2str(S_EI/S_EE,'%.2f') '*S_EE; Fr may be too high,
                %break...'])
                continue
            end
            if (S_EI/S_EE>=S_IE/S_II*LineU1(1)+LineU1(2) )
                %disp(['S_IE = ' num2str(S_IE/S_II,'%.2f') '*S_II, S_EI = '
                %num2str(S_EI/S_EE,'%.2f') '*S_EE; Fr may be too low,
                %break...'])
                continue
            end
            Fr_NoFix(:,S_EIInd,S_IEInd ) = mean(Fr_NoFixTraj{S_EIInd,S_IEInd }(:,end-SampleNum+1:end),2);
            mV_NoFix(:,S_EIInd,S_IEInd ) = mean(mV_NoFixTraj{S_EIInd,S_IEInd }(:,end-SampleNum:end),2);
            Fr_NoFixVar(:,S_EIInd,S_IEInd ) = var(Fr_NoFixTraj{S_EIInd,S_IEInd }(:,end-SampleNum:end),0,2);
            mV_NoFixVar(:,S_EIInd,S_IEInd ) = var(mV_NoFixTraj{S_EIInd,S_IEInd }(:,end-SampleNum:end),0,2);
        end
    end
%end
% save data
%Trajs = struct('Fr_NoFixTraj', Fr_NoFixTraj, 'mV_NoFixTraj',mV_NoFixTraj);
ContourData_7D = ws2struct();
% add important info to the end of filename
CommentString = sprintf('_S_EE=%.3f_S_II=%.2f_S_ElgnInd=%d_S_IlgnInd%d_rI_L6Ind%d',S_EE,S_II,S_ElgnInd,S_IlgnInd,rI_L6Ind);
save([pwd '/HPCData/FigContourL' CommentString '.mat'],'ContourData_7D')
%% Contour maps
% Fr_Plot = Fr_NoFix; S_IEBound = 8e-3; % DeleteInd =
% false(size(Fr_Plot));SampleProp % DeleteInd(:,:,S_IEtest<S_IEBound) =
% true; Fr_Plot(Fr_NoFix <=eps) = nan;
%
% h(1) = figure('Name','Rough Countour S_EI S_IE'); subplot 121
% imagesc(S_IEtest,S_EItest,squeeze(Fr_Plot(1,:,:))) xlabel('S_{IE}');
% ylabel('S_{EI}'); set(gca,'YDir','Normal') colorbar caxis([0 15]) axis
% square hold on [C1,h1]=
% contour(S_IEtest,S_EItest,squeeze(Fr_Plot(1,:,:)),[2
% 5],'ShowText','on','color','r'); clabel(C1,h1,'FontSize',10,'Color','k')
% axis square hold off %axis([min(S_EItest) max(S_EItest) min(S_IEtest)
% max(S_IEtest)])
%
% subplot 122 imagesc(S_IEtest,S_EItest,squeeze(Fr_Plot(2,:,:)))
% xlabel('S_{IE}'); ylabel('S_{EI}'); set(gca,'YDir','Normal') colorbar
% caxis([0 45]) axis square hold on [C2,h2]=
% contour(S_IEtest,S_EItest,squeeze(Fr_Plot(2,:,:)),[7
% 15],'ShowText','on','color','b'); clabel(C2,h2,'FontSize',10,'Color','k')
% axis square hold off %axis([min(S_EItest) max(S_EItest) min(S_IEtest)
% max(S_IEtest)])
%
% h(2) = figure('Name','Convergence'); imagesc(S_IEtest,S_EItest,ConvIndi)
% xlabel('S_{IE}'); ylabel('S_{EI}'); set(gca,'YDir','Normal')
%
% % h(3) = figure('Name','Trajectories-Last'); % subplot 131 %
% plot_dir(Fr_NoFixTraj{9,13}(1,end-100:end)',Fr_NoFixTraj{9,13}(2,end-100:end)');
% % xlabel('fE');ylabel('fI') % % title('last 100 iterations') % % subplot
% 132 %
% plot_dir(Fr_NoFixTraj{9,13}(1,end-50:end)',Fr_NoFixTraj{9,13}(2,end-50:end)');
% % xlabel('fE');ylabel('fI') % % title('last 50 iterations') % % subplot
% 133 %
% plot_dir(Fr_NoFixTraj{9,13}(1,end-15:end)',Fr_NoFixTraj{9,13}(2,end-15:end)');
% % xlabel('fE');ylabel('fI') % % title('last 15 iterations')
%
%
% saveas(h,fullfile(FigurePath, ['ContourFigs - S_EE=' num2str(S_EE)
% CommentString]),'fig')
end