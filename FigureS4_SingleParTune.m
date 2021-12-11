%% Table: When do single-dir parameter variation leave good area?

%% NOTE!!: Due to the limitation of Unix, I have to use ParRat as int
% So should devided by 100!
% Zhuo-Cheng Xiao 03/20/21
function [Fr_NW, mV_NW, BlowUp] = FigureS4_SingleParTune(ParInd,ParRatIn)
RatioNamesAll = {'SEEr',  'SIIr',  'SEIr',  'SIEr',...
                 'SElgnr','SIlgnr','rElgnr','rIlgnr',...
                 'SEL6r', 'SIL6r', 'rEL6r', 'rIL6r',...
                 'Sambr',          'rEambr','rIambr'};
RatioNames =  RatioNamesAll{ParInd};
ParInput = ones(size(RatioNamesAll));
ParRat = ParRatIn/100;
ParInput(ParInd) = ParRat;
disp(fprintf('Tuning %s as %.2f',RatioNames,ParRat))

% Runing it
[Fr_NW, mV_NW, BlowUp] = FigureS4_Running(ParInput(1), ParInput(2), ParInput(3), ParInput(4),...
                                          ParInput(5), ParInput(6), ParInput(7), ParInput(8),...
                                          ParInput(9), ParInput(10),ParInput(11),ParInput(12),...
                                          ParInput(13),             ParInput(14),ParInput(15),...
                                          RatioNames,  ParRat);
                                      
disp(fprintf('fE=%.2f, fI=%.2f',Fr_NW(1),Fr_NW(2)))
end

%% Concerning Specific Parameters, and run the network
% Changing Ratios for SEE   SII   SEI   SIE
%                     SElgn SIlgn rElgn rIlgn
%                     SEL6  SIL6  rEL6  rIL6
%                     SEamb SIamb rEamb rIamb
% r for ratio
% Output: Firing rates of E and I
% Zhuo-Cheng Xiao 03/20/21
%% Note!!: I was trying to incorporate previous ending states, but that makes script too complicated
function [Fr_NW, mV_NW, BlowUp] = FigureS4_Running(SEEr,  SIIr,  SEIr,  SIEr,...
                                                   SElgnr,SIlgnr,rElgnr,rIlgnr,...
                                                   SEL6r, SIL6r, rEL6r, rIL6r,...
                                                   Sambr,        rEambr,rIambr,...
                                                   RatioNames,   ParRat)
%% A Rough Estimation Contour for S_EI and S_IE
% first, setup connctivity map
CurrentFolder = pwd;
%FigurePath = [CurrentFolder '/Figures'];
addpath(CurrentFolder)
addpath([CurrentFolder '/Utils'])
addpath([CurrentFolder '/Data'])

% load saved network architecture to aviod doing it everytime
load('NtWk_Archit.mat')
% pars
S_EEInd = 3;
S_IIInd = 2;
S_ElgnInd = 2;
S_IlgnInd = 2;
rI_L6Ind = 2;

%Fixed parameters
p_EEFail = 0.2; 

tau_ampa_R = 0.5; tau_ampa_D = 3;
tau_nmda_R = 2; tau_nmda_D = 80;
tau_gaba_R = 0.5; tau_gaba_D = 5;
tau_ref = 2; % time unit is ms
%dt = 0.2;
gL_E = 1/20;  Ve = 14/3;  rhoE_ampa = 0.8; rhoE_nmda = 0.2; %S_Elgn = 2*S_EE;
gL_I = 1/15;  Vi = -2/3; %S_Ilgn = 0.084;
rhoI_ampa = 0.67;rhoI_nmda = 0.33;

%% Tweeking Pars
S_EEtest = [0.018 0.021 0.024 0.027 0.030]; 
S_IItest = [0.08  0.12  0.16  0.20];
S_EE = S_EEtest(S_EEInd);
S_II = S_IItest(S_IIInd);

lambda_E = 0.08; % ~16 LGN spike can excite a E neurons. 0.25 spike/ms makes 64 ms for such period.
lambda_I = 0.08;
%rE_amb = 0.72; rI_amb = 0.36;
rE_amb = 0.50; rI_amb = 0.50; S_amb = 0.01;

%L6 input S_IEOneTime = 0.20*S_II;
S_EL6 = 1/3*S_EE; % S_IL6 = 1/3*S_IEOneTime; Now S_IL6 is porp to S_IE
rE_L6 = 0.25; % rI_L6 to be determined

% Unlike MF Contours, we select a specific SEI and SIE:
S_EI = 1.51*S_EE;
S_IE= 0.147*S_II;%*S_EE; I only specify a vecter length here

S_IL6 = 1/3 * S_IE;
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

%% Now Multiply Every tweaking Pars with the ratios. u for use
SEEu = SEEr*S_EE;       SIIu = SIIr*S_II;       SEIu = SEIr*S_EI;         SIEu = SIEr*S_IE;
SElgnu = SElgnr*S_Elgn; SIlgnu = SIlgnr*S_Ilgn; rElgnu = rElgnr*lambda_E; rIlgnu = rIlgnr*lambda_I;
SEL6u = SEL6r*S_EL6;    SIL6u = SIL6r*S_IL6;    rEL6u = rEL6r*rE_L6;      rIL6u = rIL6r*rI_L6;
Sambu = Sambr*S_amb;                            rEambu = rEambr*rE_amb;   rIambu = rIambr*rI_amb;

%% Network Simulation
T = 4010; SimulationT = 2000;
% E-to-E delay time
dt = 0.05;
T_EEDly = 1.0; N_EEDly = floor(T_EEDly/dt);
% Compute Neuron Ind
E_sideInd = floor(0.5*n_E_HC+1):floor(2.5*n_E_HC);
[E_Ind_X,E_Ind_Y] = meshgrid(E_sideInd,E_sideInd);
E_Ind = (reshape(E_Ind_X,size(E_Ind_X,1)*size(E_Ind_X,2),1)-1)*n_E_HC*N_HC + reshape(E_Ind_Y,size(E_Ind_X,1)*size(E_Ind_X,2),1);

I_sideInd = floor(0.5*n_I_HC+1):floor(2.5*n_I_HC);
[I_Ind_X,I_Ind_Y] = meshgrid(I_sideInd,I_sideInd);
I_Ind = (reshape(I_Ind_X,size(I_Ind_X,1)*size(I_Ind_X,2),1)-1)*n_I_HC*N_HC + reshape(I_Ind_Y,size(I_Ind_X,1)*size(I_Ind_X,2),1);

% Creating Sliding Windows
sampleT = 20;
Sliding = 5;
NSlide = floor(sampleT/Sliding);
TWinBounds = 0:Sliding:T;
Wins = [];
Wins(:,1) = TWinBounds(1:end-NSlide);
Wins(:,2) = TWinBounds(NSlide+1:end);

InSs = load('Initials_L6.mat');
RefTimeE = InSs.RefTimeE; VE = InSs.VE; SpE = InSs.SpE; GE_ampa_R = InSs.GE_ampa_R; GE_nmda_R = InSs.GE_nmda_R; GE_gaba_R = InSs.GE_gaba_R;
GE_ampa_D = InSs.GE_ampa_D; GE_nmda_D = InSs.GE_nmda_D; GE_gaba_D = InSs.GE_gaba_D;
RefTimeI = InSs.RefTimeI; VI = InSs.VI; SpI = InSs.SpI; GI_ampa_R = InSs.GI_ampa_R; GI_nmda_R = InSs.GI_nmda_R; GI_gaba_R = InSs.GI_gaba_R;
GI_ampa_D = InSs.GI_ampa_D; GI_nmda_D = InSs.GI_nmda_D; GI_gaba_D = InSs.GI_gaba_D;
EEDlyRcd = sparse(N_E,N_EEDly);

%E_SpHis= []; I_SpHis=[];
E_Sp = []; I_Sp = [];
VE_T = []; VI_T = [];
sampleN = floor(Sliding/dt); % sample each 2 ms
SimulationN = floor(SimulationT/dt); % show and check every 200ms

BlowUp = false;
SampleInd = 1;

% Trace 1. Fr; 2. mV; 3. V distb
NWTrace = struct('FrEs',   zeros(size(Wins(:,1))), 'FrIs',   zeros(size(Wins(:,1))),...
                 'mVEs',   zeros(size(Wins(:,1))), 'mVIs',   zeros(size(Wins(:,1))));
NWTrace.VEDist =  cell(size(Wins(:,1)));
NWTrace.VIDist =  cell(size(Wins(:,1)));
NWTrace.Wins = Wins;

for TimeN = 1:floor(T/dt)
    [oRefTimeE,oVE,oSpE,oGE_ampa_R,oGE_nmda_R,oGE_gaba_R,... % Output
     oGE_ampa_D,oGE_nmda_D,oGE_gaba_D,...
     oRefTimeI,oVI,oSpI,oGI_ampa_R,oGI_nmda_R,oGI_gaba_R,...
     oGI_ampa_D,oGI_nmda_D,oGI_gaba_D,...
     oEEDlyRcd] = ... % A updated N*T Mat recoding the time of kicks taking effect
        V1NetworkUpdate_Ver2_EEDelay(RefTimeE,VE,SpE,GE_ampa_R,GE_nmda_R,GE_gaba_R,... % These are input kept updating
        GE_ampa_D,GE_nmda_D,GE_gaba_D,...
        RefTimeI,VI,SpI,GI_ampa_R,GI_nmda_R,GI_gaba_R,...
        GI_ampa_D,GI_nmda_D,GI_gaba_D,...
        EEDlyRcd,... % A N*T Mat recoding the time of kicks taking effect
        C_EE,C_EI,C_IE,C_II,...
        SEEu,SEIu,SIEu,SIIu,...
        tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,... % time unit is ms
        dt,p_EEFail,...
        gL_E,Ve,SElgnu,rhoE_ampa,rhoE_nmda,...
        gL_I,Vi,SIlgnu,rhoI_ampa,rhoI_nmda,...
        Sambu,rElgnu,rIlgnu,rEambu,rIambu,...
        SEL6u,SIL6u,rEL6u,rIL6u); % L6 Parameters          % The lower/upper bounds for kick waiting time
    % iteration
    RefTimeE = oRefTimeE; VE = oVE;SpE = oSpE;GE_ampa_R = oGE_ampa_R; GE_nmda_R = oGE_nmda_R; GE_gaba_R = oGE_gaba_R;
    GE_ampa_D = oGE_ampa_D; GE_nmda_D = oGE_nmda_D; GE_gaba_D = oGE_gaba_D;
    RefTimeI = oRefTimeI; VI = oVI;SpI = oSpI;GI_ampa_R = oGI_ampa_R; GI_nmda_R = oGI_nmda_R; GI_gaba_R = oGI_gaba_R;
    GI_ampa_D = oGI_ampa_D; GI_nmda_D = oGI_nmda_D; GI_gaba_D = oGI_gaba_D;
    EEDlyRcd = oEEDlyRcd;
    
    % When reaching every 500ms, demostrate and pause
    if mod(TimeN,SimulationN) == 0 && TimeN*dt>=400
        ParameterDisp = ['SEEu = ' num2str(SEEu) ', SEIu = ' num2str(SEIu) ...
            ', SIEu = ' num2str(SIEu) ', SIIu = ' num2str(SIIu) '\n' ...
            'SElgnu = ' num2str(SElgnu) ', rElgnu = ' num2str(rElgnu) ...
            ', SIlgnu = ' num2str(SIlgnu) ', rIlgnu = ' num2str(rIlgnu) '\n' ...
            'Sambu = ' num2str(Sambu) ', rEambu = ' num2str(rEambu) ', rIambu = ' num2str(rIambu)  '\n'];
        
        fprintf(ParameterDisp)
        
        WinSize = SimulationT;
        T_RateWindow = [TimeN*dt-WinSize TimeN*dt];
        E_SpInd = find(E_Sp(:,2)>=T_RateWindow(1) & E_Sp(:,2)<=T_RateWindow(2) & ismember(E_Sp(:,1),E_Ind));
        I_SpInd = find(I_Sp(:,2)>=T_RateWindow(1) & I_Sp(:,2)<=T_RateWindow(2) & ismember(I_Sp(:,1),I_Ind));
        E_Rate = length(E_SpInd)/(WinSize/1000)/length(E_Ind);
        I_Rate = length(I_SpInd)/(WinSize/1000)/length(I_Ind);
        fprintf('FrE = %.2f, FrI = %.2f\n', E_Rate,I_Rate)

    end
    
    % every 20ms time window
    if mod(TimeN,sampleN) == 0 && TimeN >= floor(sampleT/dt)
        WinSize = NWTrace.Wins(SampleInd,2) - NWTrace.Wins(SampleInd,1);
        T_RateWindow = NWTrace.Wins(SampleInd,:);
        E_SpInd = find(E_Sp(:,2)>=T_RateWindow(1) & E_Sp(:,2)<=T_RateWindow(2) & ismember(E_Sp(:,1),E_Ind));
        I_SpInd = find(I_Sp(:,2)>=T_RateWindow(1) & I_Sp(:,2)<=T_RateWindow(2) & ismember(I_Sp(:,1),I_Ind));
        
        NWTrace.FrEs(SampleInd) = length(E_SpInd)/(WinSize/1000)/length(E_Ind);
        NWTrace.FrIs(SampleInd) = length(I_SpInd)/(WinSize/1000)/length(I_Ind);
        NWTrace.mVEs(SampleInd) = nanmean(oVE(E_Ind));
        NWTrace.mVIs(SampleInd) = nanmean(oVI(I_Ind));
        
        NWTrace.VEDist{SampleInd} = oVE(E_Ind);
        NWTrace.VIDist{SampleInd} = oVI(I_Ind);
        
        SampleInd = SampleInd + 1;
    end
    
    E_Sp = [E_Sp;[find(oSpE),ones(size(find(oSpE)))*TimeN*dt]];
    I_Sp = [I_Sp;[find(oSpI),ones(size(find(oSpI)))*TimeN*dt]];
    VE_T = [VE_T;nanmean(oVE(E_Ind))];
    VI_T = [VI_T;nanmean(oVI(I_Ind))];
    
    if sum(isnan(oVE))>0.90*N_E
        BlowUp = true;
        disp('warning!: Network exploded')
        break
    end
    
    % the end of iteration
end

FrData_15 = ws2struct();
% add important info to the end of filename
CommentString = sprintf('_%s_%.2f_BU%d.mat',RatioNames,   ParRat, double(BlowUp));
save([pwd '/HPCData/FigS4ParaTune' CommentString '.mat'],'FrData_15')

% If blowup, at least save the data first...
Fr_NW = [E_Rate;I_Rate];
mV_NW = [nanmean(oVE(E_Ind));nanmean(oVI(I_Ind))];
end