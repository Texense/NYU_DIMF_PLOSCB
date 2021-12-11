% ConnectionMat.m

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

% ConnectionMat.m is a function generating connetivity matrices
% Input:
%      N_Post,NnPost,Size_Post: number, index projection, and grid size of
%                               presynaptic neurons
%      N_Pre,NnPre,Size_Pre:    Same for postsynaptic neurons
%      Peak_EE:                 Peak connectivity prob
%      SD_E:                    standard daviation of gaussian
%      Dist_UB:                 The upper bound of distances. Connections
%                               between neuron pairs with distances above this are ignored
%      SameCell:                Are pre and post cells are the same or not

% Output:C_EE, the connectivity matrix. Only includes 0 and 1. Diagonal
% ignored for the E-E and I-I matrix.

function C_EE = ConnectionMat(N_Post,NnPost,Size_Post,...
                              N_Pre,NnPre,Size_Pre,...
                              Peak_EE,SD_E,Dist_UB,SameCell)
% get location vectors for all post and pre cells
postNn_Loc = [NnPost.X(1:N_Post);NnPost.Y(1:N_Post)]*Size_Post;
preNn_Loc = [NnPre.X(1:N_Pre);NnPre.Y(1:N_Pre)]*Size_Pre;

% record all effective connection index
Ind1 = zeros(1,N_Post*1000); Ind2 = zeros(1,N_Post*1000); 
tic
subInd = 0;
for postNn = 1:N_Post % for each post synaptic cell, do
    % Distances between, and ignore too-far neurons
    Dist = sqrt(sum((repmat(postNn_Loc(:,postNn),1,N_Pre) - preNn_Loc).^2));
    [EffInd1,EffInd2] = find(Dist<=Dist_UB);
    EffDist = Dist(EffInd2);
    
    % Probability of having a projection
    P_proj = sparse(EffInd1,EffInd2,Peak_EE * exp(-EffDist.^2/(2*SD_E^2)),1,N_Pre); 
    Conn = rand(1,length(P_proj))<P_proj;
    
    % only valid connection indexes are preserved
    [~,ConnInd2] = find(Conn);
    Ind1(subInd+1:subInd+length(ConnInd2)) = postNn*ones(size(ConnInd2));
    Ind2(subInd+1:subInd+length(ConnInd2)) = ConnInd2;
    
    subInd = subInd + length(ConnInd2);
end

% eliminate ineffective indexes, and assign connections
Ind1(Ind1 ==0) = [];
Ind2(Ind2 ==0) = [];
C_EE = sparse(Ind1,Ind2,ones(size(Ind1)));

if SameCell == 1
    C_EE = C_EE-diag(diag(C_EE));
end

toc
end