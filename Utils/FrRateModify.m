% FrRateModify.m

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

%% Modify rate map
% CurrentFrE: The rough ratemap
% ConvMat: Convolution kernel
% CurrentBlowupModi: Blowup mat. Should be same dim of CurrentFrE. We get
% rid of where blowup is true
% CurrentsteadyModi: Steady mat. Should be same dim of CurrentFrE. We get
% rid of where steady is false.

% Output: FringMap
function FringMap = FrRateModify(CurrentFrE, ConvMat,CurrentBlowupModi,CurrentsteadyModi)

CurrentFrE(CurrentFrE<eps) = nan; CurrentFrE(CurrentBlowupModi) = nan;
CurrentFrE(~CurrentsteadyModi) = nan;

CurrentFrE = conv2(CurrentFrE,ConvMat);
[a,b] = size(ConvMat);

% Cut out edges after convolution
% Cutting out 'a-1' rows in all
if a>1
    a_half1 = floor((a-1)/2);
    a_half2 = (a-1) - a_half1;
    CurrentFrE(end-a+2:end-a_half1,:) = nan;
    CurrentFrE(end-a_half1+1:end,:) = [];
    CurrentFrE(a_half2+1:a-1,  :) = nan;
    CurrentFrE(1:a_half2,      :) = [];
end

if b>1
    b_half1 = floor((b-1)/2);
    b_half2 = (b-1) - a_half1;    
    CurrentFrE(:,end-b+2:end-b_half1) = nan;
    CurrentFrE(:,end-b_half1+1:end) = [];
    CurrentFrE(:,b_half2+1:b-1) = nan;
    CurrentFrE(:,1:b_half2      ) = [];
end

FringMap = CurrentFrE;
end