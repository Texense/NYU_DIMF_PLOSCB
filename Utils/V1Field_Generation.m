% V1Field_Generation.m

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

%% write a bijection between neuron index and corrdinate
% In each HC, we assume that there are 54*54 E neurons and 31*31 I neurons 
% Inut:
%     N_HC: number of HC on each line
%     Nn_Ind: index of neuron, in certain range
%     Nn_Type: 'e' or 'i'
%     varargin: [n_E_HC,n_I_HC], if appicable
% Output:
%     Nn_X, Nn_Y: the spatial intex of the given neuron(s).
function [Nn_X,Nn_Y] = V1Field_Generation(N_HC,Nn_Ind,Nn_Type,varargin)

% first get neuron spatial parameters
if nargin == 3
    [N_E,N_I,X_E,Y_E,X_I,Y_I] = V1Field_Index(N_HC);
else 
    n_E_HC = varargin{1}; n_I_HC = varargin{2};
    [N_E,N_I,X_E,Y_E,X_I,Y_I] = V1Field_Index(N_HC,n_E_HC,n_I_HC);
end
% get neuron location in field for E and I
if strcmpi(Nn_Type,'e') 
    if Nn_Ind<=0 | Nn_Ind>N_E
        disp('*** Invalid neuron index.')
        return        
    else 
        Nn_X = X_E(Nn_Ind);
        Nn_Y = Y_E(Nn_Ind);
    end
elseif strcmpi(Nn_Type,'i') 
    if Nn_Ind<=0 | Nn_Ind>N_I
        disp('*** Invalid neuron index.')
        return        
    else 
        Nn_X = X_I(Nn_Ind);
        Nn_Y = Y_I(Nn_Ind);
    end
else
    disp('*** Invalid neuron type.')
    return
end

end


% In each HC, we assume that there are 54*54 E neurons and 31*31 I neurons 
function [N_E,N_I,X_E,Y_E,X_I,Y_I] = V1Field_Index(N_HC,varargin)
% Number of E and I neurons
if nargin == 1
n_E_HC = 54; n_I_HC = 31; % per side of HC
else 
    n_E_HC = varargin{1}; n_I_HC = varargin{2};
end

N_E = n_E_HC^2 * N_HC^2; % neuron numbers In all
N_I = n_I_HC^2 * N_HC^2;
% spatial intexes of E and I neurons
[X_E,Y_E] = meshgrid(1:n_E_HC*N_HC,1:n_E_HC*N_HC);
[X_I,Y_I] = meshgrid(1:n_I_HC*N_HC,1:n_I_HC*N_HC);
end