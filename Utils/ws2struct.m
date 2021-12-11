% ws2struct.m

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

%This Script saves all the variables from the current workspace into a 
%single structure array. 
%Author: Andres Gonzalez
%Year:2012
%Version:0.2
%#########################################################################
%Summary:
%This function allows to save all the variables from the 'caller' workspace into a struct array
%Description:
%Sometimes you need to save the variables from your base workspace, but using "save" function will have them all stored individually so if you reload them into a new workspace it could be a mess, and some variables could be overwritten.
%With this function, you can save all of them into a struct array, and so they'll be nicely packaged and ready to be saved to a .mat file that, when reloaded, will be easy to identify. 
%Tags:
%save, struct, base, workspace, structure array
%Example:
% a='LALALA'
% b=[1:12:258]
% c={'cell1', 'cell2', 'cell3'}
% d=768
% e=true(3)
% theworkspace=ws2struct();
% theworkspace = 
% 
%     a: 'LALALA'
%     b: [1x22 double]
%     c: {'cell1'  'cell2'  'cell3'}
%     d: 768
%     e: [3x3 logical]
function WStruct=ws2struct()
WSVARS = evalin('caller', 'who');
for wscon=1:size(WSVARS,1)
    thisvar=evalin('caller', WSVARS{wscon});
    THEWORKSPACE.(WSVARS{wscon})=thisvar;
end
WStruct=THEWORKSPACE;