% plot_dir.m

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

function [h1, h2] = plot_dir (vX, vY)
%function [h1, h2] = plot_dir (vX, vY)
%Plotting x-y variables with direction indicating vector to the next element.
%Example
%   vX = linspace(0,2*pi, 10)';
%   vY = sin (vX);
%   plot_dir(vX, vY);
%% NOTE!! vX, vY has to be column vectors!
rMag = 0.5;
% Length of vector
lenTime = length(vX);
% Indices of tails of arrows
vSelect0 = 1:(lenTime-1);
% Indices of tails of arrows
vSelect1 = vSelect0 + 1;
% X coordinates of tails of arrows
vXQ0 = vX(vSelect0, 1);
% Y coordinates of tails of arrows
vYQ0 = vY(vSelect0, 1);
% X coordinates of heads of arrows
vXQ1 = vX(vSelect1, 1);
% Y coordinates of heads of arrows
vYQ1 = vY(vSelect1, 1);
% vector difference between heads & tails
vPx = (vXQ1 - vXQ0) * rMag;
vPy = (vYQ1 - vYQ0) * rMag;
% make plot 
h1 = plot (vX, vY, '.k-'); hold on;
% add arrows 
h2 = quiver (vXQ0,vYQ0, vPx, vPy, 0, 'r'); grid on; hold off
%axis equal
