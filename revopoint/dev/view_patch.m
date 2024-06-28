% Copyright 2010 Anand A. Joshi, David W. Shattuck and Richard M. Leahy 
% This file is part SVREG.
% 
% SVREG is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% SVREG is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with BSE.  If not, see <http://www.gnu.org/licenses/>.

function h=view_patch(FV1,th)
%function view_patch(FV)
%FV: a tessellation to view
%
%copywrite Dimitrios Pantazis, PhD student, USC

%h=figure;
%camlight
%lighting gouraud
axis equal
axis off   
axis vis3d
FV.vertices=FV1.vertices;FV.faces=FV1.faces;
nVertices = size(FV.vertices,1);
clr=hsv(1000);
clr(1:round(1000*2/6),1)=241;clr(1:round(1000*2/6),2)=241;clr(1:round(1000*2/6),3)=241;
clr(1+round(1000*2/6):round(1000*2.5/6),1)=113;clr(1+round(1000*2/6):round(1000*2.5/6),2)=188;clr(1+round(1000*2/6):round(1000*2.5/6),3)=226;
clr(1+round(1000*2.5/6):round(1000*3/6),1)=36;clr(1+round(1000*2.5/6):round(1000*3/6),2)=147;clr(1+round(1000*2.5/6):round(1000*3/6),3)=42;
clr(1+round(1000*3/6):round(1000*3.5/6),1)=240;clr(1+round(1000*3/6):round(1000*3.5/6),2)=200;clr(1+round(1000*3/6):round(1000*3.5/6),3)=40;
clr(1+round(1000*3.5/6):round(1000),1)=240;clr(1+round(1000*3.5/6):round(1000),2)=54;clr(1+round(1000*3.5/6):round(1000),3)=39;
clr=clr/256;
hpatch = patch(FV,'FaceColor','interp','EdgeColor','none','FaceVertexCData',clr(max(min(round(1000*th/6),1000),1),:),'faceAlpha',1);%,'BackFaceLighting','unlit'); %plot surface        
lighting gouraud
caxis([0,6]);colormap hsv;

