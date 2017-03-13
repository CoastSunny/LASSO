function [hf,hs,hl] = view_surface(figname,faces,verts,cdata)
%VIEW_SURFACE Convenient function to consistently plot surfaces
% function [hf,hs,hl] = view_surface(figname,faces,verts,cdata);
% figname is the name of the figure window
% faces is the triangle listing
%       - if faces is a cell array, verts needs to be a cell array of same length
%         alpha transparency is used to visualize the multiple surfaces.
%         Order in array goes from inner to outer surface.
% verts are the corresponding vertices
%       - verts is either an array or a cell array
%       
% cdata is the colordata to use.  If not given, uses random face color
%       - if cdata is a cell array, use on cell per surface for color coding
% hf is the figure handle used
% hs is the handles to the surfaces
% hl is the handles to the lighting

%<autobegin> ---------------------- 18-Nov-2003 22:02:36 -----------------------
% --------- Automatically Generated Comments Block Using AUTO_COMMENTS ---------
%
% CATEGORY: Visualization
%
% Alphabetical list of external functions (non-Matlab):
%   toolbox\windclf.m
%   toolbox\windfind.m
%
% At Check-in: $Author: Mosher $  $Revision: 13 $  $Date: 11/18/03 8:37p $
%
% Overall BrainStorm authors are:
% ** Dr. John C. Mosher,  Los Alamos National Laboratory
% ** Dr. Sylvain Baillet, CNRS Cognitive Neuroscience & Brain Imaging Laboratory
%
% Copyright (c) 2003 BrainStorm MMIII by the University of Southern California
% Principal Investigator:
% ** Professor Richard M. Leahy, USC Signal & Image Processing Institute
%
% Source code may not be distributed in original or altered form.
% For license and copyright notices, type 'bst_splashscreen' in Matlab,
%   or see BrainStorm website at http://neuroimage.usc.edu,
%   or email Professor Richard M. Leahy at leahy@sipi.usc.edu
%   for further information and permissions.
%<autoend> ------------------------ 18-Nov-2003 22:02:36 -----------------------

% /---Script Authors-------------------------------------\
% |                                                      |
% |  *** John C. Mosher, Ph.D.                           |
% |  Design Technology Group                             |
% |  Los Alamos National Laboratory                      |
% |  Los Alamos, New Mexico, USA                         |
% |  mosher@lanl.gov                                     |
% |                                                      |
% |  *** Sylvain Baillet, Ph.D.                          |
% |  Cognitive Neuroscience & Brain Imaging Laboratory   |
% |  CNRS UPR640 - LENA                                  | 
% |  Hopital de la Salpetriere, Paris, France            |
% |  sylvain.baillet@chups.jussieu.fr                    |
% |                                                      |
% \------------------------------------------------------/

% Script History ----------------------------------------------------------------------------------------
%
% SB  10-Mar-2003 : faces  and verts input arguments can now be cell arrays thereby yielding 3D plots 
%                   using alpha transparency for each surface
% SB  03-Jun-2003 : changed axis management and default view point in 3D
% JCM 19-Aug-2003 : updated comments to explain outputs
% SB  21-Oct-2003 : Basic lightning is 'none'
% --------------------------------------------------------------------------------------------------------

if iscell(faces) % Multiple plots requested on same figure window
    if length(faces) ~= length(verts) % sanity check
        errordlg('Faces and Vertices need to be cell arrays of same length', mfilename);
        return
    end
else
    faces = {faces};
    verts = {verts};
    if nargin > 3
        cdata = {cdata};
    end
    
end

if nargin == 4
    if ~iscell(cdata)
        cdata = {cdata};
    end
end
    

for k=1:length(faces) % For each requested surface
    
    if(size(verts{k},2) > 3), % assume transposed
        verts{k} = verts{k}';  % if the assumption is wrong, will crash below anyway
    end
    
    if nargin == 3 
        cdata{k} = repmat(rand(1,3),size(verts{k},1),1);
    end
    
    h = windfind(figname);
    
    figure(h)
    if isempty(h)
        windclf
        hold on
    end
    
    hs(k) = patch('faces',faces{k},'vertices',verts{k},...
        'facevertexcdata',cdata{k},'facecolor','interp','edgecolor','none');
    
    if length(faces)>1
        set(hs(k),'FaceAlpha',k/10);
    end
    
    view(77,0)
    axis equal, axis vis3d
    axis off
    if nargin == 3
        colormap(bone(256))
    end
        
    material dull
    lighting gouraud
    
    if k ==1 | isempty(findobj(h,'type','light')) % avoid accumulating too many light objects in same window
        hl(1) = camlight(-20,30);
        hl(2) = camlight(20,30);
        hl(3) = camlight(-20,-30);
        for i = 1:length(hl),
            set(hl(i),'color',[.8 1 1]/length(hl)/1.2); % mute the intensity of the lights
        end
    end
    
    if(nargout>0),
        hf = h;  % only if the user has output argument
    end
    
end
