%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT 3D WAVEFRONT
% -------
%   Created By: Brianna J. Hobert
%   Date Created: 13 Aug. 2026
%   Beta tested: Jim Scire, Sean Peters, Jared Plotkin, Avinash Somwaru
%   Optical Diagnostics Lab at the New York Institute of Technology
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function takes Zernike wavefronts and plots them in 3D. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% -------
%   wavefront ......... A 2D matrix containing the Zernike wavefront. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONAL INPUTS
% -------
%   pupil ............. Add in a mask. Must contain only 1's and 0's.
%   title ............. Add figure title (string).
%   xlabel ............ Add an x-axis label (string).
%   ylabel ............ Add a y-axis label (string).
%   zlabel ............ Add a z-axis label (string).
%   xlim .............. Add limits to the x axis.
%   ylim .............. Add limits to the y axis.
%   zlim .............. Add limits to the z axis.
%   fontsize .......... Adjust figure font size (1x1 double).
%   cmap .............. Set colourmap (string).
%   cbar .............. Toggle colourbar (logical).
%   background ........ Set figure background (string).
%   axis .............. Set figure axis.
%   grid .............. Set figure grid (logical).
%   lighting .......... Set figure lighting (logical).
%   edgecolor ......... Set edge color.
%   view .............. Change view angles of camera (1x2).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT ARGUMENTS
% -------
%   fig ............... The figure object generated in this function.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE
% -------
%   % In this example, you will use quickZern to make a wavefront and 
%   % plotzern to show it as a 3D surface.
%   
%   % Clear workspace
%       clc; clear all; close all; 
% 
%   % Make Zernike and plot it
%       [zern,wf1,wf2,pupil]=quickZern(4,4,'gridsize',1024);
%       fig=plotzern(wf1,'pupil',pupil, ...
%           'xlabel','X Axis', ...
%           'ylabel','Y Axis', ...
%           'zlabel','Z Axis', ...
%           'axis',true, ...
%           'title','Zernike Quadfoil', ...
%           'lighting','gouraud', ...
%           'cmap','jet', ...
%           'view',[-40,50], ...
%           'fontsize',24);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEE ALSO
% -------
% quickZern
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CITE THIS CODE AS
% -------
% [X] Hobert, Brianna (2026). Zernike Creator and Decomposer Toolbox
%     MATLAB Central File Exchange. Retrieved Month Day, Year. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fig=plotzern(wavefront,options)
    arguments
        wavefront                                                           % Input image
        options.title (1,:)=[];                                             % Add a title
        options.pupil (:,:) double=[];                                      % Add a pupil alpha mask in
        options.xlabel (1,:)=' ';                                           % Add an x-axis label
        options.ylabel (1,:)=' ';                                           % Add a y-axis label
        options.zlabel (1,:)=' ';                                           % Add a z-axis label
        options.zlim (1,:)=[-2,2];                                          % Add a z limit
        options.xlim (1,:)='auto';                                          % Add an x limit
        options.ylim (1,:)='auto';                                          % Add a y limit
        options.fontsize (1,1)=20;                                          % Adjust figure font size
        options.cmap (1,:)='turbo'                                          % Colormap
        options.cbar (1,1) logical=true                                     % Add a colourbar
        options.background (1,:)='white'                                    % Change figure background color
        options.axis (1,:)='off'                                            % Add or remove axis
        options.grid (1,:)='off'                                            % Add or remove grid
        options.lighting (1,:)='none'                                       % Add or remove lighting
        options.edgecolor (1,:)='none'                                      % Add or remove wireframe
        options.view (1,2)=[336.50 41.00]                                   % Change view angle of camera (1x2)
    end
    
    if ~isempty(options.pupil)
        assert(size(options.pupil,1)==size(wavefront,1) && ...
            size(options.pupil,2)==size(wavefront,2), ...
            ['Error: Pupil image and wavefront' ...
            ' image must be the same size.'])
        options.pupil(options.pupil==0)=NaN;
        wavefront=options.pupil.*wavefront;
    end
    
    % Make a surf
        fig=figure;
            mesh(wavefront,'EdgeColor',options.edgecolor, ...
                'FaceColor','interp','FaceLighting',options.lighting)
            set(gcf,'Color',options.background)
            hold on;
                view(options.view(1),options.view(2)); 
                if strcmp(options.lighting,'none')~=1
                    material('shiny')
                    % Key light
                        light_key=light('Position',[-10,0,20], ...
                            'Style','local');
                        set(light_key,'Color',[1.0, 1.0, 0.2]); 

                    % Fill light
                        light_fill=light('Position',[15,-10,10], ...
                            'Style','local');
                        set(light_fill,'Color',[0.2,0.2,0.5]); 

                    % Rim light
                        light_rim=light('Position',[0, 20, 5], ...
                            'Style','local');
                        set(light_rim,'Color',[0.4, 0.4, 0.3]); 
                end
                if ~isempty(options.title)
                    title(options.title,'Interpreter','latex', ...
                        'FontSize',options.fontsize)
                else
                    title('Wavefront','Interpreter','latex','FontSize', ...
                        options.fontsize)
                end
                xlim(options.xlim)
                ylim(options.ylim)
                zlim(options.zlim)
                axis(options.axis)
                grid(options.grid)
                if ~isempty(options.xlabel)
                   xlabel(options.xlabel,'FontSize',options.fontsize, ...
                       'Interpreter','latex')
                end
                if ~isempty(options.ylabel)
                    ylabel(options.ylabel,'FontSize',options.fontsize, ...
                        'Interpreter','latex')
                end
                if ~isempty(options.zlabel)
                    zlabel(options.zlabel,'FontSize',options.fontsize, ...
                        'Interpreter','latex')
                end
                colormap(options.cmap)
                if options.cbar==1
                    colorbar;
                end
            hold off;
end