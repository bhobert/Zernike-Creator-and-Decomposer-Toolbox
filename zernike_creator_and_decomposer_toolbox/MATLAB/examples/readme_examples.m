%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZERNIKE DECOMPOSITION TOOLBOX EXAMPLES AND README SCRIPT
% -------
%   Created By: Brianna J. Hobert
%   Date Created: 11 Aug. 2026
%   Beta tested: Jim Scire, Sean Peters, Jared Plotkin, Avinash Somwaru
%   Optical Diagnostics Lab at the New York Institute of Technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function provides examples for both quickZern and zernDecomp.
%   type 'help quickZern' and 'help decompZern' for more information.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CITE THIS CODE AS
% -------
% [X] Hobert, Brianna (2026). Zernike Creator and Decomposer Toolbox
%     MATLAB Central File Exchange. Retrieved Month Day, Year. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all; close all; 

%% Section 1: quickZern Example
    % In this example, you will use several solvers in this function
    % to generate defocus and compare each result visually. 
    help quickZern
    
    % Different approaches
      [def,Z20_def,~,~]=quickZern(2,0,"mode","default");
      [gamln,Z20_gamln,~,~]=quickZern(2,0,"mode","gammaln");
      [vp,Z20_vpa,~,~]=quickZern(2,0,"mode","vpa");
      [lstir,Z20_lstir,~,~]=quickZern(2,0,"mode","stirling"); 
      [burn,Z20_burn,~,~]=quickZern(2,0,"mode","burnside"); 
      [mort,Z20_mort,~,~]=quickZern(2,0,"mode","mortici"); 
    
    % Compare in a tiled layout
      crange=[-1,1];
      figure;
          h=tiledlayout(1,6);
              nexttile
                  imagesc(Z20_def,crange)
                  title('Default')
                  axis image; axis off;
              nexttile
                  imagesc(Z20_gamln,crange)
                  title('gammaln')
                  axis image; axis off;
              nexttile
                  imagesc(Z20_vpa,crange)
                  title('vpa')
                  axis image; axis off;
              nexttile
                  imagesc(Z20_lstir,crange)
                  title('Stirling')
                  axis image; axis off;
              nexttile
                  imagesc(Z20_burn,crange)
                  title('Burnside')
                  axis image; axis off;
              nexttile
                  imagesc(Z20_mort,crange)
                  title('Mortici')
                  axis image; axis off;
              colormap("turbo")
              colorbar
      clear all;

%% Section 2: zernDecomp Example
    % In this example, you will generate a wavefront with quickZern
    % and decompose the wavefront, proving that the weights are the same.
    help zernike_decomposition
    
    % Develop artificial wavefront
        N=512;
        [zern1,astig,~,~]=quickZern(2,2,"gridsize",N,"showzerns",1);
        [zern2,~,coma,~]=quickZern(3,1,"gridsize",N,"showzerns",1);
        [zern3,defocus,~,~]=quickZern(2,0,"gridsize",N,"showzerns",1);
        [zern4,spherical,~,~]=quickZern(4,0,"gridsize",N,"showzerns",1);
    
    % Linearly weight components
        wavefront=0.3*astig+0.3*coma+0.4*spherical;
    
    % Decomposition
        [coeffs,~,~]=zernDecomp(wavefront,"showzerns", ... 
                     true,"unwrap",false);