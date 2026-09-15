%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HOW TO ADD CUSTOM ZERNIKES TO ZERNDECOMP
% -------
%   Created By: Brianna J. Hobert
%   Date Created: 13 Aug. 2026
%   Beta tested: Jim Scire, Sean Peters, Jared Plotkin, Avinash Somwaru
%   Optical Diagnostics Lab at the New York Institute of Technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function shows how to add custom Zernikes from quickZern
%   to zernDecomp through code only. 
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

%% Step 1: Generate Custom Zernikes
    [astig2,wf5,~,~]=quickZern(4,2,'showzerns',1);
    [spherical5,wf4,~,~]=quickZern(10,0,'showzerns',1);
    [frac,~,~,~]=quickZern(2,1,'showzerns',1,'mode','burnside');

%% Step 2: Add all custom equations into cells
    customzerns={astig2(1,:), astig2(2,:), spherical5, frac(1,:), ...
        frac(2,:)};

%% Step 3: Make Dummy Wavefront
    [~,wf1,~,~]=quickZern(2,0);
    [~,wf6,~,~]=quickZern(5,1);
    [~,frac1,frac2,~]=quickZern(2,1,"mode","burnside");
    wavefront=0.2*wf1+0.5*wf4+0.4*wf5+0.1*frac1+0.6*frac2+0.01*wf6;

%% Step 4: Enable "addzerns"
    zernDecomp(wavefront,"addzerns",customzerns,"showzerns",1,"unwrap",0)