%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZERNIKE PYRAMID
% -------
%   Created By: Brianna J. Hobert
%   Date Created: 12 Aug. 2026
%   Beta tested: Jim Scire, Sean Peters, Jared Plotkin, Avinash Somwaru
%   Optical Diagnostics Lab at the New York Institute of Technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function generates the Zernike "pyramid" with extended modes
%   using the specifed approximation in line 26. 
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
wf={}; wf2={}; wf3={}; fcell={};
N=128;                                                                      % Grid size.
crange=[-1,1];                                                              % Color axis limit range.
levels=7;                                                                   % Pyramid levels to render.
approximator='burnside';                                                    % Approximator for intermittent Zernikes.
rmsmode=false;                                                              % Enables rms Zernikes.
normalized=false;                                                           % Normalize factorial-free Zernikes.

%% Rotationally Symmetric Terms
    % Higher terms
        for i=1:levels+1
            if mod(i,2)==0
                if i==2
                    [~,wf{i},~,pupil]=quickZern(i,0,"mode","default", ...
                        "gridsize",N,"rms",rmsmode,"normalize",normalized);
                else
                    [~,wf{i},~,~]=quickZern(i,0,"mode","default", ...
                        "gridsize",N,"rms",rmsmode,"normalize",normalized);
                end
            else
                [~,wf{i},~,~]=quickZern(i,0,"mode",approximator, ...
                    "gridsize",N,"rms",rmsmode,"normalize",normalized);
            end
        end

    % Add piston
        wf{:,end}=[];
        wf=[ones(N,N), wf];
        wf=wf(~cellfun(@isempty,wf));

%% Azimuthal Terms
    for m=1:levels
        for n=1:levels
            if mod(n-m,2)==0 && m<=n
                [~,wf1{n,m},wf2{n,m},~]=quickZern(n,m,"mode","default", ...
                    "gridsize",N,"rms",rmsmode,"normalize",normalized);
            elseif m<=n
                [~,wf1{n,m},wf2{n,m},~]=quickZern(n,m,"mode", ...
                    approximator,"gridsize",N,"rms",rmsmode, ...
                    "normalize",normalized);
            end
        end
    end
    for k=1:levels
        fcell{k}=0;
    end
    wf1=[fcell; wf1];
    wf2=[fcell; wf2];

%% Assemble Cells
    wf1=flip(wf1,2);
    wfs=[wf1, wf', wf2];
    wfs=rot90(wfs);

%% Make a pyramid of Zernikes
    W=levels+1; H=2*levels+1; tiles=W*H;
    figure;
        t=tiledlayout(W,H,'TileSpacing','compact','TileIndexing', ...
            'rowmajor');
        title(t,['MATLAB Zernike Creator and Decomposer Toolbox: ' ...
            'Zernike Mode Pyramid'],'Interpreter','latex','FontSize',20)
        nidx=0;
        for i=1:W
            for j=1:H
              midx=j-(levels+1);
                if isempty(wfs{j,i}) | wfs{j,i}==0
                    nexttile;
                        title(' ');
                        axis off; axis image;
                else
                    nexttile;
                        imagesc(wfs{j,i},crange)
                        Znum=sprintf('$Z_{%d}^{%d}$',nidx,midx);
                        text(N/2,N+0.17*N,Znum,'HorizontalAlignment', ...
                            'center','FontSize',12,'Interpreter','latex')
                        h=findobj(gca,'Type','Image');
                        axis off; axis image;
                        h.AlphaData=pupil;
                        colormap turbo
                end
            end
            nidx=nidx+1;
        end