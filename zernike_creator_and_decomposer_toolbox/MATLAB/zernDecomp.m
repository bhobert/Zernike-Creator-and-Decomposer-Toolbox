%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZERNDECOMP: ZERNIKE DECOMPOSITION
% -------
%   Created By: Brianna J. Hobert
%   Date Created: 26 Jan. 2026
%   Date Modified: 15 Aug. 2026
%   Beta tested: Jim Scire, Sean Peters, Jared Plotkin, Avinash Somwaru
%   Optical Diagnostics Lab at the New York Institute of Technology
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This routine conducts a Zernike decomposition on a wavefront image
%   with the Moore-Penrose psuedoinverse of the Zernike matrix on the
%   wavefront. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% -------
%   image ............. An image representing the wavefront.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONAL INPUTS
% -------
%   showzerns ......... Show Zernikes within the decomposition matrix.
%   |_ gridsize ....... Size of preview grid for showzerns.
%   |_ crange ......... Display image range for showzerns.
%
%   unwrap ............ Calculate the phase of, and spatially unwrap image. 
%   |_  showphase ..... Show unwrapped phase of image.
%
%   addzerns .......... {Nx1} cell array containing strings of custom
%                       Zernike equations developed by quickZern. 
%                       Refer to the example scripts in this toolbox. 
%
%   pupilsize ......... Adjust pupil size (default 1.0, set from 0 to 1.0).
%   |_ showpupil ...... Show pupil mask on the input wavefront. 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT ARGUMENTS
% -------
%   coeffs ............ Zernike decomposition coefficient array.
%   res ............... Residual map after Zernike decomposition.
%   zer ............... Zernike matrix from decomposition. [1]
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE
% -------
%   % In this example, you will generate a wavefront with quickZern
%   % and decompose the wavefront, proving that the weights are the same.
%
%   % Clear workspace
%         clc; close all; clear all; 
%
%   % Develop artificial wavefront
%         N=512; 
%         [zern1,astig,~,~]=quickZern(2,2,"gridsize",N,"showzerns",1);
%         [zern2,~,coma,~]=quickZern(3,1,"gridsize",N,"showzerns",1);
%         [zern3,defocus,~,~]=quickZern(2,0,"gridsize",N,"showzerns",1);
%         [zern4,spherical,~,~]=quickZern(4,0,"gridsize",N,"showzerns",1);
%
%   % Linearly weight components
%         wavefront=0.3*astig+0.3*coma+0.4*spherical;
%
%   % Decomposition
%         [coeffs,~,~]=zernDecomp(wavefront,"showzerns", ... 
%                      true,"unwrap",false);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD CUSTOM ZERNIKES TO YOUR DECOMPOSITION
% -------
% (1) Use quickZern to generate strings of desired Zernike equations.
% (2) Combine all custom Zernike strings into a {Nx1} cell array.
% (3) Enable "addzerns" in zernDecomp and provide the cell array as an
%     input.
% (4) Run zernDecomp with "showzerns" enabled to ensure that your new
%     additions are adequate. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEE ALSO
% -------
% quickZern, phase_unwrap [2]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REFERENCES
% -------
%  [1] K. Doyle, V. Genberg, G. Michels, "Integrated Optomechanical 
%       Analysis, Second Edition," SPIE Digital Library, 2012,
%       https://doi.org/10.1117/3.974624
%
%  [2] Firman (2026). 2D Weighted Phase Unwrapping 
%      (https://www.mathworks.com/matlabcentral/fileexchange/60345-2d-
%      weighted-phase-unwrapping), 
%      MATLAB Central File Exchange. Retrieved September 12, 2026.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CITE THIS CODE AS
% -------
% [X] Hobert, Brianna (2026). Zernike Creator and Decomposer Toolbox
%     MATLAB Central File Exchange. Retrieved Month Day, Year. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [coeffs, res, zer]=zernDecomp(image,options)
        arguments
            image                                                           % The image that you are decomposing
            options.showzerns (1,1) logical=false                           % An option to show the decomposition Zernikes
            options.gridsize (1,1) double=256                               % Size of preview grid for "showzerns"
            options.crange (1,2) double=[-1,1]                              % Axis limit for "showzerns"
            options.unwrap (1,1) logical=false                              % When on, this calculates the phase and unwraps the image
            options.showphase (1,1) logical=false                           % Show unwrapped phase of image
            options.addzerns (:,1)=0                                        % Optionally add a matrix with extra Zernikes generated by quickZern
            options.pupilsize (1,1) double=1.0                              % Adjust size of pupil (0 to 1)
            options.showpupil (1,1) logical=false                           % Show pupil on input wavefront.
        end
        
        assert(exist('image','var')==1, ['Error: input image does not ' ...
            'exist in the workspace.']);

    % Phase unwrap
        if options.unwrap==1
            phase=phase_unwrap(angle(image));
        else
            phase=image;
        end
    
    % Show phase
        if options.showphase==1
            figure;
                imagesc(phase)
                hold on;
                    set(gcf,'Color','white')
                    colormap(turbo)
                    colorbar
                    axis image; axis off; 
                    title('Unwrapped Phase')
                hold off;
        end

    % Set grid and pupil center up
        [H,W]=size(phase);
        [X,Y]=meshgrid(1:W,1:H);
        cx=(W+1)/2;
        cy=(H+1)/2;
        R=min(H,W)/2;  

    % Normalize coordinates
        X=(X-cx)/R;
        Y=(Y-cy)/R;
        rho=sqrt(X.^2+Y.^2);
        theta=atan2(Y,X);

    % Create a pupil
        mask=rho<=options.pupilsize;
        rho=rho(mask);
        theta=theta(mask);
        if options.showpupil==1
            figure; 
                imagesc(image,options.crange)
                hold on;
                    title('Pupil Preview on Input Wavefront')
                    h=findobj(gca,'Type','Image');
                    h.AlphaData=mask;
                    colormap("turbo");
                    colorbar; 
                    axis off; axis image; 
                hold off;
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 1: INSERT ZERNIKE EQUATIONS HERE MANUALLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zer(1,:)=ones(size(rho));                                                   % Global Piston
zer(2,:)=rho.*cos(theta);                                                   % A-Tilt
zer(3,:)=rho.*sin(theta);                                                   % B-Tilt
zer(4,:)=2*rho.^2-1;                                                        % Focus
zer(5,:)=cos(2*theta).*rho.^2;                                              % Pri Astigmatism A
zer(6,:)=sin(2*theta).*rho.^2;                                              % Pri Astigmatism B
zer(7,:)=cos(3*theta).*rho.^3;                                              % Pri-Trefoil A
zer(8,:)=sin(3*theta).*rho.^3;                                              % Pri-Trefoil B
zer(9,:)=cos(theta).*(3*rho.^3-2*rho);                                      % Pri Coma A
zer(10,:)=sin(theta).*(3*rho.^3-2*rho);                                     % Pri Coma B  
zer(11,:)=cos(4*theta).*rho.^4;                                             % Pri Tetrafoil A
zer(12,:)=sin(4*theta).*rho.^4;                                             % Pri Tetrafoil B
zer(13,:)=6*rho.^4-6*rho.^2+1;                                              % Pri Spherical
zer(14,:)=20*rho.^6-30*rho.^4+12*rho.^2-1;                                  % Sec Spherical
zer(15,:)=70*rho.^8-140*rho.^6+90*rho.^4-20*rho.^2+1;                       % Ter Spherical
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 1: INSERT ZERNIKE EQUATIONS HERE MANUALLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add custom Zernikes, as needed
    if iscell(options.addzerns)
        customzer=options.addzerns;
        customnum=[];
        if size(customzer,2)>1
            customzer=customzer';
        end
        len=size(customzer,1);
        for i=1:len
            customnum=[eval(customzer{i}), customnum];
        end
        customnum=customnum';
        zer=[zer; customnum];
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 2: INSERT CORRESPONDING LABELS FOR ZERNIKE MATRIX HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zernike_names={ ...
    'Piston', ...
    'Tip (A)', ...
    'Tilt (B)', ...
    'Defocus', ...
    'Primary Astig A', ...
    'Primary Astig B', ...
    'Trefoil A', ...
    'Trefoil B', ...
    'Coma A', ...
    'Coma B', ...
    'Tetrafoil A', ...
    'Tetrafoil B', ...
    'Primary Spherical', ...
    'Secondary Spherical', ...
    'Tertiary Spherical'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 2: INSERT CORRESPONDING LABELS TO ZERNIKE MATRIX HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Bug checks for custom zernike names
        zl=size(zernike_names,2);
        zerl=size(zer,1);
        if zl<zerl
            diff=zerl-zl;
            for i=1:diff
                zernike_names{i+zl}=sprintf('Custom Mode %d',i);
            end
        end
    
    % Moore-Penrose pseudoinverse of the Zernike matrix and wavefront
        phasecorrect=phase(mask);
        coeffs=(zer*zer')\(zer*phasecorrect);

    % Plot Zernike decomposition from Zernikes in SECTION 2 ONLY
        figure;
            set(groot, 'DefaultFigureColor', 'w');
            set(groot,'DefaultTextInterpreter','latex');
            set(groot,'DefaultAxesTickLabelInterpreter','latex');
            set(groot,'DefaultLegendInterpreter','latex');
            plot(1:size(zer,1),abs(coeffs),'ko-','LineWidth',2, ...
               'MarkerFaceColor','k');
                hold on;
                    xlabel('Zernike Mode Number','FontSize',18);
                    ylabel('Amplitude','FontSize',18);
                    title('Primary Zernike Decomposition','FontSize',20);
                    axis padded; grid on;
                    y_offset=max(abs(coeffs))*0.01;
                    for k=1:length(coeffs)
                        text(k,abs(coeffs(k))+y_offset, ...
                            zernike_names{k},'Rotation', 45,'FontSize',12);
                    end
                hold off;

    % Residuals
        res=phasecorrect-zer'*coeffs;
        residual_map=zeros(H,W);
        residual_map(mask)=res;

    % Show residuals map
        figure;
            hold on;
                set(gcf,'Color','white')
                imagesc(residual_map);
                axis image; axis off;
                colorbar; colormap(turbo)
                title('Residual Phase After Zernike Fit', ...
                    'FontSize',20);
                h=findobj(gca,'Type','Image');
                h.AlphaData=mask;
            hold off;

    %% Optional: Show residual Zernikes
        if options.showzerns==1
            clear rho; clear theta; clear X; clear Y; clear zer;

                % Establish grid
                    W=options.gridsize; H=options.gridsize;
                    [X,Y]=meshgrid(1:W,1:H);
                    cx=(W+1)/2;
                    cy=(H+1)/2;
                    R=min(H,W)/2; 

                % Normalize coordinates
                    X=(X-cx)/R;
                    Y=(Y-cy)/R;
                    rho=sqrt(X.^2+Y.^2);
                    theta=atan2(Y,X);
                    coords=[rho(:),theta(:)];

                % Sort coordinates
                    [~,idx]=sortrows(coords,'descend');
                    rho=rho(idx);
                    theta=theta(idx);     

                % Create a pupil
                    mask=rho<=options.pupilsize;
                    mask_sorted=zeros(options.gridsize^2,1);
                    mask_sorted(idx)=mask;
                    maskshape=reshape(mask_sorted,options.gridsize, ...
                        options.gridsize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 3: INSERT ZERNIKE EQUATIONS HERE MANUALLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zer(1,:)=ones(size(rho));                                                   % Global Piston
zer(2,:)=rho.*cos(theta);                                                   % A-Tilt
zer(3,:)=rho.*sin(theta);                                                   % B-Tilt
zer(4,:)=2*rho.^2-1;                                                        % Focus
zer(5,:)=cos(2*theta).*rho.^2;                                              % Pri Astigmatism A
zer(6,:)=sin(2*theta).*rho.^2;                                              % Pri Astigmatism B
zer(7,:)=cos(3*theta).*rho.^3;                                              % Pri-Trefoil A
zer(8,:)=sin(3*theta).*rho.^3;                                              % Pri-Trefoil B
zer(9,:)=cos(theta).*(3*rho.^3-2*rho);                                      % Pri Coma A
zer(10,:)=sin(theta).*(3*rho.^3-2*rho);                                     % Pri Coma B  
zer(11,:)=cos(4*theta).*rho.^4;                                             % Pri Tetrafoil A
zer(12,:)=sin(4*theta).*rho.^4;                                             % Pri Tetrafoil B
zer(13,:)=6*rho.^4-6*rho.^2+1;                                              % Pri Spherical
zer(14,:)=20*rho.^6-30*rho.^4+12*rho.^2-1;                                  % Sec Spherical
zer(15,:)=70*rho.^8-140*rho.^6+90*rho.^4-20*rho.^2+1;                       % Ter Spherical
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 3: INSERT ZERNIKE EQUATIONS HERE MANUALLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Add custom Zernikes, as needed
            if iscell(options.addzerns)
                customzer=options.addzerns;
                customnum=[];
                if size(customzer,2)>1
                    customzer=customzer';
                end
                len=size(customzer,1);
                for i=1:len
                    customnum=[eval(customzer{i}), customnum];
                end
                customnum=customnum';
                zer=[zer; customnum];
            end
    
        % Plot all decomposition Zernikes specified in SECTION 3 ONLY
            nplots=size(zer,1);
            figure;
                tx=tiledlayout('flow');
                    title(tx,'Decomposition Zernikes','Interpreter', ...
                        'latex','FontSize',20)
                    set(gcf,'Color','white')
                    for i=1:nplots
                        Z=zer(i,:);
                        Z_sorted=zeros(H*W,1);
                        Z_sorted(idx)=Z;
                        Zshape=reshape(Z_sorted,H,W);
                        t=nexttile;
                            imagesc(Zshape);
                                hold on;
                                    axis image; axis off; axis padded;
                                    h=findobj(gca,'Type','Image');
                                    h.AlphaData=maskshape;
                                    title(t,zernike_names{i},'Interpreter', ...
                                        'latex','FontSize',14)
                                    subtitle(t,sprintf('Z%d',i),'FontSize',12)
                                    clim(options.crange);
                                hold off;
                    end
                    colormap("turbo")
        end
end