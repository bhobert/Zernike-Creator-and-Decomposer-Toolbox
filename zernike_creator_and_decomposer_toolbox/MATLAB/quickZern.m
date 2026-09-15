%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUICKZERN: ZERNIKE GENERATOR
% -------
%   Created By: Brianna J. Hobert
%   Date Created: 2 Feb. 2026
%   Date Modified: 11 Aug. 2026
%   Beta tested: Jim Scire, Sean Peters, Jared Plotkin, Avinash Somwaru
%   Optical Diagnostics Lab at the New York Institute of Technology
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function generates symbolic expressions and wavefronts for 
%   close-form solutions and factorial-free* approximations of 
%   real-valued Zernikes.
% 
%                |m|           
%   ΔZ(rho,theta)       
%                 n
% 
%   Orthogonality rules [1]:
%       n,|m| are Zernike modes and n>=|m|.
%       rho,theta are the radial and angular coordinates.
%
%   *Factorial-free approximations do not need to strictly follow 
%    orthogonality rules.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% -------
%   n ............. Radial mode.
%   |m| ........... Azimuthal frequency. Always represented as abs(m). 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONAL INPUTS
% -------
%   showzerns ......... Show Zernikes within the decomposition matrix.
%   |_ gridsize ....... Size of preview grid for showzerns.
%   |_ crange ......... Display image range for showzerns.
%
%   rms ............... Use root-mean-square Zernike convention.
%
%   mode .............. Choose which solver to use. Documented in EXAMPLE.
%   |                   |_ 'default' uses the closed-form equations.
%   |                   |_ 'gammaln' uses gammaln functions.
%   |                   |_ 'vpa' uses default equations with MATLAB vpa.
%   |                   |_ 'stirling' uses the Stirling approximation* [2].
%   |                   |_ 'burnside' uses the Burnside approximation* [3].
%   |                   |_ 'mortici' uses the Mortici approximation* [4].
%   |
%   |_ normalize ...... Normalize factorial-free Zernikes. Disabled
%                       by default. 
%
%  *All approximations are are continuous, therefore the Parity theorem
%   and constraints regarding n and |m| are relaxed, allowing for the 
%   derivation of "fractional" Zernike polynomials. All approximations 
%   are also calculated in log form. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT ARGUMENTS
% -------
%   zern .............. An array containing both or one strings of the
%                       closed-form Zernike solutions. 
%   wf1 ............... Primary Zernike wavefront.
%   wf2 ............... Secondary wavefront (if azimuthal terms are used).
%   pupil ............. The pupil used as an alpha mask in figures.  
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE
% -------
%   % In this example, you will use several solvers in this function
%   % to generate defocus and compare each result visually. 
%
%   % Clear workspace
%       clc; close all; clear all; 
% 
%   % Different approaches
%       [def,Z20_def,~,~]=quickZern(2,0,"mode","default");
%       [gamln,Z20_gamln,~,~]=quickZern(2,0,"mode","gammaln");
%       [vp,Z20_vpa,~,~]=quickZern(2,0,"mode","vpa");
%       [lstir,Z20_lstir,~,~]=quickZern(2,0,"mode","stirling"); 
%       [burn,Z20_burn,~,~]=quickZern(2,0,"mode","burnside"); 
%       [mort,Z20_mort,~,~]=quickZern(2,0,"mode","mortici"); 
%
%   % Compare in a tiled layout
%       crange=[-1,1];
%       figure;
%           h=tiledlayout(1,6);
%               nexttile
%                   imagesc(Z20_def,crange)
%                   title('Default')
%                   axis image; axis off;
%               nexttile
%                   imagesc(Z20_gamln,crange)
%                   title('gammaln')
%                   axis image; axis off;
%               nexttile
%                   imagesc(Z20_vpa,crange)
%                   title('vpa')
%                   axis image; axis off;
%               nexttile
%                   imagesc(Z20_lstir,crange)
%                   title('Stirling')
%                   axis image; axis off;
%               nexttile
%                   imagesc(Z20_burn,crange)
%                   title('Burnside')
%                   axis image; axis off;
%               nexttile
%                   imagesc(Z20_mort,crange)
%                   title('Mortici')
%                   axis image; axis off;
%               colormap("turbo")
%               colorbar
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEE ALSO
% -------
% zernDecomp, vpa (MATLAB Symbolic Math Toolbox)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REFERENCES
% -------
%  [1] K. Doyle, V. Genberg, G. Michels, "Integrated Optomechanical 
%       Analysis, Second Edition," SPIE Digital Library, 2012,
%       https://doi.org/10.1117/3.974624
%
%  [2] Tweddle, I. (2003). Introduction. In: James Stirling's
%      Methodus Differentialis. Sources and Studies in the History of 
%      Mathematics and Physical Sciences. 
%      Springer, London. https://doi.org/10.1007/978-1-4471-0021-8_1
%
%  [3] W. Burnside, "A rapidly convergent series for logN!," 
%      Messenger Math., vol. 46, pp. 157–159, 1917.
%
%  [4] C. Mortici, "An ultimate extremely accurate formula 
%      for approximation of the factorial function," 
%      Arch. Math. (Basel), vol. 93, no. 1, pp. 37–45, 2009, 
%      doi: 10.1007/s00013-009-0008-5.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CITE THIS CODE AS
% -------
% [X] Hobert, Brianna (2026). Zernike Creator and Decomposer Toolbox
%     MATLAB Central File Exchange. Retrieved Month Day, Year. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [zern,wf1,wf2,pupil]=quickZern(n,m,options)
    arguments
        n                                                                   % Rotational mode
        m                                                                   % Azimuthal frequency
        options.showzerns (1,1) logical=false                               % An option to show the decomposition Zernikes
        options.gridsize (1,1) double=256                                   % Size of preview grid for "showzerns"
        options.crange (1,2) double=[-1,1]                                  % Axis limit for "showzerns"
        options.rms (1,1) logical=false                                     % Enable rms Zernikes
        options.mode (1,1) string='default'                                 % Choose between approximation modes
        options.normalize (1,1) logical=false                               % Clamp factorial-free Zernikes to [-1,1] range
    end

    % Formal warnings
        if m<0 || n<0
            warning('Negative Zernike modes detected. m=|m|.')
            m=abs(m);
            n=abs(n);
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RMN CALCULATION SWITCH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R^m_n(rho)=(A*B*term)/(D*E*F)
% A=(-1)^s
% B=(n-s)!
% term=rho^(power)
%   power=n-2s
% D=s!
% E=(((n+m)/2)-s)!
% F=(((n-m)/2)-s)!
switch options.mode
    case 'gammaln'
        % Gammaln mode
            fprintf('<strong>Gamma-ln Mode </strong>\n')
            assert(mod(n,1)==0, 'n must be an integer-valued', ... 
                'number for this mode.');
            assert(mod(m,1)==0, 'm must be an integer-valued ', ...
                'number for this mode.');
            assert(n>=m,'The condition n>=m must be met for ', ...
                'this mode.')
            if m>n||mod(n-m,2)~=0
                error('n-m must be even for this mode.');
            end
            k=(n-m)/2;
            coeffs=zeros(1,k+1);
            powers=zeros(1,k+1);
            for s=0:k
                A=(-1)^s;
                logB=gammaln(n-s+1);
                logD=gammaln(s+1);
                logE=gammaln((n+m)/2-s+1);
                logF=gammaln((n-m)/2-s+1);
                coeff=A*exp(logB-logD-logE-logF);
                coeffs(s+1)=coeff;
                powers(s+1)=n-2*s;
            end
            Rmn='';
            for i=1:length(coeffs)
                term=[num2str(coeffs(i)),'.*rho.^', ...
                    num2str(powers(i))];
                if i==1
                    Rmn=term;
                else
                    Rmn=[Rmn,'+',term];
                end
            end

    case 'default'
        % Default solver
            fprintf('<strong>Default Mode </strong>\n')
            assert(mod(n,1)==0, 'n must be an integer-valued', ... 
                'number for this mode.');
            assert(mod(m,1)==0, 'm must be an integer-valued ', ...
                'number for this mode.');
            assert(n>=m,'The condition n>=m must be met for ', ...
                'this mode.')
            if m>n||mod(n-m,2)~=0
                error('n-m must be even for this mode.');
            end
            k=(n-m)/2;
            Rmn='';
            for s=0:k
                A=(-1)^s;
                B=factorial(n-s);
                D=factorial(s);
                E=factorial(((n+m)/2)-s);
                F=factorial(((n-m)/2)-s);
                coeff=(A*B)/(D*E*F);
                power=n-2*s;
                term=[num2str(coeff),'.*rho.^(',num2str(power),')'];
                if s==0
                    Rmn=term;
                else
                    Rmn=[Rmn,'+',term];
                end
            end

    case 'vpa'
       % Variable Precision Arithmetic mode
            fprintf('<strong>VPA Mode </strong>\n')
            if m>n||mod(n-m,2)~=0
               error(['Fractional Zernike detected. |m|<=n and' ...
                   ' n-m must classically be even for this mode.']);
            end
            k=(n-m)/2;
            Rmn='';
            for s=0:k
                A=(-1)^s;
                B=factorial(vpa(n-s));
                D=factorial(vpa(s));
                E=factorial(vpa(((n+m)/2)-s));
                F=factorial(vpa(((n-m)/2)-s));
                coeff=double((A*B)/(D*E*F));
                power=n-2*s;
                term=[num2str(coeff),'.*rho.^(',num2str(power),')'];
                if s==0
                    Rmn=term;
                else
                    Rmn=[Rmn,'+',term];
                end
            end

    case 'stirling'
        % Logarithmic Stirling mode
            fprintf('<strong>Log Stirling Mode</strong>\n')
            fprintf(2,'Advisory: This mode is an approximation.\n')
            if m>n||mod(n-m,2)~=0
               warning(['Fractional Zernike detected. |m|<=n and' ...
                    ' n-m must classically be even for this mode.']);
            end
            k=(n-m)/2;
            Rmn='';
            for s=0:k
                A=(-1)^s;
                B=log(stirfac(n-s));
                D=log(stirfac(s));
                E=log(stirfac(((n+m)/2)-s));
                F=log(stirfac(((n-m)/2)-s));
                coeff=A*exp(B-D-E-F);
                power=n-2*s;
                term=[num2str(coeff),'.*rho.^(',num2str(power),')'];
                if s==0
                    Rmn=term;
                else
                    Rmn=[Rmn,'+',term];
                end
            end

    case 'burnside'
        % Burnside approximation mode
            fprintf('<strong>Log Burnside Mode</strong>\n')
            fprintf(2,'Advisory: This mode is an approximation.\n')
            if m>n||mod(n-m,2)~=0
               warning(['Fractional Zernike detected. |m|<=n and' ...
                    ' n-m must classically be even for this mode.']);
            end
            k=(n-m)/2;
            Rmn='';
            for s=0:k
                A=(-1)^s;
                B=log(burnfac(n-s));
                D=log(burnfac(s));
                E=log(burnfac(((n+m)/2)-s));
                F=log(burnfac(((n-m)/2)-s));
                coeff=A*exp(B-D-E-F);
                power=n-2*s;
                term=[num2str(coeff),'.*rho.^(',num2str(power),')'];
                if s==0
                    Rmn=term;
                else
                    Rmn=[Rmn,'+',term];
                end
            end

    case 'mortici'
        % Mortici approximation mode
            fprintf('<strong>Log Mortici Mode</strong>\n')
            fprintf(2,'Advisory: This mode is an approximation.\n')
            if m>n||mod(n-m,2)~=0
               warning(['Fractional Zernike detected. |m|<=n and' ...
                    ' n-m must classically be even for this mode.']);
            end
            k=(n-m)/2;
            Rmn='';
            for s=0:k
                A=(-1)^s;
                B=mortfac(n-s);
                D=mortfac(s);
                E=mortfac(((n+m)/2)-s);
                F=mortfac(((n-m)/2)-s);
                coeff=A*exp(B-D-E-F);
                power=n-2*s;
                term=[num2str(coeff),'.*rho.^(',num2str(power),')'];
                if s==0
                    Rmn=term;
                else
                    Rmn=[Rmn,'+',term];
                end
            end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RMN CALCULATION SWITCH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Finish developing strings based on azimuthal and rotational terms
        Rmn=strrep(Rmn,'+-','-');
        Rmn=strrep(Rmn,'.*rho.^(0)','');
        Rmn=strrep(Rmn,'.^(1)','');
        Rmn=regexprep(Rmn,'.^\+','');
        if m==0
            zern=Rmn;
            zern=strrep(zern,'(1.*','(');
            zern=strrep(zern,'(1*','(');
            if options.rms==true
                normroot=sprintf('sqrt(%d)',n+1);
                zern=[normroot,'*(',zern,')'];
            end
            disp(['<strong>Z^',num2str(m),'_',num2str(n), ...
                '(rho,theta)</strong>=',zern]);
        else
            Zmn1=['(',Rmn,').*cos(',num2str(m),'*theta)'];
            Zmn2=['(',Rmn,').*sin(',num2str(m),'*theta)'];
            Zmn1=strrep(Zmn1,'(1.*','(');
            Zmn1=strrep(Zmn1,'(1*','(');
            Zmn2=strrep(Zmn2,'(1.*','(');
            Zmn2=strrep(Zmn2,'(1*','(');
            if options.rms==true
                numroot=2*(n+1);
                normroot=sprintf('sqrt(%d)',numroot);
                Zmn1=[normroot,'*(',Zmn1,')'];
                Zmn2=[normroot,'*(',Zmn2,')'];
            end
            zern=[Zmn1;Zmn2];
            disp(['<strong>Z^',num2str(m),'_',num2str(n), ...
                '(rho,theta)</strong>=',zern(1,:)]);
            disp(['<strong>Z^',num2str(m),'_',num2str(n), ...
                '(rho,theta)</strong>=',zern(2,:)]);
        end

    % Establish grid
        [X,Y]=meshgrid(1:options.gridsize,1:options.gridsize);
        cx=(options.gridsize+1)/2;
        cy=(options.gridsize+1)/2;
        R=min(options.gridsize,options.gridsize)/2;

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
        mask=rho<=1.0;
        mask_sorted=zeros(options.gridsize^2,1);
        mask_sorted(idx)=mask;
        maskshape=reshape(mask_sorted,options.gridsize,options.gridsize);
        pupil=maskshape;
        if options.rms==true
            plttitle=sprintf('Real-valued RMS $Z^{%d}_{%d}$',m,n);
        else
            plttitle=sprintf('Real-valued $Z^{%d}_{%d}$',m,n);
        end

    % Evaluate string expressions and display generated Zernikes
        if m==0
            % Rotationally symmetric approach
                Z=eval(zern);
                Z_sorted=zeros(options.gridsize^2,1);
                Z_sorted(idx)=Z;
                Zshape=reshape(Z_sorted,options.gridsize,options.gridsize);
                if strcmp(options.mode,'stirling')==1 | ...
                        strcmp(options.mode,'burnside')==1 | ...
                        strcmp(options.mode,'mortici')==1 && ...
                        options.normalize==1
                    Zshape=Zshape/max(Zshape(:));
                end
                if options.showzerns==1
                    figure;
                        set(gcf,'defaultTextInterpreter','latex', ...
                            'Color','white');
                        imagesc(Zshape);
                            hold on;
                                title(plttitle,'FontSize',20)
                                colormap('turbo'); axis image; axis off;
                                h=findobj(gca,'Type','Image');
                                h.AlphaData=maskshape;
                                colorbar('Position', [0.82 0.1 0.035 0.8]);
                                clim(options.crange);
                            hold off;
                end
                wf1=Zshape; wf2=wf1;
        else
            % Azimuthal approach
                Z1=eval(zern(1,:));
                    Z1_sorted=zeros(options.gridsize^2,1);
                    Z1_sorted(idx)=Z1;
                    Z1shape=reshape(Z1_sorted,options.gridsize, ...
                        options.gridsize);
                Z2=eval(zern(2,:));
                    Z2_sorted=zeros(options.gridsize^2,1);
                    Z2_sorted(idx)=Z2;
                    Z2shape=reshape(Z2_sorted,options.gridsize, ...
                        options.gridsize);
                if strcmp(options.mode,'stirling')==1 | ...
                        strcmp(options.mode,'burnside')==1 | ...
                        strcmp(options.mode,'mortici')==1 && ...
                        options.normalize==1
                    Z1shape=Z1shape/max(Z1(:));
                    Z2shape=Z2shape/max(Z2(:));
                end
                if options.showzerns==1
                    fig=figure;
                        sgtitle(plttitle,'FontSize',30)
                        set(gcf,'defaultTextInterpreter','latex', ...
                            'Color','white');
                        subplot(1,2,1)
                            imagesc(Z1shape);
                                hold on;
                                    title('Cosine variant','FontSize',20)
                                    axis square; axis off; 
                                    h=findobj(gca,'Type','Image');
                                    h.AlphaData=maskshape;
                                    clim(options.crange);
                                hold off;
                        subplot(1,2,2)
                            imagesc(Z2shape);
                                hold on;
                                    title('Sine variant','FontSize',20)
                                    axis square; axis off;
                                    h=findobj(gca,'Type','Image');
                                    h.AlphaData=maskshape;
                                    clim(options.crange);
                                hold off;
                        colorbar('Position', [0.93 0.2 0.025 0.6]);
                        colormap(fig,'turbo');
                end
            wf1=Z1shape; wf2=Z2shape; 
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOCAL FUNCTIONS FOR APPROXIMATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stirling approximation
    function F=stirfac(x)
        if x==0
            F=1;
        else
            F=sqrt(2*pi*x)*(x/exp(1))^x;
        end
    end

% Burnside approximation
    function F=burnfac(x)
        F=sqrt(2*pi)*((x+0.5)/exp(1))^(x+0.5);
    end

% Mortici approximation
    function F=mortfac(x)
        F=log(sqrt(2*pi)*((x^2+x+1/6)/exp(2))^(x/2+1/4));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOCAL FUNCTIONS FOR APPROXIMATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
