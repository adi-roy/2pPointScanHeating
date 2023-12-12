classdef constants
    properties (Constant)
        wavelength = 1035; % nm
        plotWavelength = []; % nm, plot wavelength may be different from the actual wavelength, if it is invisible.
        tissueModelType = 3; %Data from Johansson, 2010 (Model 1) and Yaroslavsky, 2002 (Model 2)
        mu_a = 0.0478;  % 1/mm. from model 1
        mu_s = 5.8987; %from model 1
        g = 0.9; % anistropy factor
        nTissue = 1.36; % Refractive index of brain
        absGlass = [];
        nGlass = [];
        absWater = [];
        nWater = 1.328;
        scatterFrontCort = @(x) 10.9*((x/500).^(-0.334)); % Human data, x normalized to 500nm, Bevilacqua et al 2000
        scatterTempCort = @(x) 11.6*((x/500).^(-0.601)); % Human data, x normalized to 500nm, Bevilacqua et al 2000

        % Mechanical properties of media
        pTissue = 1.04E-6; % density: kg/mm^3
        pBlood = 1.06E-6; % density: kg/mm^3
        omegaBlood = 8.5E-3; % blood perfusion rate: /s 
        pGlass = 2.23E-6; % density: kg/mm^3
        pWater = 1E-6; % density: kg/mm^3
        pBoneCancellous = 1.178E-6; % density: kg/mm^3
        pBoneCortical = 1.908E-6; % density: kg/mm^3

        % Thermal properties of media
        T_a = 36.7; % deg C: arterial blood temperature
        cTissue = 3.65E6; % specific heat: mJ/kg/C
        kTissue = 0.527; % heat conductivity: mW/m/C
        cBlood = 3.6E6; % specific heat: mJ/kg/C
        qBrain = 9.7E-3; % Brain metabolic heat mW/mm^3
        cGlass = 0.647E6; % specific heat: mJ/kg/C
        kGlass = 0.8; % heat conductivity: mW/m/C
        cWater = 4.184E6; % specific heat: mJ/kg/C
        kWater = 25.72; % heat conductivity: mW/m/C
        cBoneCancellous = 2.274E6; % specific heat: mJ/kg/C
        kBoneCancellous = 0.31; % heat conductivity: mW/m/C
        cBoneCortical = 1.313E6; % specific heat: mJ/kg/C
        kBoneCortical = 0.32; % heat conductivity: mW/m/C
        uStep = 0.0266; % minimum length unit for heat conduction simulation: mm
        T_inf = 25;
        h = 5; %convective heat transfer coeff mW/mm2/C [nat=50e-3, forced=5]

        % Parameters for sample geometry
        dstep=[];
        zrange=[];
        rmax=[];
        d_glass = 0.14; % Thinkness of cover glass: mm
        r_glass = 2; % Radius of cover glass: mm
        d_water = []; % Thickness of immersion water: mm
        d_skull = 0.14; % Thinkness of skull: mm
        d_boneCortical = 0.01; % Thinkness of skull crust: mm

        % Parameters for focusing and scanning geometry
        NA = 0.95;
        f = 9; % Focal length of the objective: mm
        wd = 2; % Work distance of the objective: mm
        w0 = 6; % 1/e^2 radius of the beam at objective back aperture, usually matched to 70% of back aperture diameter for deep imaging. unit: mm
        FOV = 0.1; % Linear FOV of focal scan: mm
        focalDepthTissue = 0; % Focal depth in the tissue: mm 
        surfIllumRadius = []; % The raius of the instersection circle of light cone and brain tissue

        
        
    end
    
end
