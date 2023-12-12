%% MC-FDM package for predicting heating in 2p point scan imaging

% This code uses MCML based Monte Carlo simulations to 
% calculate fluence rate matrices for various gaussian beam
% focusing parameters. It then couples to an explicit finite difference
% method solver to calculate thermal profiles. 

% The highlight of this code is that it simulates heating resulting from 
% flow of the objective immersion water layer. 

% This code currently implements parallel processing with vectorization and
% runs on the CPU. GPU acceleration is work under progress.

% The heat transfer package implements Stujenske et al.'s explicit FDM
% discretization. A much faster implicit discretization scheme will be
% released soon. GPU acceleration + implicit discretization will make the
% code orders of magnitude faster.

% beamfocus.m is called by the MC module to determine beam focusing
% conditions. beamfocus.m uses a Cholesky discretization based
% algorithm to determine photon coordinates at the surface. This algorithm
% also makes it possible to simulate bivariate gaussian beams instead of
% just univariate gaussians (in xy).

% Several algorithms in this package have been adapted from the works of
% Stujenske et al. 2015, Podgorski and Ranganathan 2016, and 
% Wang et al. 2020. 
% Credits to their papers:
% 1. Stujenske JM, Spellman T, Gordon JA. Modeling the Spatiotemporal Dynamics of Light and Heat Propagation for In Vivo Optogenetics. Cell Rep. 2015 Jul 21;12(3):525-34. doi: 10.1016/j.celrep.2015.06.036. Epub 2015 Jul 9. PMID: 26166563; PMCID: PMC4512881.
% 2. Podgorski K, Ranganathan G. Brain heating induced by near-infrared
% lasers during multiphoton microscopy. Journal of Neurophysiology 2016. https://doi.org/10.1152/jn.00275.2016
% 116:3, 1012-1023. 
% 3. Tianyu Wang, Chunyan Wu, Dimitre G Ouzounov, Wenchao Gu, Fei Xia, Minsu Kim, Xusan Yang, Melissa R Warden, Chris Xu (2020) Quantitative analysis of 1300-nm three-photon calcium imaging in the mouse brain eLife 9:e53205

    

% By Aditya Roy, December 2022 %





%%
close all;
clear all;
clc;


%% Initialize model paraneters

% This code uses constants.m to obtain various simulation constants

radial_domain_dimension = 4; % depth is always set to 4
mc_stepsize = 0.01;
heat_stepsize = 0.0285; 
start_imaging_depth = 0; % in terms of attenuation lengths
end_imaging_depth = 1; % in terms of attenuation lengths
n_imaging_depths = end_imaging_depth - start_imaging_depth + 1;
imaging_depths = linspace(start_imaging_depth,end_imaging_depth,n_imaging_depths);
mu_t = constants.mu_a + constants.mu_s;
heat_starting_depth = -1*(constants.wd - imaging_depths/mu_t);

%% User inputs

heat_incubation_time = 60; % in sec
heat_simulation_time = 60; % in sec

flag_duty = 1; % Set 1 for implementing duty cycles
duty = 50; % 100% means no duty cycle
cycle_time = 0.1; 

power = [50 100];
Tcold = 25; % water temperature, both for flowing and stagnant

flag_cooling = 2; % Use flag = 1 for conventional technique, anything else for convection
hconv = 5; %mW/mm^2-K
T_inf = Tcold;


%% MC module

[weight_mat, space] = mc(mc_stepsize, radial_domain_dimension, n_imaging_depths, imaging_depths, mu_t);

%% HT module

[source, heat_x, heat_z] = fluence_rate(mc_stepsize, radial_domain_dimension, heat_stepsize,...
    heat_starting_depth, weight_mat, power, n_imaging_depths);

%
[T_space] = heat(heat_x, heat_z, heat_incubation_time, heat_simulation_time, flag_duty, duty, cycle_time, Tcold, power, ...
    source, heat_stepsize, n_imaging_depths, imaging_depths/mu_t, flag_cooling, hconv, T_inf);

%% Save
filename = 'run.mat';
save(filename);

% not saving plots by default. plot_contour.m can be used as a function to plot MC and heat maps 