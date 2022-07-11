% % % % % % % % % % % % % % % % % % % % %
% Run Doorway Camera
% 
% Function call
%   "doorway_camera(labData,scnNum,samps,facetHeight,fov,res,showPlots,bgRemoval)"
%
% Input Parameters (in order)
%   labData: 1-Laboratory Data or 0-Synthetic Data
%   scnNum: Scene Number (see below)
%   samps: Aquisition time (synthetic only, see below)
%   facetHeight: Height of facets in forward model (m)
%   fov: Size of FOV in forward model (both h and w)
%   res: resolution of measurment (square). Will downsample to match
%   lambda: tuningparam
%   showPlots: Active-Plot 1-on or 0-off
%   bgRemoval: Background Removal 1-on or 0-off
%
% Scene numbers
%   Laboratory Data
%       0,1,2,3 - Single red target (in order of appearance in paper)
%       4,5 - Two targets (in order of appearance in paper)
%       6 - Three targets
%   Synthetic Data
%       0 - Emissive Example
%       1 - Two target example
%       2 - Four target example
%       3 - Same color example
%       4 - Non-Lambertian example
%
% Aquisition Time
%   30  - Extreme noise
%   60, 125, 250, 500
%   250 - Noisy
%   500, 1000, 2000, 4000
%   8000 - Minimal noise
%
% NOTE: Several parameters pertaining to the scene, and tuning parameters
% must be changed within the function
%
% % % % % % % % % % % % % % % % % % % % %

%% Run all paper examples, in order
clc
clear all
doorway_camera(0,0,8000,2,0.15,64,2e-5,1,1);
doorway_camera(0,1,8000,2,0.15,64,1e-5,1,1);
doorway_camera(0,2,8000,2,0.15,64,5e-7,1,1);
doorway_camera(1,0,0,1,0.15,64,5e-6,1,1);
doorway_camera(1,1,0,1,0.15,64,5e-6,1,1);
doorway_camera(1,2,0,1,0.15,64,5e-6,1,1);
doorway_camera(1,3,0,1,0.15,64,5e-6,1,1);
doorway_camera(1,4,0,1,0.15,64,1e-6,1,1);
doorway_camera(1,5,0,1,0.15,64,5e-7,1,1);
doorway_camera(1,6,0,2,0.15,64,1e-6,1,1);

%% New Synthetic examples
doorway_camera(0,3,8000,2,0.15,64,1e-6,1,1); % Same color objects
doorway_camera(0,4,8000,2,0.15,64,1e-5,1,1); % Non-lambertian objects

%% Intentional Model mismatch
% object is known 2m tall
doorway_camera(0,1,8000,1,0.15,64,1e-5,1,1);
doorway_camera(0,1,8000,3,0.15,64,1e-5,1,1);

%% Degrading image quality (synthetic, noise)
samps = [15,30,60,125,250,1000,4000,8000];
lambda = [1e-5,1e-5,1e-5,1e-5,5e-4,1e-5,1e-5,1e-5];
for i = 4:5
    doorway_camera(0,1,samps(i),2,0.15,64,lambda(i),1,1);
end

%% Changing FOV size
fov = [.5,.3,.15,.075,.0375];
lambda = [1e-6,1e-5,1e-5,1e-5,1e-5];
for i = 1:5
    doorway_camera(0,1,8000,2,fov(i),64,lambda(i),1,1);
end

%% Changing Resolution
res = [4,8,16,32,64,128];
lambda = [5e-8,2e-7,1e-6,3e-6,1e-5,1e-5];
for i = 1
    doorway_camera(0,1,8000,2,0.15,res(i),lambda(i),1,1);
end