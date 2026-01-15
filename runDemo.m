% MAIN

clear all, clc, close all

tic
b = 5; 
kappa = 0.3; 
lambda = b/kappa;

%% Choose geometry-type

geometry = 'crinkled';

%funil, hourglass, crinkled, asymHourglass, swift

%% Choose model

% nonlinear only for now, (linear eventually)

%% Choose scheme??

% --------
parameterStation 

%---

%setDomain(b, geometry,domainOptions);

%---
%load(['DomainData/geometry= ', geometry, '.mat'], 'w', 'J','P');
%evolveWave(kappa, b, geometry, domainOptions, waveOptions);

tiMe = toc;

%---

processWaveData(kappa, b, geometry, waveViewOptions)


