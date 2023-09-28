function setupIzhikL5pnNetwork
% Izhikevich Quadratic Neuron spiking Modelling apical amplification
%clear
close all
clc

%%
% ========================================================================
% MODEL PARAMETERS
% ========================================================================

% Number of N x N lattice elements
params.N_grid = 70;
params.Ncells = params.N_grid^2;
% Number of time steps
params.nT = 1000*1; %20s sim plus 5 second transient
params.gs = 70/2; % BM: Example regimes 60 for NAd/10 for ach
params.threshold = 3; %normalised between -3 to 3 burst ratio
%BM: for no bursting 3, for all bursting -3 and for matched experimental bursting NAD/ACH ~ 0.3

%% Connectivity
WS = generateMexHatConnectivity(params.N_grid);
%BM: I suggest saving the connectivity

%% Somatic input
 tic
 allIcore = 9*randn(params.Ncells,params.nT);
 toc

%BM: I suggest saving a giant poisson input for reproducibility


%% Apical input and spatiotemporal smoothing
tic
apicalInput = randn(params.N_grid,params.N_grid,params.nT);
toc

%BM: I suggest saving a giant poisson input for reproducibility
%Integrator gaussian

tic
apicalInput = generateApical(params,apicalInput);
toc

%% Summing the input over 25ms windows apical dendritic timescale see MS
fbox = ones(1,1,25); % Rectangular summator 

intgrator=convn(apicalInput,fbox,'same');

intgrator = intgrator./std(intgrator(:)); %normalising the data

%Appending extra space of intgrator == timescale apical dendrites view for
%spike effects
intgrator(:,:,end+1:end+25) = zeros(params.N_grid,params.N_grid,25);

%%

tic
[firingsIndexes,firing_count,burstIndexes] = runIzhikL5pnNetwork(params,WS,intgrator,allIcore,plotflag);

filname = ['firings', num2str(PBSID)];
save(filname,'firingsIndexes','firing_count','burstIndexes','params')

toc


end




