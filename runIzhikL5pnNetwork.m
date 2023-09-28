function [firingsIndexes,firing_count,burstIndexes] = runIzhikL5pnNetwork(params,WS,intgrator,allIcore,plotflag)
%input
% params.N_grid grid size
% params.Ncells number of cells in the simulation
% params.nt simulation length in ms

%output

%% Grab parameter values
N_grid = params.N_grid;
Ncells = params.Ncells;
nT = params.nT;
threshold = params.threshold; % 90% at sd of 1.7 as intgrator without addition normalised

%%

transientTime = 0; %start recording values after Xms suggest > 1s
isplot = plotflag;

%% Initialisation

%Izhik model parameters
a=0.02*ones(Ncells,1);
b=0.2*ones(Ncells,1);
c=-65*ones(Ncells,1); % reseting voltage for the fired neurons
d=8*ones(Ncells,1); % reseting u for the fired neurons

%starting all neurons at -65mV
v=-65*ones(Ncells,1);  % Initial values of v
u=b.*v;               % Initial values of u

%Storage variables
voltage=ones(N_grid,N_grid);
firingsIndexes=[];
burstIndexes=[];
firing_count=nan(1,nT-transientTime);



%% Running the model

for tt=1:nT
    
    %Display simulation time
    if mod(tt,100)==0
        disp(['t=',num2str(tt)]);
    end
    
    
    % Find spiking neurons
    fired=find(v>=30); % indices of spikes   
    
    
    %% Modify connectivity based on spike timing
    %fired20msCollection = [fired20msCollection(:); 
    for ff = 1:numel(fired)
        
    end
    
    %% Matrix/Integrator calculations
    %Finding the cells which have turned into burst mode
    cellBurst = find( intgrator(:,:,tt) > threshold);
    burstIndexes = [burstIndexes; tt+0*cellBurst,cellBurst];
    %Changing their parameters into chattering mode
    bI = cellBurst;
    c(bI)= -55; % reseting voltage for the fired neurons
    d(bI)= 4; % reseting u for the fired neurons
    
    %IB: d=4, c = -55
    %IB: Intrinsically bursting neurons, if stimulated with a long pulse of dc current, fire an initial burst of spikes followed by shorter bursts, and then tonic spikes (Connors and Gutnick 1990). These are predominantly pyramidal neurons in layer 5. 

    
    %% Core input
    % Thalamic/core input poisson noise
     Icore = allIcore(:,tt);
    
    if ~isempty(fired)  % is not empty fired
        Icore=Icore+sum(WS(:,fired),2);
        v(fired)=c(fired);
        u(fired)=u(fired)+d(fired);
    end
    
    %% Numerical differentiation
    v=v+0.5*(0.04*v.^2+5*v+140-u+Icore); %step 0.5ms
    v=v+0.5*(0.04*v.^2+5*v+140-u+Icore); %for stability
    u=u+a.*(b.*v-u);
    
    
    %% Storage of values (at the moment just spikes could also be membraneVoltages)
    if tt > transientTime
        firingsIndexes = [firingsIndexes; tt+0*fired,fired];
        firing_count(tt-transientTime)=numel(fired);
        %         vt = v;
        %         vt(vt>10) = 10;
        %         %voltage(:,:,tt) = vt;
        %         voltage(:,:) = vt;
    end
    
    %% Optional Plotting
    if (mod(tt,5)==0) && isplot
        
        vt = reshape(v,[N_grid,N_grid]);
        vt(vt>10) = 10;
        voltage(:,:) = vt;
        
        clf
        %subplot(1,2,1);
        imagesc(voltage)
        title(['V (voltage)', 'Time: ', num2str(tt), ' ms'])
        %colormap(jet)
        caxis([-65 10])
        axis([0,N_grid,0,N_grid]);
        axis square;
        
%This shows the apical input        
%         subplot(1,2,2);
%         imagesc(reshape(d,[N_grid,N_grid]))
%         caxis([2 8])
%         title('Chattering (Blue) & Regular spiking (Red)')
%         axis([0,N_grid,0,N_grid]);
%         axis square;
         
%         imagesc(reshape(u,[100,100]))
%         colormap(jet)
%         title('Recovery variable u')
%         %imagesc(reshape(c,[N_grid,N_grid]))
%         %title('Chattering (Red) & Regular spiking (Blue)')
%         axis([0,N_grid,0,N_grid]);
%         axis square;
        
        drawnow
    end
    
    
    % Resetting back to RS
    c=-65*ones(Ncells,1); % reseting voltage for the fired neurons
    d=8*ones(Ncells,1); % reseting u for the fired neurons
    
end
end