% CONT_ATTRACTORS.m generates the continuous attractor dynamical systems 
% of EDF 1 in:
% Inagaki, Fontolan, Romani, Svoboda, Nature (2019)
%
% written by Lorenzo Fontolan

%% Clear var space and closes all figs

clear
close all

rng('shuffle')

%% Set parameters

%total number of trials
N_trials = 1000;

%%% synaptic weights
jCC = 150;
jCN = 150;
jNC = 150;
jNN = 150;
jII = 600;
jCI = 300;
jIC = 300-1.5;
jNI = 300+1.99;
jIN = 550;

% DC inputs to the 3 pops
iC = 0;
iN = 0;
iI = 0;
iI_bsl = iI;


%% Membrane time constants in ms
tauc = 20;
taun = 20;
taui = 10;

%% Synaptic time constants in ms
taucc = 100;
taunn = 100;
taucn = 100;
taunc = 100;
tauci = 10;
tauic = 25;
tauni = 10;
tauin = 25;
tauii = 10;

% time constants vector
p.taus = [tauc, taun, taui, taucc, taucn, taunc, taunn, tauci, tauic, tauni, tauin, tauii]';


%% Synaptic Matrix
p.W = [[0,0,0,jCC,jCN,0,0,-jCI,0,0,0,0];...
    [ 0,0,0,0,0,jNC,jNN,0,0,-jNI,0,0];...
    [ 0,0,0,0,0,0,0,0,jIC,0,jIN,-jII];...
    [ 1,0,0,0,0,0,0,0,0,0,0,0];...
    [ 0,1,0,0,0,0,0,0,0,0,0,0];...
    [ 1,0,0,0,0,0,0,0,0,0,0,0];...
    [ 0,1,0,0,0,0,0,0,0,0,0,0];...
    [ 0,0,1,0,0,0,0,0,0,0,0,0];...
    [ 1,0,0,0,0,0,0,0,0,0,0,0];...
    [ 0,0,1,0,0,0,0,0,0,0,0,0];...
    [ 0,1,0,0,0,0,0,0,0,0,0,0];...
    [ 0,0,1,0,0,0,0,0,0,0,0,0]];

%% time params

dt = 1;

% t_end = 4500/dt+bsln_t+t_response;
t_end = 5000/dt;
t_vec = [0:dt:t_end];

%% external inputs timepoint

bsln_t = 500/dt;

t_dur_stim = 1000/dt;
t_dur_PV = 600/dt;
t_dur_delay = 2000/dt;

% sample epoch start and end (coincides with start of delay epoch)
t_stim_start = 1000/dt + bsln_t;
t_stim_end = t_stim_start + t_dur_stim;

% Photostimulation
t_PV_start = 2000/dt + bsln_t;
t_PV_end = t_PV_start + t_dur_PV;

% End of Delay epoch
t_delay_end = (t_stim_end + t_dur_delay);

% decay time for smoothing
tauexp = 20/dt;

%% input current vectors

p.DC = [iC,iN,iI,0,0,0,0,0,0,0,0,0]';

p.PV_in = [0,0,1,0,0,0,0,0,0,0,0,0]';

p.stim_in = [1,0,0,0,0,0,0,0,0,0,0,0]';

% intensity of photostimulation
PV_vec = [0,150,300,600];

%% stimulus

% stim amplitude and variance
stim_strength = 3;
stim_sigma = 0.05;

% noise amplitude and autocorrelation decay time
sigma_noise = 5;
tau_noise = 5;

% options for integrating the ODE
opts = odeset('RelTol',1e-3,'AbsTol',1e-3);

%% Dynamics

PV_idx = 0;

for PV_intensity = PV_vec
    PV_idx = PV_idx + 1;
    
    disp(['Perturbation intensity = ', num2str(PV_intensity)])
        
    for k_trial = 1:N_trials
        
        if mod(k_trial,10)==0 %%% display the trial number every 10 trials
            k_trial
        end
        
        stim_amp_distr(PV_idx, k_trial) = (stim_strength+stim_sigma*randn(1));
        
        %% set the correct sign of the perturbation for right lick and left lick trials
        if k_trial <= N_trials/2
            p.stim_vec = @(t)(stim_amp_distr(PV_idx,k_trial).*((t>t_stim_start).*(t<=t_stim_end).*(1-exp(-(t-t_stim_start)/tauexp))+...
                (t>(t_stim_end)).*exp(-(t-(t_stim_end))/tauexp)));
            p.PV_pert_vec = @(t)(PV_intensity.*((t>t_PV_start).*(t<=t_PV_start+t_dur_PV).*(1-exp(-(t-t_PV_start)/tauexp))+...
                (t>(t_PV_start+t_dur_PV)).*exp(-(t-(t_PV_start+t_dur_PV))/tauexp)));
        else
            p.stim_vec = @(t)(-stim_amp_distr(PV_idx,k_trial).*((t>t_stim_start).*(t<=t_stim_end).*(1-exp(-(t-t_stim_start)/tauexp))+...
                (t>(t_stim_end)).*exp(-(t-(t_stim_end))/tauexp)));
            p.PV_pert_vec = @(t)(-PV_intensity.*((t>t_PV_start).*(t<=t_PV_start+t_dur_PV).*(1-exp(-(t-t_PV_start)/tauexp))+...
                (t>(t_PV_start+t_dur_PV)).*exp(-(t-(t_PV_start+t_dur_PV))/tauexp)));
        end
        
        % initial conditions
        h_in = sigma_noise.*0.00175.*randn(1,12);
        
        % noise vector
        noise_vec = zeros(3,length(t_vec));
        
        for t_idx = 1:t_end
            noise_vec(:,t_idx+1) = noise_vec(:,t_idx) + dt.*(-noise_vec(:,t_idx) + ...
                sigma_noise*randn(3,1))./tau_noise;
        end

        % ODE
        [t,h] = ode45(@(t,h) ode_cont_attractor(t, h, p, t_vec, noise_vec'), t_vec, h_in, opts);
        
        % only store currents
        h_mat(PV_idx,k_trial,:,:) = h(:,1:3);
        
    end

end

%% Plot figures

[correct_trials_right, correct_trials_left, h_right_pj_c, h_left_pj_c] = ALM_attractor_figs(h_mat, ...
    PV_vec, t_vec, dt, t_stim_start, t_stim_end, t_delay_end, N_trials );

%% Plot phase line
fr_res = 0.05;
win_ms_pp = 100;

derivative_vs_CD(h_right_pj_c{1}, h_left_pj_c{1}, length(correct_trials_right{1}),...
    length(correct_trials_left{1}), t_delay_end, win_ms_pp, fr_res, dt )
