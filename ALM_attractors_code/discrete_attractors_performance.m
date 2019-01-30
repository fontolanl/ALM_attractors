% DISCRETE_ATTRACTORS_PERFORMANCE.m generates the performance plots of EDF 9 in:
% Inagaki, Fontolan, Romani, Svoboda, Nature (2019)
%
% written by Lorenzo Fontolan

%% Clear var space and closes all figs

clear
close all

rng('shuffle')

%% Network params

N_trials = 1000;

% Membrane time constant
tauc = 100;
taun = 100;

dt = 1;

t_end = 6000;
t_vec = [0:dt:t_end];

taus = [tauc, taun, taun, tauc]';

bsln_t = 1500/dt;

t_dur_stim = 1000/dt;
t_dur_PV = 600/dt;
t_dur_delay = 2000/dt;

% sample epoch start and end (coincides with start of delay epoch)
t_stim_start = bsln_t + 1;
t_stim_end = t_stim_start + t_dur_stim - 1;

% Photostimulation
t_PV_start = t_stim_end + 1;
t_PV_end = t_PV_start + t_dur_PV - 1;

% End of Delay epoch
t_delay_end = t_stim_end + t_dur_delay;

t_ramp_start = t_stim_start;

alpha = 4.5;

%% network connectivity and parameters

% 1)'one_hemi_multi_fp_step'
% 2) 'one_hemi_single_fp';
% 3) 'two_hemi_ramp'
% 4) 'two_hemi_internal'

network_str = 'two_hemi_ramp';

switch network_str
    
    case 'one_hemi_multi_fp_step'
        % one hemi
        stim_sigma = 0.1;
        stim_amp = 0.8;
        
        % amplitude of fast noise
        sigma_noise = 9;
        
        % synaptic weights
        jCC = 5.8;
        jNN = jCC;
        jCI = 5;
        jIC = 1.08;
        jNI = jCI;
        jIN = jIC;
        jCN = 0.9;
        jNC = jCN;
        jII = 2;
        
        % input currents
        iC = 2.27;
        iN = iC;
        iI = 1;
        iI_bsl = iI;
        
        % cross-hemisphere excitation
        JH = 0.0;
        
        % static nonlinearity parameters
        tauD = 0.14;
        tauf = 0.8;
        U = 0.05;
        
        % Step-like ramp
        t_ramp_end = t_ramp_start + 10/dt;
        
        %amplitude of ramping at the end of delay
        ramp_amp = 2;
        
        %amplitude of PV perturbation
        PV_vec = [0,1,2,6];
        
    case 'one_hemi_single_fp'
        % one hemi
        stim_sigma = 0.05;
        stim_amp = 0.5;
        
        % amplitude of fast noise
        sigma_noise = 14;
        
        % synaptic weights
        jCC = 16;
        jNN = jCC;
        jCI = 8;
        jIC = 4;
        jNI = jCI;
        jIN = jIC;
        jCN = 1;
        jNC = jCN;
        jII = 1.13;
        
        % input currents
        iC = 16;
        iN = iC;
        iI = 1.07;
        iI_bsl = iI;
        
        % cross-hemisphere excitation
        JH = 0.0;
        
        % static nonlinearity parameters
        tauD = 0.0;
        tauf = 0.05;
        U = 0.05;
        
        % Step-like ramp
        t_ramp_end = t_ramp_start + 10/dt;
        
        %amplitude of ramping at the end of delay
        ramp_amp = 2;
        
        %amplitude of PV perturbation
        PV_vec = [0,0.4,2,12];
        
    case 'two_hemi_ramp'
        
        % amplitude of fast noise
        sigma_noise = 8;
        
        % synaptic weights
        jCC = 6.9;
        jNN = jCC;
        jCI = 5;
        jIC = 2.64;
        jNI = jCI;
        jIN = jIC;
        jCN = 2.2;
        jNC = jCN;
        jII = 2;
        
        % input currents
        iC = 1;
        iN = iC;
        iI = 1.2;
        iI_bsl = iI;
        
        % cross-hemisphere excitation
        JH = 0.58;
        
        % static nonlinearity parameters
        tauD = 0.1;
        tauf = 0.8;
        U = 0.1;
        
        ramp_amp = 4;
        
        ramp_type = '2s';
        
        switch ramp_type
            case '2s'
                % stimulus amplitude and variance
                stim_sigma = 0.05;
                stim_amp = 0.6;
                t_ramp_end = t_ramp_start + t_dur_delay + t_dur_stim;
                PV_vec = [0,0.5,1,3.2];
                
            case 'step'
                % stimulus amplitude and variance
                stim_sigma = 0.01;
                stim_amp = 1;
                t_ramp_end = t_ramp_start + 10/dt;
                PV_vec = [0,1.5,2.5,3.8,6.8];
                
        end
        
    case 'two_hemi_internal'
        
        stim_sigma = 0.01;
        stim_amp = 0.4;
        
        % amplitude of fast noise
        sigma_noise = 9;
        
        % change presample epoch duration
        
        bsln_t = 1500/dt;
        
        t_dur_stim = 1000/dt;
        t_dur_PV = 600/dt;
        t_dur_delay = 2000/dt;
        
        % sample epoch start and end (coincides with start of delay epoch)
        t_stim_start = bsln_t + 1;
        t_stim_end = t_stim_start + t_dur_stim - 1;
        
        % Photostimulation
        t_PV_start = t_stim_end + 1;
        t_PV_end = t_PV_start + t_dur_PV - 1;
        
        % End of Delay epoch
        t_delay_end = t_stim_end + t_dur_delay;
        t_ramp_start = t_stim_start;
        
        % synaptic weights
        jCC = 9.75;
        jNN = jCC;
        jCI = 2.5;
        jIC = 4;
        jNI = jCI;
        jIN = jIC;
        jCN = 2.5;
        jNC = jCN;
        jII = 1;
        
        % input currents
        iC = 2.3;%2;
        iN = iC;
        iI = 1;
        iI_bsl = iI;
        
        % cross-hemisphere excitation
        JH = 0.4;
        
        % static nonlinearity parameters
        tauD = 0.12;
        tauf = 0.1;
        U = 0.1;
        
        alpha = 2.5;
        
        ramp_amp = 3.15;
        
        t_ramp_end = t_ramp_start + 10/dt;
        
        PV_vec = [0,1.3,2.3,7];
        
end

%% Perturbation type

PV_type_vec{1} = 'ipsilateral';
PV_type_vec{2} = 'bilateral';

for k = 1:2
    
    PV_type = PV_type_vec{k};
    
    %% stimulus smooth
    
    sim_len = t_end/dt;
    fs = 1000/dt;
    tau_ext = 20;
    
    stim_t = zeros(sim_len,1);
    stim_t(t_stim_start+1+tau_ext/dt:t_stim_end-tau_ext/dt) = ones(length(1:(t_stim_end-t_stim_start)-2*tau_ext/dt),1);
    
    stim_smooth_vec = expsmooth(stim_t, fs, tau_ext);
    
    %% PV perturbation
    
    stim_t = zeros(sim_len,1);
    stim_t(t_PV_start+1:t_PV_end-tau_ext/dt) = ones(length(1:(t_PV_end-t_PV_start)-tau_ext/dt),1);
    
    PV_smooth_vec = expsmooth(stim_t, fs, tau_ext);
    
    %% ramping stimulus
    
    Ramp_temp = ramp_amp./(length(t_ramp_start:t_ramp_end-1)/1000)*[1,1,1,1];
    
    Ramp_vec = zeros(sim_len,4);
    for t_idx = t_ramp_start:t_ramp_end
        Ramp_vec(t_idx,:) = Ramp_temp'.*(t_idx-t_ramp_start)/1000;
    end
    Ramp_vec(t_ramp_end+1:end,:) = Ramp_vec(t_ramp_end,1);
    
    %% Dynamics
    
    PV_idx = 0;
    h_mat = zeros(length(PV_vec),N_trials,length(t_vec),4);
    r_mat = zeros(length(PV_vec),N_trials,length(t_vec),4);
    
    for PV_int = PV_vec
        PV_idx = PV_idx + 1;
        
        %% Scaled weights
        
        JCC = jCC - (jIC*jCI)/(1+jII);
        JNN = jNN - (jIN*jNI)/(1+jII);
        JCN = jCN - (jCI*jIN)/(1+jII);
        JNC = jNC - (jNI*jIC)/(1+jII);
        IC = iC - (jCI*iI)/(1+jII);
        IN = iN - (jNI*iI)/(1+jII);
        
        %% Synaptic Matrix
        
        p_matrix = [JCC, JCN, 0, 0;...
            JNC, JNN, 0, 0;...
            0, 0, JNN, JNC;...
            0, 0, JCN, JCC;];
        
        input_mat = zeros(4,t_end/dt,N_trials);
        
        % Noise
        noise_input = zeros(4,t_end/dt,N_trials);
        
        %% Dynamics
        for k_trial = 1:N_trials
            
            Stim_amp_distr(k_trial) = (stim_amp+stim_sigma*randn(1));
            
            ramp_in = Ramp_vec;
            
            if k_trial <= N_trials/2
                stim_in = Stim_amp_distr(k_trial).*[1,0,1,0].*stim_smooth_vec; % Contra trials
                if strcmp(network_str,'one_hemi_single_fp')
                    ramp_in(:,2) = 0;
                end
            else
                stim_in = Stim_amp_distr(k_trial).*[0,1,0,1].*stim_smooth_vec; % Ipsi trials
                if strcmp(network_str,'one_hemi_single_fp')
                    ramp_in(:,1) = 0;
                end
            end
            
            PV_iI = iI+PV_int*iI*PV_smooth_vec;
            
            % Scaled weights
            ICn = iC - (jCI*PV_iI)/(1+jII);
            INn = iN - (jNI*PV_iI)/(1+jII);
            
            switch PV_type
                
                case 'bilateral'
                    PV_in = [ICn,INn, INn, ICn];
                    
                case 'ipsilateral'
                    PV_in = [ICn,INn,IN*ones(sim_len,1),IC*ones(sim_len,1)];
                    
                    
                case 'contralateral'
                    PV_in = [IC*ones(sim_len,1),IN*ones(sim_len,1),INn,ICn,];
                    
            end
            
            % dynamics
            y = zeros(sim_len,4);
            h = zeros(sim_len,4);
            r = zeros(sim_len,4);
            h(1,:) = [1,1,1,1];
            
            for t_idx = 1:t_end/dt
                y(t_idx,:) = alpha.*log(1+exp(h(t_idx,:)./alpha));
                noise_input(:,t_idx,k_trial) = sigma_noise.*randn(4,1);
                
                input_mat(:,t_idx,k_trial) = p_matrix*[h_static(y(t_idx,:),tauD,tauf,U).*y(t_idx,:)]'+...
                    JH.*h_static([y(t_idx,3),y(t_idx,4),y(t_idx,1),y(t_idx,2)],tauD,tauf,U)'.*[y(t_idx,3),y(t_idx,4),y(t_idx,1),y(t_idx,2)]'...
                    + PV_in(t_idx,:)' + noise_input(:,t_idx,k_trial) + stim_in(t_idx,:)' + ramp_in(t_idx,:)';
                
                h(t_idx+1,:) = h(t_idx,:) + (-h(t_idx,:) + input_mat(:,t_idx,k_trial)')./taus'*dt;
                r(t_idx+1,:) = h_static(y(t_idx,:),tauD,tauf,U).*alpha.*log(1+exp(h(t_idx,:)./alpha));
            end
            
            %% store values
            
            h_mat(PV_idx,k_trial,:,:) = h;
            r_mat(PV_idx,k_trial,:,:) = r;
            
        end
    end
    
    r_mat_cell{k} = r_mat;           
 
end

%% Compute CD from unperturbed trials

t_vec_plt = (t_vec - t_delay_end*dt)/1000;

right_trials = 1:N_trials/2;
left_trials = N_trials/2+1:N_trials;

N_right_trials = length(right_trials);
N_left_trials = length(left_trials);

%% Performance barplot

figure
hold on

count_idx = 0;
for k = 1:2
    if k == 1
        for j=1:2
            correct_trials_right = find(r_mat_cell{k}(j,right_trials,t_delay_end,1)>r_mat_cell{k}(j,right_trials,t_delay_end,2));
            correct_trials_left = N_trials/2 + find(r_mat_cell{k}(j,left_trials,t_delay_end,1)<r_mat_cell{k}(j,left_trials,t_delay_end,2));
            correct_rate(j) = (length(correct_trials_right) + length(correct_trials_left))/N_trials;
            count_idx = count_idx + 1;
            bar(count_idx,100*correct_rate(j))
        end
        
    else
        for j=[2,4]
            correct_trials_right = find(r_mat_cell{k}(j,right_trials,t_delay_end,1)>r_mat_cell{k}(j,right_trials,t_delay_end,2));
            correct_trials_left = N_trials/2 + find(r_mat_cell{k}(j,left_trials,t_delay_end,1)<r_mat_cell{k}(j,left_trials,t_delay_end,2));
            correct_rate(j) = (length(correct_trials_right) + length(correct_trials_left))/N_trials;
            count_idx = count_idx + 1;
            bar(count_idx,100*correct_rate(j))
        end
    end
end

xticks([1:count_idx])
xticklabels({'Unperturbed','Unilateral Weak','Bilateral Weak','Bilateral Strong'})
set(gca,'TickDir','out','Ticklength',[0.025 0.025],'fontsize',10)
ylabel('% of Correct Trials')
title('Performance')
ylim([0,100])
%xlim([0,size(N_ttypes,2)+0.5])
box off

