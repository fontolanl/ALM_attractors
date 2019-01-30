function [correct_trials_right, correct_trials_left, r_right_pj_c, r_left_pj_c] = ALM_attractor_figs(r_mat,...
    PV_vec, t_vec, dt, t_stim_start, t_stim_end, t_delay_end, N_trials )
% ALM_attractor_figs allows to plot some of the model figures in Inagaki,
% Fontolan, Romani, Svoboda, Nature (2019)
%
% written by Lorenzo Fontolan

%% Trials

t_vec_plt = (t_vec - t_delay_end*dt)/1000;

bin3 = find(t_vec_plt<-1.25,1,'last');
bin4 = find(t_vec_plt<-0.05,1,'last');

right_trials = 1:N_trials/2;
left_trials = N_trials/2+1:N_trials;

N_right_trials = length(right_trials);
N_left_trials = length(left_trials);

correct_trials_right{1} = find(r_mat(1,right_trials,t_delay_end,1)>r_mat(1,right_trials,t_delay_end,2));
correct_trials_left{1} = N_trials/2 + find(r_mat(1,left_trials,t_delay_end,1)<r_mat(1,left_trials,t_delay_end,2));

%% Compute Coding direction

% smoothing window
win_ms = 100;

r_CD = [smooth(squeeze(mean(squeeze(r_mat(1,correct_trials_right{1},:,1))))...
    ,round(win_ms/dt))...
    - smooth(squeeze(mean(squeeze(r_mat(1,correct_trials_left{1},:,1)))),...
    round(win_ms/dt)),smooth(squeeze(mean(squeeze(r_mat(1,correct_trials_right{1},:,2)))),...
    round(win_ms/dt))...
    - smooth(squeeze(mean(squeeze(r_mat(1,correct_trials_left{1},:,2)))),...
    round(win_ms/dt))];

CD_delay = mean(r_CD((t_delay_end-400/dt):t_delay_end,:));

%% Compute endpoints and projected trajectories

r_right_c{1} = cat(3,squeeze(r_mat(1,correct_trials_right{1},:,1)),squeeze(r_mat(1,correct_trials_right{1},:,2)));
r_left_c{1} = cat(3,squeeze(r_mat(1,correct_trials_left{1},:,1)),squeeze(r_mat(1,correct_trials_left{1},:,2)));

r_right_pj_c{1} = zeros(length(correct_trials_right{1}),length(t_vec));
r_left_pj_c{1} = zeros(length(correct_trials_left{1}),length(t_vec));

endpoints_right = zeros(length(correct_trials_right{1}),1);
endpoints_left = zeros(length(correct_trials_left{1}),1);

for i = 1:length(correct_trials_right{1})
    
    endpoints_right(i) = mean(squeeze(r_right_c{1}(i,(t_delay_end-400/dt):t_delay_end,:))*CD_delay');
    r_right_pj_c{1}(i,:) = (squeeze(r_right_c{1}(i,:,:))*CD_delay')';
    
end

for i = 1:length(correct_trials_left{1})
    
    endpoints_left(i) = mean(squeeze(r_left_c{1}(i,(t_delay_end-400/dt):t_delay_end,:))*CD_delay');
    r_left_pj_c{1}(i,:) = (squeeze(r_left_c{1}(i,:,:))*CD_delay')';
    
end

% Median of endpoints
right_norm = median(endpoints_right);
left_norm = median(endpoints_left);

% Normalization
r_right_pj_c{1} = (r_right_pj_c{1} - left_norm)/(right_norm - left_norm);
r_left_pj_c{1} = (r_left_pj_c{1} - left_norm)/(right_norm - left_norm);

%% Smooth the vectors

r_smooth_left_pj{1} = smooth(mean(r_left_pj_c{1}),round(win_ms/dt));
r_smooth_right_pj{1} = smooth(mean(r_right_pj_c{1}),round(win_ms/dt));

r_smooth_left_pj_std{1} = smooth(std(r_left_pj_c{1}(:,:,1)),round(win_ms/dt));
r_smooth_right_pj_std{1} = smooth(std(r_right_pj_c{1}(:,:,1)),round(win_ms/dt));

Data_at_each_bin{1}(1,1) = r_smooth_right_pj{1}(bin3);
Data_at_each_bin{1}(1,2) = r_smooth_right_pj{1}(bin4);
Data_at_each_bin{1}(2,1) = r_smooth_left_pj{1}(bin3);
Data_at_each_bin{1}(2,2) = r_smooth_left_pj{1}(bin4);
mean_proj{1} = r_smooth_right_pj{1};
mean_proj{2} = r_smooth_left_pj{1};

%% generate the unperturbed correct trials fig
fig_h = figure;
set(fig_h, 'Position', [100, 100, 800, 800]);       % <? ?Height? Increased

nrow = 3;
ncol = 2;

subplot(nrow,ncol,1)
hold on
% fill([(t_vec - t_end_delay)/1000,fliplr((t_vec - t_end_delay)/1000)],...
%    [fr_vec_left-fr_vec_left_std,fliplr(fr_vec_left+fr_vec_left_std)],'r','linestyle','none');
shadedErrorBar(t_vec_plt, r_smooth_right_pj{1},...
    r_smooth_right_pj_std{1},'b',0)
shadedErrorBar(t_vec_plt, r_smooth_left_pj{1},...
    r_smooth_left_pj_std{1},'r',0)
plot(t_vec_plt, r_smooth_right_pj{1},'b','Linewidth',4)
plot(t_vec_plt, r_smooth_left_pj{1},'r','Linewidth',4)
plot(ones(2,1)*(t_stim_start - t_delay_end)*dt/1000,[-1,2],'k--','Linewidth',2)
plot(ones(2,1)*(t_stim_end - t_delay_end)*dt/1000,[-1,2],'k--','Linewidth',2)
plot(zeros(2,1),[-1,2],'k--','Linewidth',2)
set(gca,'TickDir','out','Ticklength',[0.025 0.025],'fontsize',8)
ylabel('Proj. to CD')
xlabel('Time to movement onset (s)')
box off
set(gca,'fontname','Arial','color','w','fontsize',14)
xticks([-4,-2,0])
xlim([-4.2,0.2])
ylim([-0.2,1.2])
title('Unperturbed correct trials')

%% Distribution of end points - unperturbed correct trials fig

[N_right,edges_right] = histcounts(r_right_pj_c{1}(:,t_delay_end),[-1.5:0.1:2.5]+0.05);

x_N_right = edges_right(1:end-1)+(edges_right(2)-edges_right(1))/2;

[N_left,edges_left] = histcounts(r_left_pj_c{1}(:,t_delay_end),[-1.5:0.1:2.5]+0.05);

x_N_left = edges_left(1:end-1)+(edges_left(2)-edges_left(1))/2;

subplot(nrow,ncol,2)
hold on
plot(flipud(x_N_right), N_right/N_right_trials,'b','Linewidth',3);
plot(flipud(x_N_left), N_left/N_left_trials,'r','Linewidth',3);
set(gca, 'xdir','reverse')
xticks([-1:2])
xticklabels({'','0','1'})
ylabel('Fraction of correct trials (unperturbed)')
xlabel('Proj. to CD')
xlim([-0.5,1.5])
title('Endpoint distribution - unperturbed correct trials')
set(gca,'TickDir','out','Ticklength',[0.025 0.025],'fontsize',8)
set(gca,'fontname','Arial','color','w','fontsize',14)
box off
view([90,90])

%% Projected Activity - Perturbed correct trials

subplot(nrow,ncol,3)
hold on
for j = 2:length(PV_vec)
    
    correct_trials_right{j} = find(r_mat(j,right_trials,t_delay_end,1)>r_mat(j,right_trials,t_delay_end,2));
    correct_trials_left{j} = N_trials/2 + find(r_mat(j,left_trials,t_delay_end,1)<r_mat(j,left_trials,t_delay_end,2));
    
    r_right_pj_c{j} = zeros(length(correct_trials_right{j}),length(t_vec));
    r_left_pj_c{j} = zeros(length(correct_trials_left{j}),length(t_vec));
    
    r_right_c{j} = squeeze(r_mat(j,correct_trials_right{j},:,1:2));
    r_left_c{j} = squeeze(r_mat(j,correct_trials_left{j},:,1:2));
    
    for i = 1:length(correct_trials_right{j})
        r_right_pj_c{j}(i,:) = (squeeze(r_right_c{j}(i,:,:))*CD_delay' - left_norm)/(right_norm - left_norm);
    end
    
    for i = 1:length(correct_trials_left{j})
        r_left_pj_c{j}(i,:) = (squeeze(r_left_c{j}(i,:,:))*CD_delay' - left_norm)/(right_norm - left_norm);
    end
    
    r_right_pj_c_avg{j} = smooth(mean(r_right_pj_c{j}),round(win_ms/dt));
    r_left_pj_c_avg{j} = smooth(mean(r_left_pj_c{j}),round(win_ms/dt));
    
    plot(t_vec_plt, r_left_pj_c_avg{j},...
        'Color',[1,0,j/(length(PV_vec)+3)],'Linewidth',4)
    plot(t_vec_plt, r_right_pj_c_avg{j},...
        'Color',[0,j/(length(PV_vec)+3),1],'Linewidth',4)
    
    
    Data_at_each_bin{j}(1,1) = r_right_pj_c_avg{j}(bin3);
    Data_at_each_bin{j}(1,2) = r_right_pj_c_avg{j}(bin4);
    Data_at_each_bin{j}(2,1) = r_left_pj_c_avg{j}(bin3);
    Data_at_each_bin{j}(2,2) = r_left_pj_c_avg{j}(bin4);
    
    mean_proj{2*(j-1)+1} = r_right_pj_c_avg{j};
    mean_proj{2*j}       = r_left_pj_c_avg{j};
    
end

plot(t_vec_plt, r_smooth_right_pj{1},'b--','Linewidth',4)
plot(t_vec_plt, r_smooth_left_pj{1},'r--','Linewidth',4)
plot(ones(2,1)*(t_stim_start - t_delay_end)*dt/1000,[-2,2],'k--','Linewidth',2)
plot(ones(2,1)*(t_stim_end - t_delay_end)*dt/1000,[-2,2],'k--','Linewidth',2)
plot(zeros(2,1),[-2,2],'k--','Linewidth',2)
set(gca,'TickDir','out','Ticklength',[0.025 0.025],'fontsize',8)
ylabel('Proj. to CD')
xlabel('Time to movement onset (s)')
box off
set(gca,'fontname','Arial','color','w','fontsize',14)
title('Perturbed correct trials')
xticks([-4,-2,0])
yticklabels({'-1','','0','','1'})
yticks([-1,-0.5,0,0.5,1])
xlim([-4.2,0.2])
ylim([-0.25,1.25])


%% Distribution of end points - Perturbed correct trials

r_right_pj_c_mat = cat(1,r_right_pj_c{2:end});
r_left_pj_c_mat = cat(1,r_left_pj_c{2:end});

[N_right,edges_right] = histcounts(r_right_pj_c_mat(:,t_delay_end),[-1.5:0.1:2.5]+0.05);

x_N_right = edges_right(1:end-1)+(edges_right(2)-edges_right(1))/2;

[N_left,edges_left] = histcounts(r_left_pj_c_mat(:,t_delay_end),[-1.5:0.1:2.5]+0.05);

x_N_left = edges_left(1:end-1)+(edges_left(2)-edges_left(1))/2;

%%
subplot(nrow,ncol,4)
hold on
plot(flipud(x_N_right), N_right/size(r_right_pj_c_mat,1),'b','Linewidth',3);
plot(flipud(x_N_left), N_left/size(r_left_pj_c_mat,1),'r','Linewidth',3);
set(gca, 'xdir','reverse')
xlim([-0.5,1.5])
xticks([-1:2])
xticklabels({'','0','1'})
ylabel('Fraction of correct trials')
xlabel('Proj. to CD')
title('Endpoint distribution - perturbed correct trials')
set(gca,'TickDir','out','Ticklength',[0.025 0.025],'fontsize',8)
set(gca,'fontname','Arial','color','w','fontsize',14)
box off
view([90,90])

%% Projected Activity - Perturbed error trials

subplot(nrow,ncol,5)
hold on
for j = 2:length(PV_vec)
    
    error_trials_right{j} = find(r_mat(j,right_trials,t_delay_end,1)<r_mat(j,right_trials,t_delay_end,2));
    error_trials_left{j} = N_trials/2 + find(r_mat(j,left_trials,t_delay_end,1)>r_mat(j,left_trials,t_delay_end,2));
    
    
    r_right_pj_e{j} = zeros(length(error_trials_right{j}),length(t_vec));
    r_left_pj_e{j} = zeros(length(error_trials_left{j}),length(t_vec));
    
    r_right_e{j} = squeeze(r_mat(j,error_trials_right{j},:,1:2));
    r_left_e{j} = squeeze(r_mat(j,error_trials_left{j},:,1:2));
    
    r_right_pj_e_avg{j} = zeros(length(t_vec),1);
    r_left_pj_e_avg{j} = zeros(length(t_vec),1);
    
    for i = 1:length(error_trials_right{j})
        r_right_pj_e{j}(i,:) = (squeeze(r_right_e{j}(i,:,:))*CD_delay' - left_norm)/(right_norm - left_norm);
    end
    
    for i = 1:length(error_trials_left{j})
        r_left_pj_e{j}(i,:) = (squeeze(r_left_e{j}(i,:,:))*CD_delay' - left_norm)/(right_norm - left_norm);
    end
    
    if length(error_trials_right{j})>1
        r_right_pj_e_avg{j} = smooth(mean(r_right_pj_e{j}),round(win_ms/dt));
        plot(t_vec_plt, r_right_pj_e_avg{j},...
            'Color',[0,j/(length(PV_vec)+3),1],'Linewidth',4)
        
    elseif length(error_trials_right{j})==1
        r_right_pj_e_avg{j} = smooth(r_right_pj_e{j},round(win_ms/dt));
    end
    
    if length(error_trials_left{j})>1
        r_left_pj_e_avg{j} = smooth(mean(r_left_pj_e{j}),round(win_ms/dt));
        plot(t_vec_plt, r_left_pj_e_avg{j},...
            'Color',[1,0,j/(length(PV_vec)+3)],'Linewidth',4)
    elseif length(error_trials_left{j})==1
        r_left_pj_e_avg{j} = smooth(r_left_pj_e{j},round(win_ms/dt));
    end
    
    Data_at_each_bin{j+3}(1,1) = r_right_pj_e_avg{j}(bin3);
    Data_at_each_bin{j+3}(1,2) = r_right_pj_e_avg{j}(bin4);
    Data_at_each_bin{j+3}(2,1) = r_left_pj_e_avg{j}(bin3);
    Data_at_each_bin{j+3}(2,2) = r_left_pj_e_avg{j}(bin4);
    
    mean_proj{2*(j+2)+1} = r_right_pj_e_avg{j};
    mean_proj{2*(j+3)}   = r_left_pj_e_avg{j};
    
    % decode pre stim
    bin_pre_stim = find(t_vec_plt>-2.05,1,'first');
    
    cor_right = r_right_pj_e{j}(:,bin_pre_stim)>0.5;
    cor_left   = r_left_pj_e{j}(:,bin_pre_stim)<0.5;
    
    r_right_temp = smooth(mean(r_right_pj_e{j}(cor_right,:)),round(win_ms/dt));
    r_left_temp = smooth(mean(r_left_pj_e{j}(cor_left,:)),round(win_ms/dt));
    
    right_data_tmp = r_right_temp(1:1/dt:end,:);
    left_data_tmp   = r_left_temp(1:1/dt:end,:);
    
    mean_cor_proj{2*(j+2)+1} = right_data_tmp;
    mean_cor_proj{2*(j+3)}   = left_data_tmp;
    
end

plot(t_vec_plt, r_smooth_right_pj{1},'b--','Linewidth',4)
plot(t_vec_plt, r_smooth_left_pj{1},'r--','Linewidth',4)
plot(ones(2,1)*(t_stim_start - t_delay_end)*dt/1000,[-2,2],'k--','Linewidth',2)
plot(ones(2,1)*(t_stim_end - t_delay_end)*dt/1000,[-2,2],'k--','Linewidth',2)
plot(zeros(2,1),[-2,2],'k--','Linewidth',2)
set(gca,'TickDir','out','Ticklength',[0.025 0.025],'fontsize',8)
ylabel('Proj. to CD')
xlabel('Time to movement onset (s)')
box off
set(gca,'fontname','Arial','color','w','fontsize',14)
title('Perturbed error trials')
xticks([-4,-2,0])
yticklabels({'-1','','0','','1'})
yticks([-1,-0.5,0,0.5,1])
xlim([-4.2,0.2])
ylim([-0.25,1.25])

%% Distribution of end points - perturbed error trials

r_right_pj_e_mat = cat(1,r_right_pj_e{:});
r_left_pj_e_mat = cat(1,r_left_pj_e{:});

[N_right,edges_right] = histcounts(r_right_pj_e_mat(:,t_delay_end),[-1.5:0.1:2.5]+0.05);

x_N_right = edges_right(1:end-1)+(edges_right(2)-edges_right(1))/2;

[N_left,edges_left] = histcounts(r_left_pj_e_mat(:,t_delay_end),[-1.5:0.1:2.5]+0.05);

x_N_left = edges_left(1:end-1)+(edges_left(2)-edges_left(1))/2;

subplot(nrow,ncol,6)
hold on
plot(flipud(x_N_right), N_right/size(r_right_pj_e_mat,1),'b','Linewidth',3);
plot(flipud(x_N_left), N_left/size(r_left_pj_e_mat,1),'r','Linewidth',3);
set(gca, 'xdir','reverse')
xlim([-0.5,1.5])
xticks([-1:2])
xticklabels({'','0','1'})
ylabel('Fraction of error trials')
xlabel('Proj. to CD')
title('Endpoint distribution - perturbed error trials')
set(gca,'TickDir','out','Ticklength',[0.025 0.025],'fontsize',8)
set(gca,'fontname','Arial','color','w','fontsize',14)
box off
view([90,90])

%% Trial-to-trial fluctuations

norm_var_left = std(r_left_pj_c{1}(:,:,1))/std(r_left_pj_c{1}(:,t_stim_end,1));
norm_var_right = std(r_right_pj_c{1}(:,:,1))/std(r_right_pj_c{1}(:,t_stim_end,1));

figure
hold on
plot(t_vec_plt,norm_var_left,'r','Linewidth',3)
plot(t_vec_plt,norm_var_right,'b','Linewidth',3)
set(gca,'TickDir','out','Ticklength',[0.025 0.025],'fontsize',18)
set(gca,'fontname','Arial','color','w','fontsize',14)
xlabel('Time to movement onset (s)')
ylabel({('Changes in across-trials fluc.') ; ('of proj. to CD')})
xticks([-1,0])
xlim([-1.2,0])
box off
title('Trial-averaged fluctuations')

%% Plot change in Delta

title_list = {'0.1mW correct','0.2mW correct','0.3mW correct','incorrect'};
figure
set(gcf,'Color',[1 1 1],'Position',[100 50 600 400])
for ty=1:4
    subplot(1,4,ty)
    hold on
    if ty<4
        data_1 = (abs(Data_at_each_bin{ty+1}(1,1) - Data_at_each_bin{1}(1,1))+abs(Data_at_each_bin{ty+1}(2,1) - Data_at_each_bin{1}(2,1)))/2;
        data_2 = (abs(Data_at_each_bin{ty+1}(1,2) - Data_at_each_bin{1}(1,2))+abs(Data_at_each_bin{ty+1}(2,2) - Data_at_each_bin{1}(2,2)))/2;
    else
        data_bin1_C_tmp =  nanmean([Data_at_each_bin{5}(1,1),Data_at_each_bin{6}(1,1),Data_at_each_bin{7}(1,1)]);
        data_bin1_I_tmp =  nanmean([Data_at_each_bin{5}(2,1),Data_at_each_bin{6}(2,1),Data_at_each_bin{7}(2,1)]);
        data_bin2_C_tmp =  nanmean([Data_at_each_bin{5}(1,2),Data_at_each_bin{6}(1,2),Data_at_each_bin{7}(1,2)]);
        data_bin2_I_tmp =  nanmean([Data_at_each_bin{5}(2,2),Data_at_each_bin{6}(2,2),Data_at_each_bin{7}(2,2)]);
        
        data_1 = (abs(data_bin1_C_tmp - Data_at_each_bin{1}(1,1))+abs(data_bin1_I_tmp - Data_at_each_bin{1}(2,1)))/2;
        data_2 = (abs(data_bin2_C_tmp - Data_at_each_bin{1}(1,2))+abs(data_bin2_I_tmp - Data_at_each_bin{1}(2,2)))/2;
        
    end
    title(title_list{ty})
    plot([data_1;data_2],'ko-','linewidth',5)
    ylim([0 1.4]);xlim([0.75 2.25]);set(gca,'Ytick',0:0.5:2,'Xtick',1:2)
    set(gca,'fontsize',8,'tickdir','out','TickLength',[0.025, 0.01])
    ylabel('Delta proj')
end

%% Plot phase line
for ty = 1:10
    if ty<9
        mean_proj_to_plot{ty} = mean_proj{ty}';
    else
        mean_proj_to_plot{ty} = nanmean([mean_proj{ty},mean_proj{ty+2},mean_proj{ty+4}],2)';
    end
end

modes = [1 3 5 7 10 2 4 6 8 9];

t_pl_start = -1.25; % start of time window to compute phase line in seconds

[f1, var_proj] = plot_phase_line(t_pl_start, t_vec_plt, mean_proj_to_plot, modes);

end

