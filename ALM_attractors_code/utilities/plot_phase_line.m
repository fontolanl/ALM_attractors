function [f1, var_proj] = plot_phase_line(t_end,t_axis,mean_proj,modes)
% PLOT_PHASE_LINE computes and plots phase line in Inagaki, Fontolan,
% Romani, Svoboda, Nature (2019)

    x_ax = -10:0.2:10;

    num_btstrap = size(mean_proj{1},1);
    num_modes   = numel(modes);   
    
    for n=1:num_btstrap
        proj_per_tmp = nan(num_modes,numel(t_axis));
        for md = 1:num_modes
            proj_per_tmp(md,:) = mean_proj{modes(md)}(n,:);
        end

        % speed analysis
        [speed_all_per_tmp,speed_all_sh_per_tmp] = speed_analysis_each_data_for_per(x_ax,proj_per_tmp(:,t_axis<-0.05 & t_axis>t_end),0); 
        [speed_all_per_R_tmp,speed_all_sh_per_R_tmp] = speed_analysis_each_data_for_per(x_ax,proj_per_tmp(1:num_modes/2,t_axis<-0.05 & t_axis>t_end),0);
        [speed_all_per_L_tmp,speed_all_sh_per_L_tmp] = speed_analysis_each_data_for_per(x_ax,proj_per_tmp(num_modes/2+1:num_modes,t_axis<-0.05 & t_axis>t_end),0);

        speed_all_per(n,:)    = speed_all_per_tmp';
        speed_all_sh_per(n,:) = speed_all_sh_per_tmp';
        
        speed_all_per_R(n,:)    = speed_all_per_R_tmp';
        speed_all_sh_per_R(n,:) = speed_all_sh_per_R_tmp';
        
        speed_all_per_L(n,:)    = speed_all_per_L_tmp';
        speed_all_sh_per_L(n,:) = speed_all_sh_per_L_tmp';

        % variance analysis
        var_proj{1}(n,:) = var(proj_per_tmp(1:num_modes/2,:),[],1); % contra trials
        var_proj{2}(n,:) = var(proj_per_tmp(num_modes/2+1:num_modes,:),[],1); % ipsi trials

    end
    
    %% normalize the dt to be 10ms
    diff_t_axis = t_axis(3)-t_axis(2);
    norm_factor = 0.01/diff_t_axis;
    
    %% plot phase portrait
    f1 = figure;
    hold on  
    [max_y1 min_y1] = plot_trace(x_ax,speed_all_per*norm_factor);  
    plot(x_ax,mean(speed_all_sh_per),'k:')
    ylim([min_y1 max_y1]);xlim([-0.1 1.1])
    xlabel('Projection to mode, spikes per s')
    ylabel('d(Projection to mode)/dt')
    title('Phase Line')
    set(gca,'TickDir','out','Ticklength',[0.025 0.025],'box','off','fontsize',8)
    set(gca,'TickDir','out','Ticklength',[0.025 0.025],'fontsize',18)
    hold off
    
     
end


function [speed_all_pooled,speed_all_sh_pooled] = speed_analysis_each_data_for_per(x_ax,data,thr)
 
 
    mean_all_ori  = data(:,2:end-1);
    speed_all_ori = (data(:,3:end)-data(:,1:end-2))/2;
    
    int = x_ax(2) - x_ax(1);
    
    
    % shuffle speed in trial direction
    num_A = size(speed_all_ori,1);
    num_T = size(speed_all_ori,2);
    
    for t=1:num_T
        speed_all_sh(:,t) = speed_all_ori(randperm(num_A),t);
    end

    mean_all = reshape(mean_all_ori,1,[]);
    speed_all = reshape(speed_all_ori,1,[]);
    speed_all_sh = reshape(speed_all_sh,1,[]);

    % speed density

   
    speed_all_pooled = nan(length(x_ax),1);
    for x = 1:length(x_ax)
        idx = mean_all>x_ax(x)-int/2 & mean_all<=x_ax(x)+int/2;
        if sum(idx)>thr
            speed_all_pooled(x) = nanmean(speed_all(idx));
        end
    end

    
    speed_all_sh_pooled = nan(length(x_ax),1);
    for x = 1:length(x_ax)        
        idx = mean_all>x_ax(x)-int/2 & mean_all<=x_ax(x)+int/2;
        if sum(idx)>thr
            speed_all_sh_pooled(x) = nanmean(speed_all_sh(idx));
        end
    end
   
end


function [max_y min_y]=plot_trace(x_ax,data)
    hold on 
    mean_data = data;
    %sem       = std(data);
    mask      = ~isnan(mean_data);

    %fill([x_ax(mask),fliplr(x_ax(mask))],[mean_data(mask) + sem(mask),fliplr(mean_data(mask) - sem(mask))],[0.5 0.5 0.5],'edgecolor','none')
    plot(x_ax(mask), mean_data(mask),'k')

    max_y = max(mean_data(20:end-20)) *1.05;
    min_y = min(mean_data(20:end-20));
    if min_y<0
        min_y = min_y *1.05;
    else
        min_y = min_y*0.95; 
    end


    set(gca,'fontsize',8,'tickdir','out','TickLength',[0.025, 0.01])  
   
end



