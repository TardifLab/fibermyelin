%% Create synthetic data from the equation to run the simulations

%create a grid with different combinations of AFD and MTR


% now create an array of combinations of afds and mtrs
afd1=0:0.1:1;
afd2=0:0.1:1;
mtr1=0.2:0.05:0.4;
mtr2=0.2:0.05:0.4;

afd = zeros(length(afd1),length(mtr1)*length(mtr2),1,2); %time dimension is 2 (number of fibers)
mtr = zeros(length(afd1),length(mtr1)*length(mtr2),1,2); %time dimension is 2 (number of fibers)

%afd array: each row is a different afd combination, each column is the
%same
for a=1:length(afd1)
    afd(a,1:end,1,1) = afd1(a); 
    afd(a,1:end,1,2) = afd2(end+1-a);
end

count=1;
for m1=1:length(mtr1)
    for m2=1:length(mtr2)
        mtr(1:end,count,1) = mtr1(m1);
        mtr(1:end,count,2) = mtr2(m2);
        count=count+1;
    end
end

%create a new directions file which has the combination of 2 fibers (this
%is from the directions inmage loaded in sim_signals
dirs_v1 = dirs.img(46,67,36,1:3); %d1=[0.0577,-0.6954,0.7163]; but it has to be a 1x1x3 matrix 
dirs_v2 = dirs.img(46,65,36,1:3); %d2=[0.9882,-0.1046,-0.1118]
new_dirs = zeros(length(afd1),length(mtr1)*length(mtr2),1,size(dirs.img,4)); %the time dimension has to be bigger than 3 or python crashes...
new_dirs(:,:,:,1:3)=repmat(dirs_v1,length(afd1),length(mtr1)*length(mtr2)); %this is the orientation of the 1st fiber
new_dirs(:,:,:,4:6)=repmat(dirs_v2,length(afd1),length(mtr1)*length(mtr2)); % this is the orientation of the 2nd fiber

%new icvf file (it's about the same for all 3 voxels, ~0.44)
new_icvf = 0.44*ones(length(afd1),length(mtr1)*length(mtr2),2); %also make 2 in time of else python complains

save_nii(make_nii(new_icvf),'sim_equation/new-icvf.nii');
save_nii(make_nii(new_dirs),'sim_equation/new-directions.nii');
save_nii(make_nii(afd),'sim_equation/new-afd.nii');
save_nii(make_nii(mtr),'sim_equation/new-mtr.nii');

%% now plot the results

%read in output
%no noise
mtr_out=load_nii('sim_equation/sim-equation-out-no-noise/mtr_.nii'); % get raw data

figure
for a=1:length(afd1)
    randc = rand(1,3);
    plot(mtr(a,:,1),mtr_out.img(a,:,1),'o','color',randc)
    text(0.3,0.4-0.01*a,strcat('AFD1=',num2str(afd(a,1,1)),' AFD2=',num2str(afd(a,1,2))),'color',randc)
    hold on

end
title ('MTR1: No noise case')
xlabel('Real')
ylabel('Estimated')
plot(mtr(a,:,1),mtr(a,:,1),'k') %unity line

figure
for a=1:length(afd1)
    randc = rand(1,3);
    plot(mtr(a,:,2),mtr_out.img(a,:,2),'o','color',randc)
    text(0.3,0.4-0.01*a,strcat('AFD1=',num2str(afd(a,1,1)),' AFD2=',num2str(afd(a,1,2))),'color',randc)
    hold on
    title ('MTR2: No noise case')
end
title ('MTR2: No noise case')
xlabel('Real')
ylabel('Estimated')
plot(mtr(a,:,1),mtr(a,:,1),'k') %unity line

%3% noise
mtr_out=load_nii('sim_equation/sim-equation-out/mtr_.nii'); % get raw data
figure
for a=1:length(afd1)
    randc = rand(1,3);
    plot(mtr(a,:,1),mtr_out.img(a,:,1),'o','color',randc)
    text(0.3,0.6-0.01*a,strcat('AFD1=',num2str(afd(a,1,1)),' AFD2=',num2str(afd(a,1,2))),'color',randc)
    hold on
    
end
title ('MTR1 3% noise case')
xlabel('Real')
ylabel('Estimated')
plot(mtr(a,:,1),mtr(a,:,1),'k')
ylim([0 0.65])

figure
for a=1:length(afd1)
    randc = rand(1,3);
    plot(mtr(a,:,2),mtr_out.img(a,:,2),'o','color',randc)
    text(0.3,0.6-0.01*a,strcat('AFD1=',num2str(afd(a,1,1)),' AFD2=',num2str(afd(a,1,2))),'color',randc)
    hold on
    
end
title ('MTR2 3% noise case')
xlabel('Real')
ylabel('Estimated')
plot(mtr(a,:,1),mtr(a,:,1),'k')
ylim([0 0.65])

%% now with many iterations
%this volume has each iteration in the z dimension
mtr_iter = load_nii('sim_equation/sim-equation-iters/mtr-merge-100.nii');
% use simulate_MTRs.py then average on the command line as per notes
mtr_avg=load_untouch_nii('sim_equation/sim-equation-iters/mtr-merge-100-avg_.nii');
mtr_sd=load_untouch_nii('sim_equation/sim-equation-iters/mtr-merge-100-sd_.nii');

figure
for a=1:length(afd1)
    randc = rand(1,3);
    plot(mtr(a,:,1),mtr_avg.img(a,:,1),'o','color',randc,'MarkerSize',10)
    text(0.3,0.6-0.01*a,strcat('AFD1=',num2str(afd(a,1,1)),' AFD2=',num2str(afd(a,1,2))),'color',randc)
    hold on
    errorbar(mtr(a,:,1),mtr_avg.img(a,:,1),mtr_sd.img(a,:,1),'LineStyle','none')
    
end
title ('MTR1-average 3% noise case')
xlabel('Real')
ylabel('Estimated')
plot(mtr(a,:,1),mtr(a,:,1),'k')
ylim([0 0.65])

%% plot the bias
%MTR1
figure
for a=1:length(afd1)
    randc = rand(1,3);
    plot(mtr(a,:,1),(mtr_avg.img(a,:,1)-mtr(a,:,1)),'o','color',randc,'MarkerSize',10)
    text(0.26,0.2-0.01*a,strcat('AFD1=',num2str(afd(a,1,1)),' AFD2=',num2str(afd(a,1,2))),'color',randc,'Fontsize',18)
    hold on
    errorbar(mtr(a,:,1),mtr_avg.img(a,:,1)-mtr(a,:,1),mtr_sd.img(a,:,1),'LineStyle','none','color',randc)
    %now write average for each AFD combo at each mtr
    for t=0:5:20
        avg_mtr=mean(mean(mtr_avg.img(a,t+1:t+5,1)-mtr(a,t+1:t+5,1)));
        std_mtr=mean(mean(mtr_sd.img(a,t+1:t+5,1)));
        plot(mean(mean(mtr(:,t+1:t+5,1))),avg_mtr,'x','MarkerSize',25)
        text(mean(mean(mtr(:,t+1:t+5,1)))+0.005,0.05-0.01*a,strcat(sprintf('%2.1E',avg_mtr),' std=',sprintf('%2.1E',std_mtr)),...
            'color',randc,'Fontsize',15); %(this will display the bias in MTR1 for all combinations of MTR2 for each afd combo
    end
end
title ('MTR1 BIAS 3% noise case, 100 iterations')
xlabel('Real')
ylabel('Estimated')
plot(mtr(a,:,1),mtr(a,:,1)-mtr(a,:,1),'k')
ylim([-0.2 0.2])

% now add the average over all AFD combos for each MTR1 (this will display the bias in MTR1 for all combinations of MTR2 and AFDs)
for t=0:5:20
    avg_mtr=mean(mean(mtr_avg.img(:,t+1:t+5,1)-mtr(:,t+1:t+5,1)));
    plot(mean(mean(mtr(:,t+1:t+5,1))),avg_mtr,'*','MarkerSize',35)
    text(mean(mean(mtr(:,t+1:t+5,1))),-0.15,strcat('avg over all:',num2str(avg_mtr)),'Fontsize',16);%(this will display the bias in MTR1 for all combinations of MTR2 and AFDs)
end

%MTR2
figure
for a=1:length(afd1)
    randc = rand(1,3);
    plot(mtr(a,:,2),(mtr_avg.img(a,:,2)-mtr(a,:,2)),'o','color',randc,'MarkerSize',10)
    text(0.26,0.2-0.01*a,strcat('AFD1=',num2str(afd(a,1,1)),' AFD2=',num2str(afd(a,1,2))),'color',randc)
    hold on
    errorbar(mtr(a,:,2),mtr_avg.img(a,:,2)-mtr(a,:,2),mtr_sd.img(a,:,2),'LineStyle','none','color',randc)

end
title ('MTR2 BIAS 3% noise case, 100 iterations')
xlabel('Real')
ylabel('Estimated')
plot(mtr(a,:,2),mtr(a,:,2)-mtr(a,:,2),'k')
ylim([-0.2 0.2])

%% 1% noise
%% now with many iterations
%this volume has each iteration in the z dimension
mtr_iter = load_nii('sim_equation/sim-equation-iters-1pc/mtr-merge-100.nii');
% use simulate_MTRs.py then average on the command line as per notes
mtr_avg=load_untouch_nii('sim_equation/sim-equation-iters-1pc/mtr-merge-100-avg_.nii');
mtr_sd=load_untouch_nii('sim_equation/sim-equation-iters-1pc/mtr-merge-100-sd_.nii');

figure
for a=1:length(afd1)
    randc = rand(1,3);
    plot(mtr(a,:,1),mtr_avg.img(a,:,1),'o','color',randc,'MarkerSize',10)
    text(0.3,0.6-0.01*a,strcat('AFD1=',num2str(afd(a,1,1)),' AFD2=',num2str(afd(a,1,2))),'color',randc)
    hold on
    errorbar(mtr(a,:,1),mtr_avg.img(a,:,1),mtr_sd.img(a,:,1),'LineStyle','none')
    
end
title ('MTR1-average 1% noise case')
xlabel('Real')
ylabel('Estimated')
plot(mtr(a,:,1),mtr(a,:,1),'k')
ylim([0 0.65])

%% plot the bias
%MTR1
figure
%for a=1:length(afd1)
for a =1:6 %only certain combinations of AFD
    randc = rand(1,3);
    plot(mtr(a,:,1),(mtr_avg.img(a,:,1)-mtr(a,:,1)),'o','color',randc,'MarkerSize',10)
    text(0.26,0.08-0.005*a,strcat('AFD1=',num2str(afd(a,1,1)),' AFD2=',num2str(afd(a,1,2))),'color',randc,'Fontsize',18)
    hold on
    errorbar(mtr(a,:,1),mtr_avg.img(a,:,1)-mtr(a,:,1),mtr_sd.img(a,:,1),'LineStyle','none','color',randc)
    %now write average for each AFD combo at each mtr
    for t=0:5:20
        avg_mtr=mean(mean(mtr_avg.img(a,t+1:t+5,1)-mtr(a,t+1:t+5,1)));
        std_mtr=mean(mean(mtr_sd.img(a,t+1:t+5,1)));
        plot(mean(mean(mtr(:,t+1:t+5,1))),avg_mtr,'x','MarkerSize',25)
        text(mean(mean(mtr(:,t+1:t+5,1)))+0.005,0.05-0.01*a,strcat(sprintf('%2.1E',avg_mtr),' std=',sprintf('%2.1E',std_mtr)),...
            'color',randc,'Fontsize',12); %(this will display the bias in MTR1 for all combinations of MTR2 for each afd combo
    end
end
title ('MTR1 BIAS 1% noise case, 100 iterations')
xlabel('Real')
ylabel('Estimated')
plot(mtr(a,:,1),mtr(a,:,1)-mtr(a,:,1),'k')
ylim([-0.07, 0.07])

% now add the average over all AFD combos for each MTR1 (this will display the bias in MTR1 for all combinations of MTR2 and AFDs)
for t=0:5:20
    avg_mtr=mean(mean(mtr_avg.img(:,t+1:t+5,1)-mtr(:,t+1:t+5,1)));
    plot(mean(mean(mtr(:,t+1:t+5,1))),avg_mtr,'*','MarkerSize',35)
    text(mean(mean(mtr(:,t+1:t+5,1))),-0.05,strcat('avg over all:',num2str(avg_mtr)),'Fontsize',12);%(this will display the bias in MTR1 for all combinations of MTR2 and AFDs)
end

%now plot difference in MTR1 and MTR2, i.e. withi this amount of noise, can
%we tell them apart. plot true difference vs what was found

figure
%for a=1:length(afd1)
for a =1:6 %only certain combinations of AFD
    randc = rand(1,3);
    plot(mtr(a,:,1)-mtr(a,:,2),(mtr_avg.img(a,:,1)-mtr_avg.img(a,:,2)),'o','color',randc,'MarkerSize',10) %mtr1-mtr2
    text(0.26,0.08-0.005*a,strcat('AFD1=',num2str(afd(a,1,1)),' AFD2=',num2str(afd(a,1,2))),'color',randc,'Fontsize',18)
    hold on
    errorbar(mtr(a,:,1)-mtr(a,:,2),mtr_avg.img(a,:,1)-mtr_avg.img(a,:,2),mtr_sd.img(a,:,1),'LineStyle','none','color',randc)
    %now write average for each AFD combo at each mtr
%     for t=0:5:20
%         avg_mtr=mean(mean(mtr_avg.img(a,t+1:t+5,1)-mtr(a,t+1:t+5,1)));
%         std_mtr=mean(mean(mtr_sd.img(a,t+1:t+5,1)));
%         plot(mean(mean(mtr(:,t+1:t+5,1))),avg_mtr,'x','MarkerSize',25)
%         text(mean(mean(mtr(:,t+1:t+5,1)))+0.005,0.05-0.01*a,strcat(sprintf('%2.1E',avg_mtr),' std=',sprintf('%2.1E',std_mtr)),...
%             'color',randc,'Fontsize',12); %(this will display the bias in MTR1 for all combinations of MTR2 for each afd combo
%     end
end
title ('MTR1-MTR2 difference 1% noise case, 100 iterations')
xlabel('Real difference MTR1-MTR2')
ylabel('Estimated MTR1-MTR2')
%plot(mtr(a,:,1),mtr(a,:,1)-mtr(a,:,1),'k')
ylim([-0.1, 0.1])



