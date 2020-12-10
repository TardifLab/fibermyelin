%% Create synthetic data from the equation to run the simulations

% create a grid with different combinations of tensors


% now create an array of combinations of AD adn RD (axial and radial diffusivity)
% the average tensor for only b=1000 shell is:
%
% mrstats ad.mif -output mean -mask single.nii 
%AD = 0.00164679;
% std = 0.000227582
% mrstats rd.mif -output mean -mask single.nii 
%RD = 0.000363085;
% std = 0.0001014
% assuming normal distribution
% 
% Actually, let's read the values in from a textfile
%mrdump ad.mif AD.txt -mask single.nii
%mrdump rd.mif RD.txt -mask single.nii
%ADtxt=importdata('AD.txt');
%RDtxt=importdata('RD.txt');
%outdir = 'sim_equation/';

% let's use the whole brain to get a more representative FAs
%mrdump ad.mif AD-all.txt -mask ../all-singles.nii
%mrdump rd.mif RD-all.txt -mask ../all-singles.nii
%ADtxt=importdata('AD-all.txt');
%RDtxt=importdata('RD-all.txt');
%outdir = 'sim_equation-all/';

%ad=[prctile(ADtxt,10),prctile(ADtxt,25),prctile(ADtxt,50),prctile(ADtxt,75) prctile(ADtxt,90)];
%rd=[prctile(RDtxt,10),prctile(RDtxt,25),prctile(RDtxt,50),prctile(RDtxt,75) prctile(RDtxt,90)];

%let's make rounder FA values 0.3-0.9 
% ad = 0.0016*ones(7,1);
% rd = [9.8 8.1 6.7 5.4 4.1 2.8 1.4].*1E-4;
% outdir = 'sim_equation-diff-T1/';

%let's do this again, but keep MD at 7.5E-4
ad = [10.2 11 12.2 13.4 14.9 16.6 19].*1E-4;
rd = [6.2 5.7 5.2 4.6 3.8 3.0 1.8].*1E-4;
outdir = 'sim_equation-diff-T1-same-MD/';


AD = zeros(length(ad),length(rd),1,2); %time dimension is 2 (number of fibers)
RD = zeros(length(ad),length(rd),1,2); %time dimension is 2 (number of fibers)

% now make the array that we will input to the code
% there are length(ad)xlength(rd) combinations of tensors, with the info
% about each tensor being in a different time dimension
% make the combos first
for i=1:length(ad)
    tensor_combo(i,:) = [ad(i), rd(i)];
end

% now make all combinations of each tensor
count=1;
for m1=1:length(ad)
    for m2=1:length(ad)
        AD(m1,m2,1,1) = tensor_combo(m1,1);
        AD(m1,m2,1,2) = tensor_combo(m2,1);
        RD(m1,m2,1,1) = tensor_combo(m1,2);
        RD(m1,m2,1,2) = tensor_combo(m2,2);
        count=count+1;
    end
end

%create a new directions file which has the combination of 2 fibers (this
%is from the directions image loaded in sim_signals
dirs=load_nii('/home/bic/ilana/links/tardif/mt-diff/mt_diff_64_20191025_104945/directions_voxel_strides.nii');
dirs_v1 = dirs.img(46,67,36,1:3); %d1=[0.0577,-0.6954,0.7163]; but it has to be a 1x1x3 matrix, that's why it's tough to create 
dirs_v2 = dirs.img(46,65,36,1:3); %d2=[0.9882,-0.1046,-0.1118]
new_dirs = zeros(length(ad),length(rd),1,size(dirs.img,4)); %the time dimension has to be bigger than 3 or python crashes...
%new_dirs(:,:,:,1:3)=repmat(dirs_v1,length(ad),length(rd)); %this is the orientation of the 1st fiber
%new_dirs(:,:,:,4:6)=repmat(dirs_v2,length(ad),length(rd)); % this is the orientation of the 2nd fiber
dirs_v1(:,:,:,1:3)=[1 0 0];
dirs_v2(:,:,:,1:3)=[0 0 1];
new_dirs(:,:,:,1:3)=repmat(dirs_v1,length(ad),length(rd)); %this is the orientation of the 1st fiber
new_dirs(:,:,:,4:6)=repmat(dirs_v2,length(ad),length(rd)); % this is the orientation of the 2nd fiber


save_nii(make_nii(new_dirs),strcat(outdir,'/new-directions.nii'));
save_nii(make_nii(AD),strcat(outdir,'/new-AD.nii'));
save_nii(make_nii(RD),strcat(outdir,'/new-RD.nii'));

afd = 0.5*ones(length(ad),length(rd),1,2); %all 0.5 for now
save_nii(make_nii(afd),strcat(outdir,'/new-afd.nii'));
%80-20
afd(:,:,:,1) = 0.8*ones(length(ad),length(rd),1,1); 
afd(:,:,:,2) = 0.2*ones(length(ad),length(rd),1,1); 
save_nii(make_nii(afd),strcat(outdir,'/new-afd-p8p2.nii'));
%35-65
afd(:,:,:,1) = 0.35*ones(length(ad),length(rd),1,1); 
afd(:,:,:,2) = 0.65*ones(length(ad),length(rd),1,1); 
save_nii(make_nii(afd),strcat(outdir,'/new-afd-p35p65.nii'));

%what is the FA of these tensors
FA=sqrt((AD-RD).^2./(AD.^2+2.*RD.^2));
save_nii(make_nii(FA),strcat(outdir,'/new-FA.nii'));

%figure out the combinations of AD and RDfor FA = 0.4 to 0.8
lambda1=0.0016;
lambda2 = 1E-4:1E-5:11E-4;
fa = sqrt((lambda1-lambda2).^2./(lambda1.^2+2.*lambda2.^2));
figure
plot(lambda2,fa)
grid on
hold on
plot(lambda2,md.*100)
title('What lambda2 for specific FA?')

%now keep MD the same
MD  = 7.5E-4;
y1 = 7.5E-4:1E-5:20E-4;
y2 = (3*MD-y1)./2;
fa = sqrt((y1-y2).^2./(y1.^2+2.*y2.^2));
figure
plot3(y1,y2,fa)
xlabel('lamda1')
ylabel('lamda2')
zlabel('FA')
grid on
md = (y1+2.*y2)./3;
title('What  lambda2 lambda2 for MD=7.5E-4')


%% now plot the results

%read in output
%no noise
t1_out=load_nii(strcat(outdir,'/0pc/t1.nii')); % get raw data
figure
%hanky panky so it looks exactly like mrtrix
bias_t1 = t1_out.img-750;
imagesc(flip(fliplr(permute(bias_t1(:,:,1)./7.5,[2 1 ]))),[-20 20])
xticks([1:5])
xticklabels({'fatest','~fat','average','~skinny','skiniest'})
yticks([1:5])
yticklabels({'fatest','~fat','average','~skinny','skiniest'})
colormap(bluewhitered), colorbar
title ('%T1 (bias) for Fiber 1: No noise case')
xlabel('Tensor shape for fiber 1')
ylabel('Tensor shape for fiber 2')
print(strcat(outdir,'/0pc/T1-fiber1'),'-dpng','-r0')

%fiber 2 T1
figure
imagesc(flip(fliplr(permute(bias_t1(:,:,2)/7.5,[2 1 ]))),[-20 20])
xticks([1:5])
%xticklabels({'fatest','~fat','average','~skinny','skiniest'})
%xticklabels({'0.55','0.66','0.73','0.82','0.89'})
xticklabels({'0.7063','0.5799','0.4371','0.2951','0.1554'}) %whole brain
yticks([1:5])
yticklabels({'fatest','~fat','average','~skinny','skiniest'})
colormap(bluewhitered), colorbar
title ('%T1 (bias) for Fiber 2: No noise case')
xlabel('Tensor shape for fiber 1')
ylabel('Tensor shape for fiber 2')
print(strcat(outdir,'/0pc/T1-fiber2'),'-dpng','-r0')

%% if I hold So fixed at 500 (what I know it is)
%read in output
%no noise
t1_out=load_nii(strcat(outdir,'/0pc-So500/t1.nii')); % get raw data
T1 = [750 900];
nomT1_1 = T1(1).*ones(size(t1_out.img(:,:,1)));
nomT1_2 = T1(2).*ones(size(t1_out.img(:,:,2)));
figure
%hanky panky so it looks exactly like mrtrix
bias_t1_1 = (t1_out.img(:,:,1)-nomT1_1)./nomT1_1*100;
imagesc(flip(fliplr(permute(bias_t1_1,[2 1 ]))),[-20 20])
xticks([1:7])
%xticklabels({'fatest','~fat','average','~skinny','skiniest'})
xticklabels({'0.9','0.8','0.7','0.6','0.5','0.4','0.3'}) %whole brain
yticks([1:7])
%yticklabels({'fatest','~fat','average','~skinny','skiniest'})
yticklabels({'0.9','0.8','0.7','0.6','0.5','0.4','0.3'})
colormap(bluewhitered), colorbar
title ('%T1 (bias) for Fiber 1: No noise case [So fixed]')
xlabel('Tensor shape for fiber 1')
ylabel('Tensor shape for fiber 2')
print(strcat(outdir,'/0pc-So500/T1-fiber1'),'-dpng','-r0')

%fiber 2 T1
figure
bias_t1_2= (t1_out.img(:,:,2)-nomT1_2)./nomT1_2*100;
imagesc(flip(fliplr(permute(bias_t1_2,[2 1 ]))),[-20 20])
xticks([1:7])
%xticklabels({'fatest','~fat','average','~skinny','skiniest'})
%xticklabels({'0.55','0.66','0.73','0.82','0.89'})
%xticklabels({'0.1554','0.2951','0.4371','0.5799','0.7063'}) %whole brain
xticklabels({'0.9','0.8','0.7','0.6','0.5','0.4','0.3'}) 
yticks([1:7])
%yticklabels({'fatest','~fat','average','~skinny','skiniest'})
yticklabels({'0.9','0.8','0.7','0.6','0.5','0.4','0.3'})
colormap(bluewhitered), colorbar
title ('%T1 (bias) for Fiber 2: No noise case [So fixed]')
xlabel('Tensor shape for fiber 1')
ylabel('Tensor shape for fiber 2')
print(strcat(outdir,'/0pc-So500/T1-fiber2'),'-dpng','-r0')

%% now with many iterations
%this volume has each iteration in the z dimension
%t1_iter = load_nii('sim_equation/3pc/t1-merge-100.nii');
% use simulate_t1s.py then average on the command line as per notes
%t1_avg=load_untouch_nii(strcat(outdir,'/3pc/t1-merge-100-avg_.nii'));
%t1_sd=load_untouch_nii(strcat(outdir,'/3pc/t1-merge-100-sd_.nii'));
t1_avg=load_untouch_nii(strcat(outdir,'/3pc-So500/t1-merge-100-avg.nii'));
t1_sd=load_untouch_nii(strcat(outdir,'/3pc-So500/t1-merge-100-std.nii'));
% = t1_avg.img-750;
%bias_t1std = t1_sd.img-750;
T1 = [750 900];
nomT1_1 = T1(1).*ones(size(t1_avg.img(:,:,1)));
nomT1_2 = T1(2).*ones(size(t1_avg.img(:,:,2)));
bias_t1_1 = (t1_avg.img(:,:,1)-nomT1_1)./nomT1_1*100;
bias_t1_2 = (t1_avg.img(:,:,2)-nomT1_2)./nomT1_2*100;
figure
imagesc(flip(fliplr(permute(bias_t1_1,[2 1 ]))),[-20 20])
%xticks([1:5])
%xticklabels({'fatest','~fat','average','~skinny','skiniest'})
%yticks([1:5])
%yticklabels({'fatest','~fat','average','~skinny','skiniest'})
xticks(1:7)
xticklabels({'0.9','0.8','0.7','0.6','0.5','0.4','0.3'}) 
yticks([1:7])
yticklabels({'0.9','0.8','0.7','0.6','0.5','0.4','0.3'}) 
colormap(bluewhitered), colorbar
title ('% T1 (bias) for Fiber 1: 3% noise case')
xlabel('Tensor shape for fiber 1')
ylabel('Tensor shape for fiber 2')
print(strcat(outdir,'/3pc-So500/T1-fiber1'),'-dpng','-r0')

%fiber 2
figure
imagesc(flip(fliplr(permute(bias_t1_2,[2 1 ]))),[-20 20])
%xticks([1:5])
%xticklabels({'fatest','~fat','average','~skinny','skiniest'})
%yticks([1:5])
%yticklabels({'fatest','~fat','average','~skinny','skiniest'})
xticks(1:7)
xticklabels({'0.9','0.8','0.7','0.6','0.5','0.4','0.3'}) 
yticks([1:7])
yticklabels({'0.9','0.8','0.7','0.6','0.5','0.4','0.3'}) 
colormap(bluewhitered), colorbar
title ('% T1 (bias) for Fiber 2: 3% noise case')
xlabel('Tensor shape for fiber 1')
ylabel('Tensor shape for fiber 2')
print(strcat(outdir,'/3pc-So500/T1-fiber2'),'-dpng','-r0')

%% standard deviation
pcstd_t1_1 = t1_sd.img(:,:,1)./nomT1_1*100; %in percent
pcstd_t1_2 = t1_sd.img(:,:,2)./nomT1_2*100;
figure
imagesc(flip(fliplr(permute(pcstd_t1_1,[2 1 ]))),[0 30])
%xticks([1:5])
%xticklabels({'fatest','~fat','average','~skinny','skiniest'})
%yticks([1:5])
%yticklabels({'fatest','~fat','average','~skinny','skiniest'})
xticks(1:7)
xticklabels({'0.9','0.8','0.7','0.6','0.5','0.4','0.3'}) 
yticks([1:7])
yticklabels({'0.9','0.8','0.7','0.6','0.5','0.4','0.3'})
colormap(bluewhitered), colorbar
title ('% T1 STD for Fiber 1: 3% noise case')
xlabel('Tensor shape for fiber 1')
ylabel('Tensor shape for fiber 2')
print(strcat(outdir,'/3pc-So500/T1std-fiber1'),'-dpng','-r0')

%fiber 2
figure
imagesc(flip(fliplr(permute(pcstd_t1_2,[2 1 ]))),[0 30])
%xticks([1:5])
%xticklabels({'fatest','~fat','average','~skinny','skiniest'})
%yticks([1:5])
%yticklabels({'fatest','~fat','average','~skinny','skiniest'})
xticks(1:7)
xticklabels({'0.9','0.8','0.7','0.6','0.5','0.4','0.3'}) 
yticks([1:7])
yticklabels({'0.9','0.8','0.7','0.6','0.5','0.4','0.3'})
colormap(bluewhitered), colorbar
title ('% T1 STD for Fiber 2: 3% noise case')
xlabel('Tensor shape for fiber 1')
ylabel('Tensor shape for fiber 2')
print(strcat(outdir,'/3pc-So500/T1std-fiber2'),'-dpng','-r0')


%% 1.5% noise
t1_avg=load_untouch_nii(strcat(outdir,'/1p5-So500/t1-merge-100-avg.nii'));
t1_sd=load_untouch_nii(strcat(outdir,'/1p5-So500/t1-merge-100-std.nii'));
% = t1_avg.img-750;
%bias_t1std = t1_sd.img-750;
T1 = [750 900];
nomT1_1 = T1(1).*ones(size(t1_avg.img(:,:,1)));
nomT1_2 = T1(2).*ones(size(t1_avg.img(:,:,2)));
bias_t1_1 = (t1_avg.img(:,:,1)-nomT1_1)./nomT1_1*100;
bias_t1_2 = (t1_avg.img(:,:,2)-nomT1_2)./nomT1_2*100;
figure
imagesc(flip(fliplr(permute(bias_t1_1,[2 1 ]))),[-20 20])
%xticks([1:5])
%xticklabels({'fatest','~fat','average','~skinny','skiniest'})
%yticks([1:5])
%yticklabels({'fatest','~fat','average','~skinny','skiniest'})
xticks(1:7)
xticklabels({'0.9','0.8','0.7','0.6','0.5','0.4','0.3'}) 
yticks([1:7])
yticklabels({'0.9','0.8','0.7','0.6','0.5','0.4','0.3'}) 
colormap(bluewhitered), colorbar
title ('% T1 (bias) for Fiber 1: 1.5% noise case')
xlabel('Tensor shape for fiber 1')
ylabel('Tensor shape for fiber 2')
print(strcat(outdir,'/1p5-So500/T1-fiber1'),'-dpng','-r0')

%fiber 2
figure
imagesc(flip(fliplr(permute(bias_t1_2,[2 1 ]))),[-20 20])
%xticks([1:5])
%xticklabels({'fatest','~fat','average','~skinny','skiniest'})
%yticks([1:5])
%yticklabels({'fatest','~fat','average','~skinny','skiniest'})
xticks(1:7)
xticklabels({'0.9','0.8','0.7','0.6','0.5','0.4','0.3'}) 
yticks([1:7])
yticklabels({'0.9','0.8','0.7','0.6','0.5','0.4','0.3'}) 
colormap(bluewhitered), colorbar
title ('% T1 (bias) for Fiber 2: 1.5% noise case')
xlabel('Tensor shape for fiber 1')
ylabel('Tensor shape for fiber 2')
print(strcat(outdir,'/1p5-So500/T1-fiber2'),'-dpng','-r0')

%% standard deviation
pcstd_t1_1 = t1_sd.img(:,:,1)./nomT1_1*100; %in percent
pcstd_t1_2 = t1_sd.img(:,:,2)./nomT1_2*100;
figure
imagesc(flip(fliplr(permute(pcstd_t1_1,[2 1 ]))),[0 30])
%xticks([1:5])
%xticklabels({'fatest','~fat','average','~skinny','skiniest'})
%yticks([1:5])
%yticklabels({'fatest','~fat','average','~skinny','skiniest'})
xticks(1:7)
xticklabels({'0.9','0.8','0.7','0.6','0.5','0.4','0.3'}) 
yticks([1:7])
yticklabels({'0.9','0.8','0.7','0.6','0.5','0.4','0.3'})
colormap(bluewhitered), colorbar
title ('% T1 STD for Fiber 1: 1.5% noise case')
xlabel('Tensor shape for fiber 1')
ylabel('Tensor shape for fiber 2')
print(strcat(outdir,'/1p5-So500/T1std-fiber1'),'-dpng','-r0')

%fiber 2
figure
imagesc(flip(fliplr(permute(pcstd_t1_2,[2 1 ]))),[0 30])
%xticks([1:5])
%xticklabels({'fatest','~fat','average','~skinny','skiniest'})
%yticks([1:5])
%yticklabels({'fatest','~fat','average','~skinny','skiniest'})
xticks(1:7)
xticklabels({'0.9','0.8','0.7','0.6','0.5','0.4','0.3'}) 
yticks([1:7])
yticklabels({'0.9','0.8','0.7','0.6','0.5','0.4','0.3'})
colormap(bluewhitered), colorbar
title ('% T1 STD for Fiber 2: 1.5% noise case')
xlabel('Tensor shape for fiber 1')
ylabel('Tensor shape for fiber 2')
print(strcat(outdir,'/1p5-So500/T1std-fiber2'),'-dpng','-r0')


%% realistic SNR
meanS=5433.47;
stdNoise = 112.772;
SNR = meanS/stdNoise;
%Generate values from a normal distribution with meanS and standard deviation stdNoise.
S = meanS + stdNoise.*randn(200,1);
Re = meanS/sqrt(2) + stdNoise.*randn(200,1);
Im = meanS/sqrt(2) + stdNoise.*randn(200,1);
Scomb = sqrt(Re.^2 + Im.^2);
% do they end up the same-ish?
figure
histogram(S)
hold on
histogram(Scomb,'FaceColor',[1 0 0])
legend({'signal','combined'})
title('SNR of magnitude and Re and Im')



%% plot the bias
%t11
figure
for a=1:length(afd1)
    randc = rand(1,3);
    plot(t1(a,:,1),(t1_avg.img(a,:,1)-t1(a,:,1)),'o','color',randc,'MarkerSize',10)
    text(0.26,0.2-0.01*a,strcat('AFD1=',num2str(afd(a,1,1)),' AFD2=',num2str(afd(a,1,2))),'color',randc,'Fontsize',18)
    hold on
    errorbar(t1(a,:,1),t1_avg.img(a,:,1)-t1(a,:,1),t1_sd.img(a,:,1),'LineStyle','none','color',randc)
    %now write average for each AFD combo at each t1
    for t=0:5:20
        avg_t1=mean(mean(t1_avg.img(a,t+1:t+5,1)-t1(a,t+1:t+5,1)));
        std_t1=mean(mean(t1_sd.img(a,t+1:t+5,1)));
        plot(mean(mean(t1(:,t+1:t+5,1))),avg_t1,'x','MarkerSize',25)
        text(mean(mean(t1(:,t+1:t+5,1)))+0.005,0.05-0.01*a,strcat(sprintf('%2.1E',avg_t1),' std=',sprintf('%2.1E',std_t1)),...
            'color',randc,'Fontsize',15); %(this will display the bias in t11 for all combinations of t12 for each afd combo
    end
end
title ('t11 BIAS 3% noise case, 100 iterations')
xlabel('Real')
ylabel('Estimated')
plot(t1(a,:,1),t1(a,:,1)-t1(a,:,1),'k')
ylim([-0.2 0.2])

% now add the average over all AFD combos for each t11 (this will display the bias in t11 for all combinations of t12 and AFDs)
for t=0:5:20
    avg_t1=mean(mean(t1_avg.img(:,t+1:t+5,1)-t1(:,t+1:t+5,1)));
    plot(mean(mean(t1(:,t+1:t+5,1))),avg_t1,'*','MarkerSize',35)
    text(mean(mean(t1(:,t+1:t+5,1))),-0.15,strcat('avg over all:',num2str(avg_t1)),'Fontsize',16);%(this will display the bias in t11 for all combinations of t12 and AFDs)
end

%t12
figure
for a=1:length(afd1)
    randc = rand(1,3);
    plot(t1(a,:,2),(t1_avg.img(a,:,2)-t1(a,:,2)),'o','color',randc,'MarkerSize',10)
    text(0.26,0.2-0.01*a,strcat('AFD1=',num2str(afd(a,1,1)),' AFD2=',num2str(afd(a,1,2))),'color',randc)
    hold on
    errorbar(t1(a,:,2),t1_avg.img(a,:,2)-t1(a,:,2),t1_sd.img(a,:,2),'LineStyle','none','color',randc)

end
title ('t12 BIAS 3% noise case, 100 iterations')
xlabel('Real')
ylabel('Estimated')
plot(t1(a,:,2),t1(a,:,2)-t1(a,:,2),'k')
ylim([-0.2 0.2])

%% 1% noise
%% now with many iterations
%this volume has each iteration in the z dimension
t1_iter = load_nii('sim_equation/sim-equation-iters-1pc/t1-merge-100.nii');
% use simulate_t1s.py then average on the command line as per notes
t1_avg=load_untouch_nii('sim_equation/sim-equation-iters-1pc/t1-merge-100-avg_.nii');
t1_sd=load_untouch_nii('sim_equation/sim-equation-iters-1pc/t1-merge-100-sd_.nii');

figure
for a=1:length(afd1)
    randc = rand(1,3);
    plot(t1(a,:,1),t1_avg.img(a,:,1),'o','color',randc,'MarkerSize',10)
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



