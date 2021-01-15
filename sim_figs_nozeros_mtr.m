function sim_figs_nozeros_mtr(in,outdir,SNR)
%%
%% This is used to plot the results of simulations with several iterations with noise

% in_avg : input mtr file name (iterations are in z dimension)
% outdir : output directory
% SNR: SNR level

%the fake data is made in sim_equation.m, these are the combos
afd1=0:0.1:1;
afd2=0:0.1:1;
mtr1=0.2:0.05:0.4;
mtr2=0.2:0.05:0.4;

AFD = zeros(length(afd1),length(mtr1)*length(mtr2),1,2); %time dimension is 2 (number of fibers)
MTR = zeros(length(afd1),length(mtr1)*length(mtr2),1,2); %time dimension is 2 (number of fibers)

%afd array: each row is a different afd combination, each column is the
%same
for a=1:length(afd1)
    AFD(a,1:end,1,1) = afd1(a); 
    AFD(a,1:end,1,2) = afd2(end+1-a);
end

count=1;
for m1=1:length(mtr1)
    for m2=1:length(mtr2)
        MTR(1:end,count,1) = mtr1(m1);
        MTR(1:end,count,2) = mtr2(m2);
        count=count+1;
    end
end

%% now with many iterations
% use simulate_MTRs.py then combine on the command line as per notes
mtr=load_untouch_nii(in);


%average only non-zero values
for i=1:size(mtr.img,1)
    for j=1:size(mtr.img,2)
        mtr_avg.img(i,j,1) = mean(rmoutliers(permute(mtr.img(i,j,:,1),[3 2 1]))); %need permute because rmoutliers only works on 1D arrays
        mtr_sd.img(i,j,1) = std(rmoutliers(permute(mtr.img(i,j,:,1),[3 2 1])));
        mtr_avg.img(i,j,2) = mean(rmoutliers(permute(mtr.img(i,j,:,2),[3 2 1])));
        mtr_sd.img(i,j,2) = std(rmoutliers(permute(mtr.img(i,j,:,2),[3 2 1])));
        [h,p.img(i,j)] = ttest2(rmoutliers(permute(mtr.img(i,j,:,1),[3 2 1])),rmoutliers(permute(mtr.img(i,j,:,2),[3 2 1])));
    end
end

bias_mtr_1 = (mtr_avg.img(:,:,1)-MTR(:,:,:,1))./MTR(:,:,:,1)*100;
bias_mtr_2 = (mtr_avg.img(:,:,2)-MTR(:,:,:,2))./MTR(:,:,:,2)*100;

% fiber 1
figure
dat = flip(fliplr(permute(bias_mtr_1,[2 1 ])));
count=5;
for p = 1:5:25
    subplot(5,1,6-count)
    imagesc(dat(p:p+4,:),[-20 20])
    xticks(1:11)
    xticklabels({'1-0','0.9-0.1','0.8-0.2','0.7-0.3','0.6-0.4','0.5-0.5','0.4-0.6','0.7-0.3','0.8-0.2','0.9-0.1','1-0'})
    xtickangle(45)
    yticks(1:5)
    ylabel('MTR2')
    yticklabels({'0.4','0.35','0.3','0.25','0.2'})
    title(strcat('MTR1 = ',num2str(mtr1(count))))
    count=count-1;
    colormap(bluewhitered), colorbar
end
suptitle (strcat('% mtr (bias) for Fiber 1; SNR = ',num2str(SNR)))
xlabel('AFD1 - AFD2')

print(strcat(outdir,'/mtr-fiber1'),'-dpng','-r0')
biasmtr(:,:,1) = bias_mtr_1;
biasmtr(:,:,2) = bias_mtr_2;
nii = make_nii(biasmtr);
save_nii(nii, strcat(outdir,'/mtr-biasFiber.nii'))



%fiber 2
figure
dat = flip(fliplr(permute(bias_mtr_2,[2 1 ])));
count=5;
for p = 1:5:25
    subplot(5,1,6-count)
    imagesc(dat(p:p+4,:),[-20 20])
    xticks(1:11)
    xticklabels({'1-0','0.9-0.1','0.8-0.2','0.7-0.3','0.6-0.4','0.5-0.5','0.4-0.6','0.7-0.3','0.8-0.2','0.9-0.1','1-0'})
    xtickangle(45)
    yticks(1:5)
    ylabel('MTR2')
    yticklabels({'0.4','0.35','0.3','0.25','0.2'})
    title(strcat('MTR1 = ',num2str(mtr1(count))))
    count=count-1;
    colormap(bluewhitered), colorbar
end
suptitle (strcat('% mtr (bias) for Fiber 2; SNR = ',num2str(SNR)))
xlabel('AFD1 - AFD2')

print(strcat(outdir,'/mtr-fiber2'),'-dpng','-r0')

%% standard deviation
pcstd_mtr_1 = mtr_sd.img(:,:,1)./MTR(:,:,:,1)*100; %in percent
pcstd_mtr_2 = mtr_sd.img(:,:,2)./MTR(:,:,:,2)*100;
figure
dat=flip(fliplr(permute(pcstd_mtr_1,[2 1 ])));
count=5;
for p = 1:5:25
    subplot(5,1,6-count)
    imagesc(dat(p:p+4,:),[0 30])
    xticks(1:11)
    xticklabels({'1-0','0.9-0.1','0.8-0.2','0.7-0.3','0.6-0.4','0.5-0.5','0.4-0.6','0.7-0.3','0.8-0.2','0.9-0.1','1-0'})
    xtickangle(45)
    yticks(1:5)
    ylabel('MTR2')
    yticklabels({'0.4','0.35','0.3','0.25','0.2'})
    title(strcat('MTR1 = ',num2str(mtr1(count))))
    count=count-1;
    colormap(bluewhitered), colorbar
end
suptitle (strcat('% mtr STD for Fiber 1, SNR = ',num2str(SNR)))
xlabel('AFD1 - AFD2')
print(strcat(outdir,'/mtrstd-fiber1'),'-dpng','-r0')

%fiber 2
figure
dat=flip(fliplr(permute(pcstd_mtr_2,[2 1 ])));
count=5;
for p = 1:5:25
    subplot(5,1,6-count)
    imagesc(dat(p:p+4,:),[0 30])
    xticks(1:11)
    xticklabels({'1-0','0.9-0.1','0.8-0.2','0.7-0.3','0.6-0.4','0.5-0.5','0.4-0.6','0.7-0.3','0.8-0.2','0.9-0.1','1-0'})
    xtickangle(45)
    yticks(1:5)
    ylabel('MTR2')
    yticklabels({'0.4','0.35','0.3','0.25','0.2'})
    title(strcat('MTR1 = ',num2str(mtr1(count))))
    count=count-1;
    colormap(bluewhitered), colorbar
end

suptitle (strcat('% mtr STD for Fiber 2, SNR = ',num2str(SNR)))
xlabel('AFD1 - AFD2')

print(strcat(outdir,'/mtrstd-fiber2'),'-dpng','-r0')
pcstdmtr(:,:,1) = pcstd_mtr_1;
pcstdmtr(:,:,2) = pcstd_mtr_2;
nii = make_nii(pcstdmtr);
save_nii(nii, strcat(outdir,'/mtr-pcstd.nii'))

% %% Can we we tell the values apart? p value
% figure
% imagesc(flip(fliplr(permute(p.img,[2 1 ]))),[0 5E-2])
% xticks(1:7)
% xticklabels({'0.9','0.8','0.7','0.6','0.5','0.4','0.3'}) 
% yticks([1:7])
% yticklabels({'0.9','0.8','0.7','0.6','0.5','0.4','0.3'})
% colormap(bluewhitered), colorbar
% title (strcat('p-value, SNR = ',num2str(SNR)))
% xlabel('Tensor shape for fiber 1')
% ylabel('Tensor shape for fiber 2')
% %text(0.5,1,strcat('fit fails= ',num2str(fail*100),'%'))
% print(strcat(outdir,'/pval'),'-dpng','-r0')
% nii = make_nii(p.img);
% save_nii(nii, strcat(outdir,'/pval.nii'))



