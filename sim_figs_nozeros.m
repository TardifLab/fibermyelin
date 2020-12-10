function sim_figs_nozeros(in,outdir,T1,SNR)
%%
%% This is used to plot the results of simulations with several iterations with noise

% in_avg : input t1 file name (iterations are in z dimension)
% outdir : output directory
% T1: array of T1 nominalvalues [T1_1 T1_2]
% SNR: SNR level

%% now with many iterations
% use simulate_tensors.py then combine on the command line as per notes

% t1_avg=load_untouch_nii(in_avg);
% t1_sd=load_untouch_nii(in_std);
%nomT1_1 = T1(1).*ones(size(t1_avg.img(:,:,1)));
%nomT1_2 = T1(2).*ones(size(t1_avg.img(:,:,2)));

% average and std, remove 0s first
t1=load_untouch_nii(in);
% failure rate T1>3000ms amd T1=0ms
%fail=(nnz(t1.img==0)+nnz(t1.img>3000))/nnz(t1.img(:,:,:,1));

%remove excessively high values
%t1.img(t1.img>2000)=0;

%average only non-zero values
for i=1:size(t1.img,1)
    for j=1:size(t1.img,2)
%         t1_avg.img(i,j,1) = mean(nonzeros(t1.img(i,j,:,1)));
%         t1_sd.img(i,j,1) = std(nonzeros(t1.img(i,j,:,1)));
%         t1_avg.img(i,j,2) = mean(nonzeros(t1.img(i,j,:,2)));
%         t1_sd.img(i,j,2) = std(nonzeros(t1.img(i,j,:,2)));
%         t1_avg.img(i,j,1) = mean(rmoutliers(permute(t1.img(i,j,:,1),[3 2 1]),'percentiles',[10 90]));
%         t1_sd.img(i,j,1) = std(rmoutliers(permute(t1.img(i,j,:,1),[3 2 1]),'percentiles',[10 90]));
%         t1_avg.img(i,j,2) = mean(rmoutliers(permute(t1.img(i,j,:,2),[3 2 1]),'percentiles',[10 90]));
%         t1_sd.img(i,j,2) = std(rmoutliers(permute(t1.img(i,j,:,2),[3 2 1]),'percentiles',[10 90]));
        t1_avg.img(i,j,1) = mean(rmoutliers(permute(t1.img(i,j,:,1),[3 2 1]))); %need permute because rmoutliers only works on 1D arrays
        t1_sd.img(i,j,1) = std(rmoutliers(permute(t1.img(i,j,:,1),[3 2 1])));
        t1_avg.img(i,j,2) = mean(rmoutliers(permute(t1.img(i,j,:,2),[3 2 1])));
        t1_sd.img(i,j,2) = std(rmoutliers(permute(t1.img(i,j,:,2),[3 2 1])));
        [h,p.img(i,j)] = ttest2(rmoutliers(permute(t1.img(i,j,:,1),[3 2 1])),rmoutliers(permute(t1.img(i,j,:,2),[3 2 1])));
%         

    end
end

bias_t1_1 = (t1_avg.img(:,:,1)-T1(1))./T1(1)*100;
bias_t1_2 = (t1_avg.img(:,:,2)-T1(2))./T1(2)*100;
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
title (strcat('% T1 (bias) for Fiber 1, T1nom=',num2str(T1(1)),'ms; SNR = ',num2str(SNR)))
xlabel('Tensor shape for fiber 1')
ylabel('Tensor shape for fiber 2')
print(strcat(outdir,'/T1-fiber1'),'-dpng','-r0')
biasT1(:,:,1) = bias_t1_1;
biasT1(:,:,2) = bias_t1_2;
nii = make_nii(biasT1);
save_nii(nii, strcat(outdir,'/T1-biasFiber.nii'))



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
title (strcat('% T1 (bias) for Fiber 2, T1nom=',num2str(T1(2)),'ms; SNR = ',num2str(SNR)))
xlabel('Tensor shape for fiber 1')
ylabel('Tensor shape for fiber 2')
print(strcat(outdir,'/T1-fiber2'),'-dpng','-r0')

%% standard deviation
pcstd_t1_1 = t1_sd.img(:,:,1)./T1(1)*100; %in percent
pcstd_t1_2 = t1_sd.img(:,:,2)./T1(2)*100;
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
title (strcat('% T1 STD for Fiber 1, T1nom=',num2str(T1(1)),'ms; SNR = ',num2str(SNR)))
xlabel('Tensor shape for fiber 1')
ylabel('Tensor shape for fiber 2')
print(strcat(outdir,'/T1std-fiber1'),'-dpng','-r0')

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
title (strcat('% T1 STD for Fiber 2, T1nom=',num2str(T1(2)),'ms; SNR = ',num2str(SNR)))
xlabel('Tensor shape for fiber 1')
ylabel('Tensor shape for fiber 2')
%text(0.5,1,strcat('fit fails= ',num2str(fail*100),'%'))
print(strcat(outdir,'/T1std-fiber2'),'-dpng','-r0')
pcstdT1(:,:,1) = pcstd_t1_1;
pcstdT1(:,:,2) = pcstd_t1_2;
nii = make_nii(pcstdT1);
save_nii(nii, strcat(outdir,'/T1-pcstd.nii'))

%% Can we we tell the values apart? p value
figure
imagesc(flip(fliplr(permute(p.img,[2 1 ]))),[0 5E-2])
xticks(1:7)
xticklabels({'0.9','0.8','0.7','0.6','0.5','0.4','0.3'}) 
yticks([1:7])
yticklabels({'0.9','0.8','0.7','0.6','0.5','0.4','0.3'})
colormap(bluewhitered), colorbar
title (strcat('p-value, SNR = ',num2str(SNR)))
xlabel('Tensor shape for fiber 1')
ylabel('Tensor shape for fiber 2')
%text(0.5,1,strcat('fit fails= ',num2str(fail*100),'%'))
print(strcat(outdir,'/pval'),'-dpng','-r0')
nii = make_nii(p.img);
save_nii(nii, strcat(outdir,'/pval.nii'))



