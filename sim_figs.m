function sim_figs(in_avg,in_std,outdir,T1,SNR)
%%
%% This is used to plot the results of simulations with several iterations with noise

% in_avg : input t1 file name (iterations are in z dimension)
% outdir : output directory
% T1: array of T1 nominalvalues [T1_1 T1_2]
% SNR: SNR level

%% now with many iterations
% use simulate_tensors.py then combine on the command line as per notes

t1_avg=load_untouch_nii(in_avg);
t1_sd=load_untouch_nii(in_std);


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
print(strcat(outdir,'/T1std-fiber2'),'-dpng','-r0')