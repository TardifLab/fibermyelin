function sim_figs_nonoise(in,outdir,T1)

% in : input t1 file name
% outdir : output directory
% T1: array of T1 values [T1_1 T1_2]
SNR=inf;

%read in output
t1_out=load_nii(in); % get raw data
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
title (strcat('%T1 (bias) for Fiber 1, T1nom=',num2str(T1(1)),'ms; SNR=',num2str(SNR)))
xlabel('Tensor shape for fiber 1')
ylabel('Tensor shape for fiber 2')
print(strcat(outdir,'/T1-fiber1'),'-dpng','-r0')

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
title (strcat('%T1 (bias) for Fiber 2, T1nom=',num2str(T1(2)),'ms: SNR=',num2str(SNR)))
xlabel('Tensor shape for fiber 1')
ylabel('Tensor shape for fiber 2')
print(strcat(outdir,'/T1-fiber2'),'-dpng','-r0')

