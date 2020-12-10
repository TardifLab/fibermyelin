function sim_figs_nonoise_mtr(result,in_mtr,outdir)

% result : mtr result from simulations
% in_mtr : nominal MTR values image (see sim_equation.m, it's many combos of MTR)
% outdir : output directory
SNR=inf;

%read in output of sim
mtr_out=load_untouch_nii(result); % get raw data
mtr_nom=load_untouch_nii(in_mtr); % get raw data
nomMTR_1 = mtr_nom.img(:,:,1).*ones(size(mtr_out.img(:,:,1)));
nomMTR_2 = mtr_nom.img(:,:,2).*ones(size(mtr_out.img(:,:,2)));
figure
%hanky panky so it looks exactly like mrtrix
bias_mtr_1 = (mtr_out.img(:,:,1)-nomMTR_1)./nomMTR_1*100;
imagesc(flip(fliplr(permute(bias_mtr_1,[2 1 ]))),[-20 20])
xticks(1:11)
xticklabels({'1-0','0.9-0.1','0.8-0.2','0.7-0.3','0.6-0.4','0.5-0.5','0.4-0.6','0.7-0.3','0.8-0.2','0.9-0.1','1-0'}) 
xtickangle(45)
yticks([1:25])
yticklabels({'0.4','0.35','0.3','0.25','0.2','0.4','0.35','0.3','0.25','0.2','0.4','0.35','0.3','0.25','0.2','0.4','0.35','0.3','0.25','0.2','0.4','0.35','0.3','0.25','0.2'}) 
colormap(bluewhitered), colorbar
title (strcat('% mtr (bias) for Fiber 1; SNR = ',num2str(SNR)))
xlabel('AFD1 - AFD2')
ylabel('MTR1 { MTR2')
print(strcat(outdir,'/T1-fiber1'),'-dpng','-r0')

%fiber 2 T1
figure
bias_mtr_2= (mtr_out.img(:,:,2)-nomMTR_2)./nomMTR_2*100;
imagesc(flip(fliplr(permute(bias_mtr_2,[2 1 ]))),[-20 20])
xticks(1:11)
xticklabels({'1-0','0.9-0.1','0.8-0.2','0.7-0.3','0.6-0.4','0.5-0.5','0.4-0.6','0.7-0.3','0.8-0.2','0.9-0.1','1-0'}) 
xtickangle(45)
yticks([1:25])
yticklabels({'0.4','0.35','0.3','0.25','0.2','0.4','0.35','0.3','0.25','0.2','0.4','0.35','0.3','0.25','0.2','0.4','0.35','0.3','0.25','0.2','0.4','0.35','0.3','0.25','0.2'}) 
colormap(bluewhitered), colorbar
title (strcat('% mtr (bias) for Fiber 2; SNR = ',num2str(SNR)))
xlabel('AFD1 - AFD2')
ylabel('MTR1 { MTR2')
print(strcat(outdir,'/T1-fiber2'),'-dpng','-r0')

