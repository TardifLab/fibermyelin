
% Function for plotting the bar graphs for ISMRM 2019 CC and cingulum mean
% T1s
% October 31, 2019

function func_T1_bar_plot_human_groups(cc_single_T1_list,cc_double_T1_list,cing_single_T1_list,cing_double_T1_list)

cc_single_meanT1 = mean(cc_single_T1_list);
cc_single_sigmaT1 = std(cc_single_T1_list);

cc_double_meanT1 = mean(cc_double_T1_list);
cc_double_sigmaT1 = std(cc_double_T1_list);

cing_single_meanT1 = mean(cing_single_T1_list);
cing_single_sigmaT1 = std(cing_single_T1_list);

cing_double_meanT1 = mean(cing_double_T1_list);
cing_double_sigmaT1 = std(cing_double_T1_list);

nan1 = 0;
nan2 = 0;

sigmas = [cc_single_sigmaT1 cing_single_sigmaT1 nan1 cc_double_sigmaT1 cing_double_sigmaT1];

figure;

xpos = [1.5 2 2.8 3.6 4.1];

bar1 = bar(xpos(1),cc_single_meanT1,0.4);
hold on
bar2 = bar(xpos(2),cing_single_meanT1,0.4);
hold on
bar3 = bar(xpos(3),0.4);
hold on
bar4 = bar(xpos(4),cc_double_meanT1,0.4);
hold on
bar5 = bar(xpos(5),cing_double_meanT1,0.4);

set(bar1,'FaceColor','black');
set(bar2,'FaceColor','red');
set(bar4,'FaceColor','black');
set(bar5,'FaceColor','red');

errorbar(xpos(1),cc_single_meanT1,sigmas(1),'k.');
errorbar(xpos(2),cing_single_meanT1,sigmas(2),'k.');
errorbar(xpos(4),cc_double_meanT1,sigmas(4),'k.');
errorbar(xpos(5),cing_double_meanT1,sigmas(5),'k.');

legend(sprintf('C.C. single fibre voxels: T1 = %g ± %g ms\n',round(cc_single_meanT1,0),round(cc_single_sigmaT1,0)), sprintf('Cing. single fibre voxels: T1 = %g ± %g ms\n',round(cing_single_meanT1,0),round(cing_single_sigmaT1,0)), sprintf('C.C. crossing fibre voxels: T1 = %g ± %g ms',round(cc_double_meanT1,0),round(cc_double_sigmaT1,0)), sprintf('Cing. crossing fibre voxels: T1 = %g ± %g ms\n',round(cing_double_meanT1,0),round(cing_double_sigmaT1,0))) % set the legend, title, and axis labels 

%legend(sprintf('C.C.: Single fiber voxels T_1 = %g ± %g ms\n Crossing fibre voxels T_1 = %g ± %g ms',round(cc_single_meanT1,0),round(cc_single_sigmaT1,0),round(cc_double_meanT1,0),round(cc_double_sigmaT1,0)), sprintf('Cing. single fibre voxels: T_1 = %g ± %g ms\n Crossing fibre voxels T_1 = %g ± %g ms',round(cing_single_meanT1,0),round(cing_single_sigmaT1,0)),round(cing_double_meanT1,0),round(cing_double_sigmaT1,0), sprintf('C.C. crossing fibre voxels: T1 = %g +/- %g ms',round(cc_double_meanT1,0),round(cc_double_sigmaT1,0)), sprintf('Cing. crossing fibre voxels: T1 = %g +/- %g ms\n',round(cing_double_meanT1,0),round(cing_double_sigmaT1,0))) % set the legend, title, and axis labels 

%legend(sprintf('C.C.: Single fiber voxels T_1 = %g ± %g ms\n Crossing fibre voxels T_1 = %g ± %g ms',round(cc_single_meanT1,0),round(cc_single_sigmaT1,0),round(cc_double_meanT1,0),round(cc_double_sigmaT1,0)), sprintf('Cing. single fibre voxels: T_1 = %g ± %g ms\n Crossing fibre voxels T_1 = %g ± %g ms',round(cing_single_meanT1,0),round(cing_single_sigmaT1,0)),round(cing_double_meanT1,0),round(cing_double_sigmaT1,0)) % set the legend, title, and axis labels 


title('T1 in single and crossing fibre voxels');
ylabel('T1 (ms)');
xticks(xpos);
xticklabels({'pons','cst','','pons','cst'})
xtickangle(0);
ylim([0 1300]);
grid on;
grid minor;

end
