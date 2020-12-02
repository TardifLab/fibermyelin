function T1_bars(sorted_T1,str1,str2)

% Function for plotting the T1 bar graphs for ISMRM 2019, MRM paper, and MEng thesis 
%%
% sorted_T1 : text file output from fibermyelin, with list of values for
%       each population and corssing
% str1 : string to search for to determine population (e.g. '~x', '~y', or '~z')
% str2 : string to search for to determine population (e.g. '~x', '~y', or '~z')


%%
% November 13, 2019
% Use for the IRDWI paper human experiments
% 
% June 2020 modified by Ilana to be a bit easier to use
%
% e.g. T1_bars('fixel-dir-output-CST-PONS-Dpar1000/sorted_T1s.txt','~x','~z')

%%%%%%%%%
% this is not a sexy way of reading but I don't want to waste too much time
str=fileread(sorted_T1);
lines = regexp( str, '\n', 'split' ); % each line of the file in a cell

ind = find(contains(lines,str1)); %this should return 2 indeces, the 1st will be crossings, 2nd singles
indx = ind(1);
inds1 = ind(2); %index singles of 1st direction
ind2 = find(contains(lines,str2)); %singles of other direction
inds2 = ind2(2);

%read data for crossings
c = lines(1,indx+2:inds1-5); %-5 to remove all \n and words ... not sexy
pop1 = str2double(lines(1,inds1+2:inds2-2));
pop2 = str2double(lines(1,inds2+2:end-1));

% now split the crossings, there a smarter way of doing this...
for i=1:size(c,2)
    cross(i,:) = str2double(strsplit(c{1,i},' '));
end

%get rid of 0s
pop1 = nonzeros(pop1);
pop2 = nonzeros(pop2);

s
pop1_meanT1 = mean(pop1); pop1_sigmaT1 = std(pop1); %calculate the mean T1s and the standard deviation of the T1s
pop1_crossings_meanT1 = mean(cross(:,1)); pop1_crossings_sigmaT1 = std(cross(:,1));
pop2_meanT1 = mean(pop2); pop2_sigmaT1 = std(pop2);
pop2_crossings_meanT1 = mean(cross(:,2)); pop2_crossings_sigmaT1 = std(cross(:,2));

sigmas = [round(pop1_sigmaT1,0) round(pop2_sigmaT1,0) round(pop1_crossings_sigmaT1,0) round(pop2_crossings_sigmaT1,0)]; % sigma list so the bars are put at the correct position on the plot
xpos = [1.5 2 2.8 3.6 4.1 4.9 5.7 6.2]; % set the position of the bars so that like ones are grouped together 

figure;
bar1 = bar(xpos(1),round(pop1_meanT1,0),0.4); hold on % put the bars in specific positions on the plot so that like data are grouped 
bar2 = bar(xpos(2),round(pop2_meanT1,0),0.4); hold on
bar3 = bar(xpos(3), 0, 0.4); hold on
bar4 = bar(xpos(4),round(pop1_crossings_meanT1,0),0.4); hold on
bar5 = bar(xpos(5),round(pop2_crossings_meanT1,0),0.4); hold on

% First stat significance bar
y = 1200;
line([1.5,2],[y,y],'Color','k');
x11 = 1.5;
line([x11,x11],[1175,1200],'Color','k');
x12 = 2;
line([x12,x12],[1175,1200],'Color','k');
text(1.7,1230,'*'); % add the star 

% Second stat significance bar
line([3.6,4.1],[y,y],'Color','k');
x21 = 3.6;
line([x21,x21],[1175,1200],'Color','k');
x22 = 4.1;
line([x22,x22],[1175,1200],'Color','k');
text(3.8,1230,'*','FontSize',10);

h = zeros(3, 1);    % for making a custom legend (plots nothing)
h(1) = bar(NaN,NaN);
h(2) = bar(NaN,NaN);
h(3) = bar(NaN,NaN);
set(h(1),'FaceColor','black');
set(h(2),'FaceColor','white');

set(bar1,'FaceColor','black'); % set the visible bar colours
set(bar2,'FaceColor','white');
set(bar4,'FaceColor','black');
set(bar5,'FaceColor','white');

errorbar(xpos(1),pop1_meanT1,sigmas(1),'Color',[.5 .5 .5],'LineWidth',1);   % put the error / sigma bars with the corresponding correct bars 
errorbar(xpos(2),pop2_meanT1,sigmas(2),'Color',[.5 .5 .5],'LineWidth',1);
errorbar(xpos(4),pop1_crossings_meanT1,sigmas(3),'Color',[.5 .5 .5],'LineWidth',1);
errorbar(xpos(5),pop2_crossings_meanT1,sigmas(4),'Color',[.5 .5 .5],'LineWidth',1);

set(gcf,'color','w'); % set the figure background colour to white 
set(gca,'FontSize',8.5); % set the figure background colour to white 
set(gca,'FontName','Helvetica');
set(gca,'TickLength',[0 0]);

%lg = legend(h, sprintf(' Genu of CC\n Single fibre T_1 = %g ± %g ms\n Crossing fibre T_1 = %g ± %g ms',round(asp1_single_meanT1,0),round(asp1_single_sigmaT1,0),round(asp1_crossings_meanT1,0),round(asp1_crossings_sigmaT1,0)),sprintf(' Cingulum \n Single fibre T_1 = %g ± %g ms\n Crossing fibre T_1 = %g ± %g ms',round(asp2_single_meanT1,0),round(asp2_single_sigmaT1,0),round(asp2_crossings_meanT1,0),round(asp2_crossings_sigmaT1,0))); % custom legend
lg = legend(h, sprintf(' Corpus callosum\n Single fibre T_1 = %g ± %g ms\n Crossing fibre T_1 = %g ± %g ms',round(asp1_single_meanT1,0),round(asp1_single_sigmaT1,0),round(asp1_crossings_meanT1,0),round(asp1_crossings_sigmaT1,0)),sprintf(' Cingulate bundles \n Single fibre T_1 = %g ± %g ms\n Crossing fibre T_1 = %g ± %g ms',round(asp2_single_meanT1,0),round(asp2_single_sigmaT1,0),round(asp2_crossings_meanT1,0),round(asp2_crossings_sigmaT1,0))); % custom legend
lg.Location = 'eastoutside';

%title('afdthresh = 0.2; T_1 of genu of CC and cingulum in single and crossing fibre voxels');
%title('afdthresh = 0.2; T_1 of pons and CST in single and crossing fibre voxels');

lg.FontSize = 8.5;
legend boxoff;
ylabel('T_1 (ms)');
xticks(xpos);

xticks([1.74 3.85 5.95])
%xticklabels({sprintf('Single fibres       '),sprintf('Crossing fibres'),sprintf('IGNORE THIS       ')})
%set(gca,'XTickLabel',{sprintf('Single fibres (IR-DWI)'),sprintf('Crossings (IR-DWI)'),sprintf('  Crossings     (IR)  ')});
set(gca,'XTickLabel',{sprintf('Single fibers'),sprintf('Crossings')});
%fix_xticklabels(gca,0.1,{'FontSize',8.5});
ylim([0 1400]);

% Add the p-value star to the legend
text(5.6,1040,'*','FontSize',10); % add the star to the legend 
text(6.05,1040,'p < 0.001','FontSize',8.5); % add the star to the legend

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6.9 3.4];
print('T1_bars_CC-CING-ROI_6p9x4p_simplified','-dpng','-r0')

end





