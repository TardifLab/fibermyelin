function box_plot_T1(sorted_T1,str1,str2,fiber,out)

% Box plot version of the function for plotting the T1 bar graphs for ISMRM 2019, MRM paper, and MEng thesis 
% 
% sorted_T1 : text file output from fibermyelin, with list of values for
%       each population and corssing
% str1 : string to search for to determine population (e.g. '~x', '~y', or '~z')
% str2 : string to search for to determine population (e.g. '~x', '~y', or '~z')
% fiber: string either 'cc' or 'cst' 
% out: output dir for figure
%%
% November 13, 2019
% Use for the IRDWI paper human experiments
% 
% June 2020 modified by Ilana to be a bit easier to use
%
% e.g. box_plot_T1('fixel-dir-output-CST-PONS-Dpar1000/sorted_T1s.txt','~x','~z','cst')

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

%get rid of 0s and 4000s
%pop1 = nonzeros(pop1);pop1 = pop1(pop1<4000);
%pop2 = nonzeros(pop2);pop1 = pop1(pop1<4000);
%cross1 = nonzeros(cross(:,1)); cross1=cross1(cross1<4000);
%cross2 = nonzeros(cross(:,2)); cross2=cross2(cross2<4000);
pop1 = rmoutliers(pop1);
pop2 = rmoutliers(pop2);
cross1 = rmoutliers(cross(:,1));
cross2 = rmoutliers(cross(:,2));

pop1_meanT1 = mean(pop1); pop1_sigmaT1 = std(pop1); %calculate the mean T1s and the standard deviation of the T1s
pop1_crossings_meanT1 = mean(cross1); pop1_crossings_sigmaT1 = std(cross1);
pop2_meanT1 = mean(pop2); pop2_sigmaT1 = std(pop2);
pop2_crossings_meanT1 = mean(cross2); pop2_crossings_sigmaT1 = std(cross2);

sigmas = [round(pop1_sigmaT1,0) round(pop2_sigmaT1,0) round(pop1_crossings_sigmaT1,0) round(pop2_crossings_sigmaT1,0)]; % sigma list so the bars are put at the correct position on the plot

xpos = [1.5 1.6 1.8 1.9]; % set the position of the boxes so that like ones are grouped together 
%xpos = [1.5 1.6 1.75 1.85]; % set the position of the boxes so that like ones are grouped together 
    % redo with these gaps

figure; % Make the figure 

box1 = pop1'; % Make and fill the boxes 
box2 = pop2';
box3 = cross1;
box4 = cross2;

boxes = [box1; box2; box3; box4]; % Store the boxes
g = [zeros(length(box1),1); ones(length(box2),1); 2*ones(length(box3),1); 3*ones(length(box4),1)]; % Place the content of the boxes

boxplot(boxes,g,'positions',xpos,'Symbol','k.'); % Display the boxes

set(gcf,'color','w'); % set the figure background colour to white 
set(gca,'FontSize',8.5); 
set(gca,'FontName','Helvetica');
%set(gca,'TickLength',[0 0]);

%text(1.542,550,'*','FontSize',15); % add the star for single fibers in the plot
%text(1.842,550,'*','FontSize',15); % add the star for crossing fibers in the plot


% Set the colors of the boxes
color = ['w', 'k', 'w', 'k'];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',0.3); % fills in the boxes with color; FaceAlpha means opcaity
end
c = get(gca, 'Children'); % Would be used to match colours with legend. Done manually instead. 

if strcmp(fiber,'cc')
    % CC-CING ROI - Custom legend with colours matching box colours
    [~,h_legend] = legend(sprintf(' Genu of the corpus callosum \n Single fiber mean T_1 = %g ± %g ms\n Crossings mean T_1 = %g ± %g ms',...
       round(pop1_meanT1,0),round(pop1_sigmaT1,0),round(pop1_crossings_meanT1,0),round(pop1_crossings_sigmaT1,0)),...
       sprintf(' Cingulate bundles \n Single fiber mean T_1 = %g ± %g ms\n Crossings mean T_1 = %g ± %g ms',...
       round(pop2_meanT1,0),round(pop2_sigmaT1,0),round(pop2_crossings_meanT1,0),round(pop2_crossings_sigmaT1,0)),...
       'Location','eastoutside');
else
    % PONS-CST ROI - Custom legend with colours matching box colours
    [~,h_legend] = legend(sprintf(' Pons \n Single fiber mean T_1 = %g ± %g ms\n Crossings mean T_1 = %g ± %g ms',...
        round(pop1_meanT1,0),round(pop1_sigmaT1,0),round(pop1_crossings_meanT1,0),round(pop1_crossings_sigmaT1,0)),...
        sprintf(' Corticospinal tract \n Single fiber mean T_1 = %g ± %g ms\n Crossings mean T_1 = %g ± %g ms',...
        round(pop2_meanT1,0),round(pop2_sigmaT1,0),round(pop2_crossings_meanT1,0),round(pop2_crossings_sigmaT1,0)),...
        'Location','eastoutside');
end
    
PatchInLegend = findobj(h_legend,'type','patch');
set(PatchInLegend(1), 'FaceAlpha', 0.3,'FaceColor',[0 0 0]); % Set the opacity equal to the box opacity
set(PatchInLegend(2), 'FaceAlpha', 0.3,'FaceColor',[1 1 1]);
PatchInLegend(1).Vertices = [0.05 0.6; 0.15 0.6; 0.15 0.83; 0.05 0.83]; % Set the position of the legend boxes 
PatchInLegend(2).Vertices = [0.05 -0.05; 0.15 -0.05; 0.15 0.18; 0.05 0.18];

TextInLegend = findobj(h_legend,'type','text');
set(TextInLegend,'FontSize',8.5);
TextInLegend(1).Position = [0.2052 0.65 0]; % Set the position of the legend text (space it out)
TextInLegend(2).Position = [0.2052 0.01 0];

%text(2.003,1227,'*','FontSize',15); % add the star to the legend 
%text(2.05,1240,'p < 0.001','FontSize',8.5); % add the star to the legend
%text(2.01,1227,'*','FontSize',15); % add the star to the legend 
%text(2.07,1240,'p < 0.001','FontSize',8.5); % add the star to the legend


legend boxoff;

ylabel('T_1 (ms)');
ylim([500 1500]); % Set the y limit
set(gca,'xtick',[mean(xpos(1:2)) mean(xpos(3:4))]); % Define the xticks and labels
set(gca,'xticklabel',{'Single fibers','Crossings'});

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 5]; % Rectangular figure size for MRM
%fig.PaperPosition = [0 0 3.42 3]; % square figure size for MRM

%print('T1_BOX_PONS-CST-ROI_6p9x4p_simplified','-dpng','-r0') % Save the figure as a .png image
thresh = 0.05;

[h_singles,p_singles] = ttest2(pop1,pop2);  % do the t-test // If the crossings are the same as the singles - is good!
[h_crossings,p_crossings] = ttest2(cross1,cross2);

disp(' ');
disp('Singles:');
disp(h_singles); % If h =1, statistically independent (i.e. significant difference). If 0, not independent
disp(p_singles);  % If p<0.05, stat sig
disp('Crossings:');
disp(h_crossings);
disp(p_crossings);

ystart = 1400;
tt = strcat('pop1 vs pop2 singles: ',num2str(p_singles));
text(2,ystart,tt,'FontSize',11); % add the star to the legend 
tt = strcat('pop1 vs pop2 crossings: ',num2str(p_crossings));
text(2,ystart-50,tt,'FontSize',11); % add the star to the legend 


% now compare singles and crossings of the same pop
disp('Population 1');
[h,p] = ttest2(pop1,cross1)  % do the t-test // If the crossings are the same as the singles - is good!
tt = strcat('pop1 (s vx x): ',num2str(p));
text(2,ystart-100,tt,'FontSize',11); 
disp('Population 2');
[h,p] = ttest2(pop2,cross2)
tt = strcat('pop2 (s vx x): ',num2str(p));
text(2,ystart-150,tt,'FontSize',11); 


name = strcat(out,'/T1_BOX-',fiber);
print(name,'-dpng','-r0') % Save the figure as a .png image 
end





