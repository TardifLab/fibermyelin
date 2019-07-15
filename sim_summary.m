%%

T1nom = 750; %nominal T1
gap = 110:10:280;
shuffle0=[];
shuffle1_x=[];
shuffle2_x=[];
%% singles
for i=gap
    for s=0:2 %shuffle
        f_shuffle = strcat('vals-singles_',num2str(s),'_gap',num2str(i),'.txt');
        tmp=importdata(f_shuffle);
        vals=nonzeros(tmp(4:5,:));
        bias = mean(vals-T1nom);
        bias_std = std(vals-T1nom);
        switch s
            case 0
                shuffle0=[shuffle0; bias bias_std];
            case 1
                shuffle1_x=[shuffle1_x; bias bias_std];
            case 2
                shuffle2_x=[shuffle2_x; bias bias_std];            
        end
        
    end
end

%% crossings
shuffle0_x=[];
shuffle1_x=[];
shuffle2_x=[];
for i=gap
    for s=0:2 %shuffle
        f_shuffle = strcat('vals-crossings_',num2str(s),'_gap',num2str(i),'.txt');
        tmp=importdata(f_shuffle);
        vals=nonzeros(tmp(4:5,:));
        bias = mean(vals-T1nom);
        bias_std = std(vals-T1nom);
        switch s
            case 0
                shuffle0_x=[shuffle0_x; bias bias_std];
            case 1
                shuffle1_x=[shuffle1_x; bias bias_std];
            case 2
                shuffle2_x=[shuffle2_x; bias bias_std];            
        end
        
    end
end
%%plot
figure
subplot(2,2,1)
plot(gap,shuffle0(:,1),'b')
hold on
plot(gap,shuffle1(:,1),'g')
plot(gap,shuffle2(:,1),'r')
legend('shuffle 0','shuffle 1','shuffle 2','Location','northwest')
title('Bias; single fibres','Fontsize',14)
xlabel('gap [ms]')
ylabel('bias [ms]')
xlim([min(gap) max(gap)])
ylim([0 25])

subplot(2,2,3)
plot(gap,shuffle0(:,2),'b')
hold on
plot(gap,shuffle1(:,2),'g')
plot(gap,shuffle2(:,2),'r')
legend('shuffle 0','shuffle 1','shuffle 2','Location','northwest')
title('Standard deviation; single fibres','Fontsize',14)
xlabel('gap [ms]')
ylabel('standard deviation [ms]')
xlim([min(gap) max(gap)])
ylim([10 50])

%crossings
subplot(2,2,2)
plot(gap,shuffle0_x(:,1),'b')
hold on
plot(gap,shuffle1_x(:,1),'g')
plot(gap,shuffle2_x(:,1),'r')
legend('shuffle 0','shuffle 1','shuffle 2','Location','northwest')
title('Bias; crossing fibres','Fontsize',14)
xlabel('gap [ms]')
ylabel('bias [ms]')
xlim([min(gap) max(gap)])
ylim([0 25])

subplot(2,2,4)
plot(gap,shuffle0_x(:,2),'b')
hold on
plot(gap,shuffle1_x(:,2),'g')
plot(gap,shuffle2_x(:,2),'r')
legend('shuffle 0','shuffle 1','shuffle 2','Location','northwest')
title('Standard deviation; crossing fibres','Fontsize',14)
xlabel('gap [ms]')
ylabel('standard deviation [ms]')
xlim([min(gap) max(gap)])    
ylim([10 50])