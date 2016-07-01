names = {'alessandro_data3','anna_data3_a2','matias_data3_a2','pius_data3','thomas_data3','ziqiang_data3'};

h1 = figure;
h1.PaperType = 'a4';
h1.PaperUnits = 'centimeters';
h1.PaperPosition = [0, 0, 24, 29];

for i=1:6
    load(['C:\Users\linus\Google Drive\pitch shift experiment data\rec_vowel_exp_data\processed data\', names{i}, '.mat']);
    data3_const_shift.condition_name = 'const.';
    data3_var_shift.condition_name = 'var.';
    
    data3_const_shift.time_before_shift_ms = 50;
    data3_var_shift.time_before_shift_ms = 50;
    data3_const_shift.time_after_shift_ms = 350;
    data3_var_shift.time_after_shift_ms = 350;
    
    leg1=vowel_exp_plot_2c(data3_const_shift,data3_var_shift,subplot(6,2,i*2-1),subplot(6,2,i*2),'title',0);
    subplot(6,2,i*2-1);
    %ylim([-50 50]);
    
    m(i) = uicontrol('style','text');
    m(i).String = sprintf('%i:',i);
    m(i).FontSize = 16;
    m(i).BackgroundColor = 'white';
    m(i).FontWeight = 'bold';
    m(i).Units = 'normalized';
    m(i).Position = [0.02,1.01-i*0.142,0.05,0.02];
    %title(sprintf('Subject %i (reference: %i Hz)',i,data3_const_shift.piano_freq));
end
subplot(6,2,12);
xlabel('time [ms]');
subplot(6,2,11)
xlabel('time [ms]'); 
legend(leg1,'Orientation','horizontal','FontSize', 7,'Position',[0.37,-0.01,0.25,0.1]);%'Location', 'southoutside');

%subplot(6,2,10); ylim([0,5]);

print('Plot','-dpng', '-r300');