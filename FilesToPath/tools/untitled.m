open('../AR_test.fig')
copyfig(gcf)
xlabel('')
ylabel('')
figure_to_publish('../AR_test2')


%% transient_time
open('../transient_time.fig')
copyfig(gcf)

figure_to_publish('../transient_time')

%% kuromoto
open('../kuromoto.fig')
copyfig(gcf)
xlabel('')
ylabel('')
figure_to_publish('../kuromoto2')

%% NS 
open('../NS_sigma.fig')
copyfig(gcf)
xlabel('')
ylabel('')
legend('off')
figure_to_publish('../NS_sigma3')
set(gca,'XScale','linear');

