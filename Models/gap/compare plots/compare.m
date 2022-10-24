%load experimental data for manu

load('ExpDataManu.mat')
params = xb;
[expNC13, expNC14]=manu_model(params);

% x = 1:40;
% tiledlayout(2,2)
% nexttile
% plot(x,tgt_val(:,1)) %experimental
% hold on
% plot(x,expNC13(:,1)) %model
% % hold on
% % plot(x,NC13(:,1))
% legend('Expt','ISRES+','ISRES')
% title('Hb')
% hold off
% 
% x = 1:40;
% nexttile
% plot(x,tgt_val(:,2)) %experimental
% hold on
% plot(x,expNC13(:,2)) %model
% % hold on
% % plot(x,NC13(:,2))
% legend('Expt','ISRES+','isres')
% title('Kr')
% hold off
% 
% x = 1:40;
% nexttile
% plot(x,tgt_val(:,3)) %experimental
% hold on
% plot(x,expNC13(:,3)) %model
% % hold on
% % plot(x,NC13(:,3))
% legend('Expt','ISRES+','isres')
% title('Gt')
% hold off
% 
% x = 1:40;
% nexttile
% plot(x,expNC13(:,4)) %model
% % hold on
% % plot(x,NC13(:,4))
% hold off
% legend('ISRES+','isres')
% title('Kni')

%experimental
Hbb_E = NC14(:,8);
Kr_E = NC14(:,14);
Gt_E = NC14(:,22);
Kni_E = NC14(:,28);

%model
Hbb_M = expNC14(:,8);
Kr_M = expNC14(:,14);
Gt_M = expNC14(:,22);
Kni_M = expNC14(:,28);

figure;
t = tiledlayout(2,2);


x = 35:91;

nexttile
plot(x,Hbb_E,'Linewidth',2)
hold on
plot(x,Hbb_M,'Linewidth',2)
% hold on
% plot(x,Hbb)
legend('Data','ISRES+','Location','best')
title('\it hunchback')
set(gca,'FontSize',14)
hold off

nexttile
plot(x,Kr_E,'Linewidth',2)
hold on
plot(x,Kr_M,'Linewidth',2)
% hold on
% plot(x,Kr)
% legend('Data','ISRES+')
title('\it Kruppel')
set(gca,'FontSize',14)
hold off

nexttile
plot(x,Gt_E,'Linewidth',2)
hold on
plot(x,Gt_M,'Linewidth',2)
% hold on
% plot(x,Gt)
% legend('Data','ISRES+')
title('\it giant')
set(gca,'FontSize',14)
hold off

nexttile
plot(x,Kni_E,'Linewidth',2)
hold on
plot(x,Kni_M,'Linewidth',2)
% hold on
% plot(x,Kni)
% legend('Data','ISRES+')
title('\it knirps')
set(gca,'FontSize',14)
hold off

t.YLabel.String = 'Relative protein concentration';
t.XLabel.String = 'A-P position (% EL)';
t.Title.String = 'Gap gene model';
set(gca,'FontSize',14)