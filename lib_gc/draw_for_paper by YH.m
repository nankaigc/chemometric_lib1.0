%fig_Y_pred(cal_y,cvpredpls(:,:,LVs),'C')
set(gcf,'Position',[100 100 260 220]);

figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',10),'FontSize',figure_FontSize);

clear figure_FontSize