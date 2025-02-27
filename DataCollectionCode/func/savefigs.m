function savefig(h, filename, outdir, fontsize, aspect)
set(gca,'FontSize',fontsize);
set(gca,'TickDir','out');
set(gca,'LineWidth',1.4); 
xl = get(gca,'XLabel');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', fontsize-4)
set(xl, 'FontSize', fontsize);
yl = get(gca,'YLabel');
yAX = get(gca,'YAxis');
set(yAX,'FontSize', fontsize-4)
set(yl, 'FontSize', fontsize);
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 aspect];
saveas(h,fullfile(outdir,sprintf('%s.eps',filename)),'epsc2');
saveas(h,fullfile(outdir,sprintf('%s.png',filename)),'png');