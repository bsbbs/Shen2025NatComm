function mysavefig(h, filename, outdir, fontsize, aspect)
set(gca,'FontSize',fontsize);
set(gca,'TickDir','out');
ax = gca;
ticklength = 1.6;
ax.TickLength = ax.TickLength*ticklength;
% Set the axes line color to black
ax.XColor = 'black'; % For x-axis
ax.YColor = 'black'; % For y-axis

% Set the tick color to black (this will change the color of the tick labels)
ax.TickLabelInterpreter = 'tex';
ax.XAxis.TickLabelFormat = '\\color{black} %g';
ax.YAxis.TickLabelFormat = '\\color{black} %g';

%set(gca,'LineWidth',1); 
% xl = get(gca,'XLabel');
% xAX = get(gca,'XAxis');
% set(xAX,'FontSize', fontsize-2)
% set(xl, 'FontSize', fontsize);
% yl = get(gca,'YLabel');
% yAX = get(gca,'YAxis');
% set(yAX,'FontSize', fontsize-2)
% set(yl, 'FontSize', fontsize);
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 aspect];
saveas(h,fullfile(outdir,sprintf('%s.pdf',filename)),'pdf');