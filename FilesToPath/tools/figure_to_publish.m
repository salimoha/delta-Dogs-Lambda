function figure_to_publish(cmfile,FigTitle)
%make the figures ready to publ;ish
% try to save the figure in eps form
 set(gca, 'XGrid', 'off', 'YGrid', 'off',  'FontSize', 14);
 
% set(gca,   'FontSize', 18);
% set(gcf, 'Position' , [60 60 700 500], 'color' , [1 1 1], 'PaperPositionMode' , 'auto', 'InverthardCopy' , 'off' );
% set(gcf, 'color' , [1 1 1], 'PaperPositionMode' , 'auto', 'InverthardCopy' , 'on' );
set(gcf, 'color' , [1 1 1], 'PaperPositionMode' , 'auto', 'InverthardCopy' , 'on' );
% set(gcf,'units','normalized', 'outerposition', [0,0,0.5,0.5])

% if nargin <= 3
%     saveOpt=1;
% end
% if nargin <= 1
%     cmfile='./'
% end

if nargin >=2
%      FigTitle = char(['Time: ' dayname ' with cloud fraction limit of '  fixVar(2)  num2str(int32(cfHigh),'%02d') ]);
 title(FigTitle)
end  


if nargin >=1
% cmfile = [startTime '_Cf_limits_upper']
             saveas(gcf,cmfile, 'fig');       
             saveas(gcf,cmfile, 'png');
             saveas(gcf,cmfile, 'epsc2'); 
%             close(gcf);
end

% close 
