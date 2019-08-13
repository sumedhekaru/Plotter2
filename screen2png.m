function screen2png(fig,filename)
%SCREEN2JPEG Generate a JPEG file of 'fig' figure with
% dimensions consistent with the figure's screen dimensions.
%
% SCREEN2PNG('filename') saves the current figure to the
% JPEG file "filename".
%
% Sean P. McCarthy
% Copyright (c) 1984-98 by MathWorks, Inc. All Rights Reserved


oldscreenunits = get(fig,'Units');
oldpaperunits = get(fig,'PaperUnits');
oldpaperpos = get(fig,'PaperPosition');
set(fig,'Units','pixels');
scrpos = get(fig,'Position');
newpos = scrpos/100;
set(fig,'PaperUnits','inches',...
'PaperPosition',newpos)
print('-dpng', filename, '-r100');
drawnow
set(fig,'Units',oldscreenunits,...
'PaperUnits',oldpaperunits,...
'PaperPosition',oldpaperpos)