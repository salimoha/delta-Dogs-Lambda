function cm = naturalskymap(n, withbackground)
% naturalskymap is a colormap for cloudy skies
%	naturalskymap(N) is an N-by-3 matrix containing a colormap appropriate
%	for drawing a sky.  By default, the first item in the colormap is all
%	white in order to give a white background to non-square cloudmaps.
%
%	naturalskymap(N, withbackground) allows the disabling of this white
%	background.
%
%	If no N is specified, N will be inferred from the current figure's colormap.
%
%	See Also: colormap, prettyplot.

% Input handling
if(nargin < 1 || isempty(n))
	cf = get(0,'CurrentFigure');
	if(~isempty(cf))
		n = size(get(cf,'Colormap'),1);
	else
		n = 64;
	end
end
if(nargin < 2)
	withbackground = true;
end

% base colors for background, clear sky, bright clouds, and dark clouds
backgroundcolor = [1  1  1];
basecolors = [ 0 .4 .8; .9 .9 .9; .7 .7 .7];
% take care of the background
if(withbackground)
	cm = backgroundcolor;
	n = n-1;
else
	cm = zeros(0,3);
end
% add other colors
if(n<=3)
	cm = [cm; basecolors(1:n,:)];
else
	% number of colors to insert in each range
	nfirst = floor((n-3)/2);
	nsecond = ceil((n-3)/2);
	% interpolate
	nfirst = (0:nfirst)'/(nfirst+1);
	cleartothin = (1-nfirst)*basecolors(1,:) + nfirst*basecolors(2,:);
	nsecond = (0:(nsecond+1))'/(nsecond+1); % go all the way to 1 this time
	thintothick = (1-nsecond)*basecolors(2,:) + nsecond*basecolors(3,:);
	cm = [cm; cleartothin; thintothick];
end

end
