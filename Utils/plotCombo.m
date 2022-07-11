function [comboView] = plotCombo(leftView,rightView,c1, c2,t)
mask1 = zeros(size(leftView));
mask2 = zeros(size(rightView));



mv = t*max(abs([c1(:); c2(:)]));
ifilt1 = find(c1>mv);
ifilt2 = find(c2>mv);

% 
% mv = t*max([leftView(:); rightView(:)]);
% ifilt1 = find(leftView>mv);
% ifilt2 = find(rightView>mv);

% mask1(ifilt1) = 1;
% mask2(ifilt2) = 1;

mask1(ifilt1,:) = 1;
mask2(ifilt2,:) = 1;

mask = mask1.*mask2;

maxView = max(leftView, rightView);

comboView = mask.*maxView;


end

