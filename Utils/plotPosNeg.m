function [disp_pos, disp_neg] = plotPosNeg(disp)

disp_pos = disp;
disp_neg = disp;
disp_pos(disp_pos<0)=0;
disp_neg(disp_neg>0)=0;
disp_neg = abs(disp_neg);

end

