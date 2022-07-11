function im_out = myimresize(im_in, res_out)
    out_y = res_out(1);
    out_x = res_out(2);
    [in_x, in_y, in_c] = size(im_in);
    samp_x = round(linspace(1,in_x,out_x+1)); samp_x = samp_x(1:end-1)+diff(samp_x(1:2))/2;
    samp_y = round(linspace(1,in_y,out_y+1)); samp_y = samp_y(1:end-1)+diff(samp_y(1:2))/2;
    im_in_LPF = cat(3,...
        imgaussfilt(im_in(:,:,1),'FilterSize',2*floor(min(in_y/out_y,in_x/out_x)/2)+1),...
        imgaussfilt(im_in(:,:,2),'FilterSize',2*floor(min(in_y/out_y,in_x/out_x)/2)+1),...
        imgaussfilt(im_in(:,:,3),'FilterSize',2*floor(min(in_y/out_y,in_x/out_x)/2)+1));
    im_out = im_in_LPF(samp_y,samp_x,:);
end