function imagecyl(img, recon_grid)
if (size(img,ndims(img)) == 1 || ndims(img) == 1 || size(img,1) == 1)
    img = cat(3,img, img, img);
end

img_thresh = 0.05;

img = squeeze(reshape(img, [], 1,3));
if (size(recon_grid,3) == 4)
    for i = 1:length(recon_grid)
        if (max(img(i))>0)
            patch(squeeze(recon_grid(1,i,:)), squeeze(recon_grid(2,i,:)), [min(1,max(0,img(i,1))),min(1,max(0,img(i,2))),min(1,max(0,img(i,3)))],'EdgeColor','none');
            hold on;
        end
    end
else
    for i = 1:length(recon_grid)
        if (max(img(i,:,:))>img_thresh)
%             line(squeeze(recon_grid(1,i,:)), squeeze(recon_grid(2,i,:)), 'Color',...
%                 [1,1,1],...
%                 'LineWidth', 7);
            line(squeeze(recon_grid(1,i,:)), squeeze(recon_grid(2,i,:)), 'Color',...
                [min(1,max(0,img(i,1))),min(1,max(0,img(i,2))),min(1,max(0,img(i,3)))],...
                'LineWidth', 5);
            hold on;
        end
    end
end
meas_size = .15;
if (length(recon_grid) < length(img))
    rectangle('Position', [1.6-meas_size/2;-meas_size;meas_size;meas_size], ...
        'EdgeColor', 'none', 'FaceColor', [min(1,max(0,img(end,1))),min(1,max(0,img(end,2))),min(1,max(0,img(end,3)))])
    rectangle('Position', [2.4-meas_size/2;-meas_size;meas_size;meas_size], ...
        'EdgeColor', 'none', 'FaceColor', [min(1,max(0,img(end,1))),min(1,max(0,img(end,2))),min(1,max(0,img(end,3)))])
end
set(gca,'Color','k')
axis square
axis equal
end