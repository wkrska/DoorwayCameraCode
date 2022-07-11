function imagesccyl(img, recon_grid, varargin)
    if (size(img,ndims(img)) == 1 || ndims(img) == 1 || size(img,1) == 1)
        img = cat(3,img, img, img);
    end
    img = mean(squeeze(reshape(img, [], 1,3)),2);
    imgnew = img-min(img(:));
    imgnew = imgnew/max(imgnew(:));
    colors = parula(1001);
    if (size(recon_grid,3) == 4)
        for i = 1:length(recon_grid)
            patch(squeeze(recon_grid(1,i,:)), squeeze(recon_grid(2,i,:)), colors(floor(max(imgnew(i),0)*1000)+1,:),'EdgeColor','none');
            hold on;
        end
    else
        for i = 1:length(recon_grid)
            if (img(i))
            line(squeeze(recon_grid(1,i,:)), squeeze(recon_grid(2,i,:)), 'Color', colors(floor(max(imgnew(i),0)*1000)+1,:), 'LineWidth', 5);
%             r = pdist(squeeze(recon_grid(:,i,:))');
%             xy = mean(recon_grid(:,i,:),3) - [r;r]/2;
%             rectangle('Position', [xy;r;r], 'Curvature', [1 1], 'EdgeColor', 'none', 'FaceColor', [min(1,max(0,img(i,1))),min(1,max(0,img(i,2))),min(1,max(0,img(i,3)))])
            hold on;
            end
        end
    end
    
    meas_size = .15;
    if (length(recon_grid) < length(img))
        rectangle('Position', [1.6-meas_size/2;-meas_size;meas_size;meas_size], ...
            'EdgeColor', 'none', 'FaceColor', colors(floor(max(imgnew(end),0)*1000)+1,:))
        rectangle('Position', [2.4-meas_size/2;-meas_size;meas_size;meas_size], ...
            'EdgeColor', 'none', 'FaceColor', colors(floor(max(imgnew(end),0)*1000)+1,:))
    end
    
    vals = linspace(min(img),max(img),11);
    colorbar('Ticks',[0:.1:1],'TickLabels',...
        {num2str(vals(1)), num2str(vals(2)), num2str(vals(3)), num2str(vals(4)), num2str(vals(5)), ...
        num2str(vals(6)), num2str(vals(7)), num2str(vals(8)), num2str(vals(9)), num2str(vals(10)), num2str(vals(11))});
end