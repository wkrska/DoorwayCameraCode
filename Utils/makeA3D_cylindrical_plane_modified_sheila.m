function [A, recon_grid, S] = makeA3D_cylindrical_plane_modified_sheila(pos,oversample,discr,p,plotstuff)
% flip = 0, wall is from left side, flip = 1 wall is from right
mp = p.meas_pix * oversample;

% Define plane representing each pix
plane_height = 2;

% Pixels of hidden space: Closest vertex, left vertex, furthest vertex, right vertex
a_left = linspace(0, pi, p.rec_size(1)+1);
remove = find(a_left > (pi-p.min_angle));
a_left(remove) = [];
a_right = linspace(0, pi, p.rec_size(2)+1);
remove = find(a_right < p.min_angle);
a_right(remove) = [];

disp("Constructing grid");
rec_pix=0; % Number of reconstruction pixels total
dead_radius = .1; % Radius from each corner to exclude
rx = zeros(1,2);
ry = zeros(1,2);

% Save pixels as groups, to be used to form S
indicies = zeros(p.rec_size(1),p.rec_size(2));

% % % % % % % Build recon grid % % % % % % %
for i = 1:length(a_left)-1
    for j = 1:length(a_right)-1
        if (a_right(j) > a_left(i))
            verts=zeros(2,2);
            verts(1,:)=[pos(2,1)+cot(a_right(j+1))*(pos(1,1)-pos(2,1))/(cot(a_right(j+1))-cot(a_left(i+1))),(pos(1,1)-pos(2,1))/(cot(a_right(j+1))-cot(a_left(i+1)))];
            verts(2,:)=[pos(2,1)+cot(a_right(j  ))*(pos(1,1)-pos(2,1))/(cot(a_right(j  ))-cot(a_left(i  ))),(pos(1,1)-pos(2,1))/(cot(a_right(j  ))-cot(a_left(i  )))];
            
            % Condition 1: Within the box
            midpnt = mean(verts);
            condition_1 = (midpnt(1)>=0 && midpnt(1)<=p.room_dim(1) && midpnt(2)>=0 && midpnt(2)<=p.room_dim(2));
            
            %                 % Condition 2: Within the size constraint (m^2)
            %                 threshold = 2.8; % Length in cm
            %                 condition_2 = pdist(verts, 'euclidean') >= (threshold / 100);
            
            % Condition 3: Minimum distance away
            condition_3 = (min(table(pdist([verts;pos(1,:)])).Var1(2:end)) >= dead_radius) && ...
                (min(table(pdist([verts;pos(2,:)])).Var1(2:end)) >= dead_radius);
            
            %                 % Condition 4: Cannot be between corners
            %                 condition_4 = ((sum(verts(:,1)>=pos(2,1))==2) || (sum(verts(:,1)<=pos(1,1))==2)|| (sum(verts(:,2)>=(pos(2,2)+dead_radius))==2));
            
            % Condition 5: Not within a certain ellipse
            threshold = 1.3;
            condition_5 = (pdist([midpnt; pos(1,:)])+pdist([midpnt; pos(2,:)])) > threshold;
            
            if (condition_1 && condition_3 && condition_5)
                rec_pix = rec_pix+1;
                rx(rec_pix,:) = verts(:,1)';
                ry(rec_pix,:) = verts(:,2)';
                
                % Record both angle groups in groups matrix
                indicies(i,j)=rec_pix;
                
                if (j==1)
                    rec_pix = rec_pix+1;
                    rx(rec_pix,:) = [rx(rec_pix-1,2),rx(rec_pix-1,2)];
                    ry(rec_pix,:) = [ry(rec_pix-1,2),0];
                    indicies(i,length(a_right))=rec_pix;
                end
                if (i==length(a_left)-1)
                    rec_pix = rec_pix+1;
                    rx(rec_pix,:) = [rx(rec_pix-1,1),rx(rec_pix-1,1)];
                    ry(rec_pix,:) = [0,ry(rec_pix-1,1)];
                    indicies(length(a_left),j)=rec_pix;
                end
            end
        end
    end
end

% Reconstruction grid to be used, returnable value
recon_grid(1,:,:) = rx(:,:);
recon_grid(2,:,:) = ry(:,:);
fprintf("%d reconstruction pixels\n", rec_pix);


% % % % % % % Construct sparsity matricies % % % % % % %
S = zeros(rec_pix*2);
    
style = "binary";

% Binary on/off
if (strcmp("binary",style))
    disp("Constructing sparsity matrix S (binary)");
    for i = 1:max(p.rec_size) % for each angle 
        elements_l = nonzeros(indicies(i,:)); % Get list of non-zero elements in each row i.e every pixel in that radial column
        elements_r = nonzeros(indicies(:,i)); % Get list of non-zero elements in each column i.e every pixel in that radial column
        for j = elements_l % For each element
           S(j,elements_l) = 1; % Set intersection of given element and every other element to 1 (in row)
           S(elements_l,j) = 1; % Set intersection of given element and every other element to 1 (in column)
        end
        for j = elements_r
           S(j,elements_r) = 1;
           S(elements_r,j) = 1;
        end
    end
    for i = 1:rec_pix % Make diagonal 0's
        S(i,i) = 0;
    end
end

% Tapered
if (strcmp("tapered",style))
    disp("Constructing sparsity matrix S (tapered)");
    % Weight by distance [1 away, 2 away, ... n+ away]
    dist_weight = [0,.25,.75, 1];
    for i = 1:rec_pix % for each angle
        [l,r] = find((indicies==i)); % current angular coords
        elements_l = nonzeros(indicies(l,:)); % Get list of non-zero elements in each row i.e every pixel in that radial column
        elements_r = nonzeros(indicies(:,r)); % Get list of non-zero elements in each column i.e every pixel in that radial column
        for j = 1:length(elements_l) % For each element
            [cl,cr]=find((indicies==elements_l(j)));
            S(i,elements_l(j)) = dist_weight(max(min(length(dist_weight),abs(cr-r)),1)); 
            S(elements_l(j),i) = dist_weight(max(min(length(dist_weight),abs(cr-r)),1));
        end
        for j = 1:length(elements_r) % For each element
            [cl,cr]=find((indicies==elements_r(j)));
            S(i,elements_r(j)) = dist_weight(max(min(length(dist_weight),abs(cl-l)),1)); 
            S(elements_r(j),i) = dist_weight(max(min(length(dist_weight),abs(cl-l)),1));
        end
    end
    for i = 1:rec_pix % Make diagonal 0's 
        S(i,i) = 0;
    end
end

S(end/2+1:end,end/2+1:end) = S(1:end/2,1:end/2);

% S matrix shared across corners
for i = 1:rec_pix
    [l1,r1]=find((indicies==i)); % Left angle, right angle of i
    for j = 1:rec_pix
        [l2,r2]=find((indicies==j)); % Left angle, right angle of j
        % Left side perspective
        if (l1==l2 && r1<r2)
            S(i,end/2+j) = 1;
            S(end/2+j,i) = 1;
        end
        % Right side perspective
        if (r1==r2 && l1<l2)
            S(i,end/2+j) = 1;
            S(end/2+j,i) = 1;
        end
    end
end

% S matrix shared across corners
for i = 1:rec_pix
    [l1,r1]=find((indicies==i));
    for j = 1:rec_pix
        [l2,r2]=find((indicies==j));
        % Left side perspective
        if (l1==l2 && r1>r2)
            S(i,rec_pix+j) = 1;
        end
        % Right side perspective
        if (r1==r2 && l1<l2)
            S(rec_pix+j,i) = 1;
        end
    end
end

% % % % % % % Calculate A % % % % % % %
% Pixels of left measurement * oversample
[cx_left,cy_left] = meshgrid(linspace(pos(1,1)-p.meas_size(1)/2, pos(1,1)+p.meas_size(1)/2, mp),linspace(pos(1,2), pos(1,2)-p.meas_size(1), mp));

% Makes single vector
cx_left = cx_left(:);
cy_left = cy_left(:);
cz = zeros(mp^2,1);
% cpos_left = gpuArray([cx_left, cy_left, cz]);
cpos_left = [cx_left, cy_left, cz];

% Pixels of right measurement * oversample
[cx_right,cy_right] = meshgrid(linspace(pos(2,1)-p.meas_size(1)/2, pos(2,1)+p.meas_size(1)/2, mp),linspace(pos(2,2), pos(2,2)-p.meas_size(1), mp));

% Makes single vector
cx_right = cx_right(:);
cy_right = cy_right(:);
% cpos_right = gpuArray([cx_right, cy_right, cz]);
cpos_right = [cx_right, cy_right, cz];

% Angle from each meas pix to both corners
cangle_left = atan((cx_left-pos(1,1))./(-eps+cy_left-pos(1,2)));
cangle_left2 = atan((cx_left-pos(2,1))./(-eps+cy_left-pos(2,2)));
cangle_right = atan((cx_right-pos(2,1))./(-eps+cy_right-pos(2,2)));
cangle_right2 = atan((cx_right-pos(1,1))./(-eps+cy_right-pos(1,2)));

A_left = zeros(p.meas_pix*p.meas_pix,length(rx()));
A_right = zeros(p.meas_pix*p.meas_pix,length(rx()));

% % % % % % % For each hidden pixel... % % % % % % %
parfor po = 1:length(rx)
%     meas_left = zeros(p.meas_pix * oversample * p.meas_pix * oversample,1,'gpuArray');
%     meas_right = zeros(p.meas_pix * oversample * p.meas_pix * oversample,1,'gpuArray');
    meas_left = zeros(p.meas_pix * oversample * p.meas_pix * oversample,1);
    meas_right = zeros(p.meas_pix * oversample * p.meas_pix * oversample,1);  
    % Calculate normal
    normal_surface = [-(ry(po,1)-ry(po,2)), rx(po,1)-rx(po,2),0];
    normal_surface = normal_surface/sqrt(sum(normal_surface.^2));
    
    % Calculate plane area of each patch
    plane_area_eff = (1/discr^2)*plane_height*pdist([rx(po,:)',ry(po,:)'], 'euclidean');
    
    % Organized Discretization: discr*discr grid
    discr_x = linspace(rx(po,1),rx(po,2),discr+1);
    discr_x = discr_x(1:end-1)+(discr_x(2)-discr_x(1))/2;
    discr_y = linspace(ry(po,1),ry(po,2),discr+1);
    discr_y = discr_y(1:end-1)+(discr_y(2)-discr_y(1))/2;
    discr_z = linspace(0,plane_height,discr+1);
    discr_z = discr_z(1:end-1)+(discr_z(2)-discr_z(1))/2;
    [px, pz] = meshgrid(discr_x,discr_z);
    [py, pz] = meshgrid(discr_y,discr_z);
    
    %         % Random Discretization: random distribution of discr^2 points
    %         rand_nums = rand(discr^2,1)
    %         px = rand_nums*(rx(po,4)-rx(po,2))+rx(po,2);
    %         py = rand_nums*(ry(po,4)-ry(po,2))+ry(po,2);
    %         pz = rand(discr^2,1)*plane_height;
    
    % render each patch
    for pt = 1:length(px(:))
        % Determine if it is visible
        rangle = atan((px(pt)-pos(1,1))./(py(pt)-pos(1,2)));
        rangle2 = atan((px(pt)-pos(2,1))./(py(pt)-pos(2,2)));
        
        vis_left = (cangle_left<rangle);
        vis_left2 = (cangle_left2>rangle2);
%         vis_left = gpuArray(vis_left(:).*vis_left2(:));
        vis_left = vis_left(:).*vis_left2(:);
        
        vis_right = (cangle_right>rangle2);
        vis_right2 = (cangle_right2<rangle);
%         vis_right = gpuArray(vis_right(:).*vis_right2(:));
        vis_right = vis_right(:).*vis_right2(:);
        
        
        % Intensity
        vec_left = cpos_left - [px(pt), py(pt), pz(pt)];
        m_left = plane_area_eff*max(0,sum(vec_left.*normal_surface,2)).*max(0,(-vec_left(:,3)))./sum(vec_left.^2,2).^(3/2);
        
        vec_right = cpos_right - [px(pt), py(pt), pz(pt)];
        m_right = plane_area_eff*max(0,sum(vec_right.*normal_surface,2)).*max(0,(-vec_right(:,3)))./sum(vec_right.^2,2).^(3/2);
        
        % Add contribution to meas from pixel patch
        meas_left = meas_left + m_left.*vis_left;
        meas_right = meas_right + m_right.*vis_right ;
    end
    
    if (sum(meas_left<0)>0 || sum(meas_left<0)>0)
        disp("Negative!")
    end
    
    
    %         A_left(:,po)  = gather(reshape(imresize(reshape( meas_left,mp,mp),1/oversample),[],1));
    %         A_right(:,po) = gather(reshape(imresize(reshape(meas_right,mp,mp),1/oversample),[],1));
%     A_left(:,po)  = reshape(imresize(gather(reshape( meas_left,mp,mp)),1/oversample,'bilinear','Antialiasing',false),[],1);
%     A_right(:,po) = reshape(imresize(gather(reshape(meas_right,mp,mp)),1/oversample,'bilinear','Antialiasing',false),[],1);

    A_left(:,po)  = reshape(imresize(reshape( meas_left,mp,mp),1/oversample,'bilinear','Antialiasing',false),[],1);
    A_right(:,po) = reshape(imresize(reshape(meas_right,mp,mp),1/oversample,'bilinear','Antialiasing',false),[],1);
    
    if (plotstuff)
        % Plot progress
        subplot(1,3,1)
        hold on
        line(squeeze(rx(po,:)),squeeze(ry(po,:)), 'Color', [.3,1,.3], 'LineWidth', 5);
        %             rectangle('Position', [[mean(rx(po,:)) mean(ry(po,:))],pdist([rx(po,:)',ry(po,:)']),pdist([rx(po,:)',ry(po,:)'])], 'Curvature', [1 1], 'EdgeColor', 'none', 'FaceColor', [.3,1,.3])
        subplot(1,3,2)
        imagesc(reshape(A_left(:,po),p.meas_pix,p.meas_pix))
        subplot(1,3,3)
        imagesc(reshape(A_right(:,po),p.meas_pix,p.meas_pix))
        drawnow
    end
    
    % Display progress
    if (mod(po,50)==0)
        fprintf("po=%d (%.2f%%)\n",po,po/rec_pix*100)
    end
end
%
%     % Discard nearest pixels within N pixels
%     dead_radius = 10;
%     mask = zeros(p.meas_pix);
%     for i = 1:dead_radius
%         mask(i,:) = sqrt(([1:p.meas_pix]-p.meas_pix/2).^2+i^2) < dead_radius;
%     end
%     mask = ~mask(:);
%     A_left = A_left.*mask;
%     A_right = A_right.*mask;
A = [A_left;A_right];
delete(gcp('nocreate'));
end