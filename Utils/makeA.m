%% Function Details:
% 2.5 Dimensional (fixed height objects)
% Uses custom "cylindrical" scene grid
% Planar scene facets lying on intersections of grid
% Discards facets too close to corner, and at too shallow an angle
% Accounts for wall thickness
% Calculates effects of doorframe reflections
% DOES NOT calculate MEGS
% DOES NOT account for wall reflections


function [A, recon_grid, nAmb] = makeA(pos,oversample,discr,p,plotstuff)

%%%% Definition of parameters
% pos: [x_1 y_1; x_2 y_2] defines left far corner, and right near corner of doorframe, under which the measurements are centered
% oversample: distretization of pixles
% discr: distretization of scene facets (in both x and z)
% plotstuff: currently useless bcs parallellized

% Define corners as c l/r n/f(x,y)
pos_lf = pos(1,:);
pos_ln = [pos(1,1),pos(2,2)];
pos_rn = pos(2,:);
pos_rf = [pos(2,1),pos(1,2)];

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
nAmb = 0; % Number of "ambient" or non-standard fascets in scene
dead_radius = .1; % Radius from each corner to exclude
rx = zeros(1,2);
ry = zeros(1,2);

% % % % % % % Build recon grid % % % % % % %
for i = 1:length(a_left)-1
    for j = 1:length(a_right)-1
        if (a_right(j) > a_left(i))
            
            % Build grid starting at far corners
            verts=zeros(2,2);
            verts(1,:)=[pos_rf(1)+cot(a_right(j+1))*(pos_lf(1)-pos_rf(1))/(cot(a_right(j+1))-cot(a_left(i+1))),(pos_lf(1)-pos_rf(1))/(cot(a_right(j+1))-cot(a_left(i+1)))];
            verts(2,:)=[pos_rf(1)+cot(a_right(j  ))*(pos_lf(1)-pos_rf(1))/(cot(a_right(j  ))-cot(a_left(i  ))),(pos_lf(1)-pos_rf(1))/(cot(a_right(j  ))-cot(a_left(i  )))];
            
            % Condition 1: Within the box
            midpnt = mean(verts);
            condition_1 = (midpnt(1)>=0 && midpnt(1)<=p.room_dim(1) && midpnt(2)>=0 && midpnt(2)<=p.room_dim(2));
            
            %                 % Condition 2: Within the size constraint (m^2)
            %                 threshold = 2.8; % Length in cm
            %                 condition_2 = pdist(verts, 'euclidean') >= (threshold / 100);
            
            % Condition 3: Minimum distance away
            condition_3 = (min(table(pdist([verts;pos_lf])).Var1(2:end)) >= dead_radius) && ...
                (min(table(pdist([verts;pos_rf])).Var1(2:end)) >= dead_radius);
            
            %                 % Condition 4: Cannot be between corners
            %                 condition_4 = ((sum(verts(:,1)>=pos_rf(1))==2) || (sum(verts(:,1)<=pos_lf(1))==2)|| (sum(verts(:,2)>=(pos_rf(2)+dead_radius))==2));
            
            % Condition 5: Not within a certain ellipse
            threshold = 1.3;
            condition_5 = (pdist([midpnt; pos_lf])+pdist([midpnt; pos_rf])) > threshold;
            
            if (condition_1 && condition_3 && condition_5)
                rec_pix = rec_pix+1;
                rx(rec_pix,:) = verts(:,1)';
                ry(rec_pix,:) = verts(:,2)';
                
                if (j==1)
                    rec_pix = rec_pix+1;
                    rx(rec_pix,:) = [rx(rec_pix-1,2),rx(rec_pix-1,2)];
                    ry(rec_pix,:) = [ry(rec_pix-1,2),pos_lf(2)];
                end
                if (i==length(a_left)-1)
                    rec_pix = rec_pix+1;
                    rx(rec_pix,:) = [rx(rec_pix-1,1),rx(rec_pix-1,1)];
                    ry(rec_pix,:) = [pos_lf(2),ry(rec_pix-1,1)];
                end
            end
        end
    end
end

% % Calculate Doorframe
% rec_pix = rec_pix+1;    % Left
% nAmb = nAmb+1;
% rx(rec_pix,:) = [pos_ln(1),pos_ln(1)];
% ry(rec_pix,:) = [pos_ln(2),pos_lf(2)];
% rec_pix = rec_pix+1;    % Right
% nAmb = nAmb+1;
% rx(rec_pix,:) = [pos_rn(1),pos_rn(1)];
% ry(rec_pix,:) = [pos_rn(2),pos_rf(2)];


% Reconstruction grid to be used, returnable value
recon_grid(1,:,:) = rx(:,:);
recon_grid(2,:,:) = ry(:,:);
fprintf("%d reconstruction pixels\n", rec_pix);


% % % % % % % Construct sparsity matricies % % % % % % %

% % % % % % % Calculate A % % % % % % %
% Pixels of measurement * oversample
[cx_left,cy_left] = meshgrid(linspace(pos_ln(1)-p.meas_size(1)/2, pos_ln(1)+p.meas_size(1)/2, mp),linspace(pos_ln(2), pos_ln(2)-p.meas_size(1), mp));
[cx_right,cy_right] = meshgrid(linspace(pos_rn(1)-p.meas_size(1)/2, pos_rn(1)+p.meas_size(1)/2, mp),linspace(pos_rn(2), pos_rn(2)-p.meas_size(1), mp));

% Makes single vector
cx_left = cx_left(:);
cy_left = cy_left(:);

cx_right = cx_right(:);
cy_right = cy_right(:);

cz = zeros(mp^2,1);

cpos_left = gpuArray([cx_left, cy_left, cz]);
cpos_right = gpuArray([cx_right, cy_right, cz]);

% Angle from each meas pix to both corners, both front and back
ca_ln = atan((cx_left-pos_ln(1))./(-eps+cy_left-pos_ln(2)));  % e.g. angle from left meas to left near corner
ca_lf = atan((cx_left-pos_lf(1))./(-eps+cy_left-pos_lf(2)));
ca_ln_2 = atan((cx_left-pos_rn(1))./(-eps+cy_left-pos_rn(2)));  % e.g. angle from left meas to right near corner
ca_lf_2 = atan((cx_left-pos_rf(1))./(-eps+cy_left-pos_rf(2)));
ca_rn = atan((cx_right-pos_rn(1))./(-eps+cy_right-pos_rn(2)));
ca_rf = atan((cx_right-pos_rf(1))./(-eps+cy_right-pos_rf(2)));
ca_rn_2 = atan((cx_right-pos_ln(1))./(-eps+cy_right-pos_ln(2)));
ca_rf_2 = atan((cx_right-pos_lf(1))./(-eps+cy_right-pos_lf(2)));

A_left = zeros(p.meas_pix*p.meas_pix,length(rx()));
A_right = zeros(p.meas_pix*p.meas_pix,length(rx()));

% % % % % % % For each hidden pixel... % % % % % % %
parfor po = 1:length(rx)
    meas_left = zeros(p.meas_pix * oversample * p.meas_pix * oversample,1,'gpuArray');
    meas_right = zeros(p.meas_pix * oversample * p.meas_pix * oversample,1,'gpuArray');

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
    if (length(rx)-po<nAmb)
        discr_x = linspace(rx(po,1),rx(po,2),discr*2+1);
        discr_x = discr_x(1:end-1)+(discr_x(2)-discr_x(1))/2;
        discr_y = linspace(ry(po,1),ry(po,2),discr*2+1);
        discr_y = discr_y(1:end-1)+(discr_y(2)-discr_y(1))/2;
        discr_z = linspace(0,plane_height,discr/2+1);
        discr_z = discr_z(1:end-1)+(discr_z(2)-discr_z(1))/2;
    end
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
        
        % angle between corner and fascet discr
        ra_ln = atan((px(pt)-pos_ln(1))./(py(pt)-pos_ln(2)));
        ra_lf = atan((px(pt)-pos_lf(1))./(py(pt)-pos_lf(2)));
        ra_rn = atan((px(pt)-pos_rn(1))./(py(pt)-pos_rn(2)));
        ra_rf = atan((px(pt)-pos_rf(1))./(py(pt)-pos_rf(2)));
        
        vis_left = (max(ca_ln,ca_lf)<min(ra_ln,ra_lf));
        vis_left2 = (min(ca_ln_2,ca_lf_2)>max(ra_rn,ca_rf));
        vis_left = gpuArray(vis_left(:).*vis_left2(:));
        
        vis_right = (min(ca_rn,ca_rf)>max(ra_rn,ra_rf));
        vis_right2 = (max(ca_rn_2,ca_rf_2)<min(ra_ln,ra_lf));
        vis_right = gpuArray(vis_right(:).*vis_right2(:));
        
        
        % Intensity
        vec_left = cpos_left - [px(pt), py(pt), pz(pt)];
        m_left = plane_area_eff*max(0,sum(vec_left.*normal_surface,2)).*max(0,(-vec_left(:,3)))./sum(vec_left.^2,2).^(3/2);
        
        vec_right = cpos_right - [px(pt), py(pt), pz(pt)];
        m_right = plane_area_eff*max(0,sum(vec_right.*normal_surface,2)).*max(0,(-vec_right(:,3)))./sum(vec_right.^2,2).^(3/2);
        
        % Add contribution to meas from pixel patch
        meas_left = meas_left + m_left.*vis_left;
        meas_right = meas_right + m_right.*vis_right ;
    end
    
    if (sum(meas_left<0)>0 || sum(meas_right<0)>0)
        disp("Negative!")
    end
    
    A_left(:,po)  = reshape(imresize(gather(reshape( meas_left,mp,mp)),1/oversample,'bilinear','Antialiasing',false),[],1);
    A_right(:,po) = reshape(imresize(gather(reshape(meas_right,mp,mp)),1/oversample,'bilinear','Antialiasing',false),[],1);
    
    % Display progress
    if (mod(po,50)==0)
        fprintf("po=%d (%.2f%%)\n",po,po/rec_pix*100)
    end
end

% % % % % % % Calculate Doorframe % % % % % % %
% rec_pix = rec_pix+1;    % Left
% rx(rec_pix,:) = [pos_ln(1),pos_ln(1)];
% ry(rec_pix,:) = [pos_ln(2),pos_lf(2)];
% rec_pix = rec_pix+1;    % Right
% rx(rec_pix,:) = [pos_rn(1),pos_rn(1)];
% ry(rec_pix,:) = [pos_rf(2),pos_rn(2)];
% % rec_pix = rec_pix+1;    % Top
% % rx(rec_pix,:) = [rx(rec_pix-1,2),rx(rec_pix-1,2)];
% % ry(rec_pix,:) = [ry(rec_pix-1,2),pos_lf(2)];


for po = length(rx)-2:length(rx)
    meas_left = zeros(p.meas_pix * oversample * p.meas_pix * oversample,1,'gpuArray');
    meas_right = zeros(p.meas_pix * oversample * p.meas_pix * oversample,1,'gpuArray');
    
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
        
        % angle between corner and fascet discr
        ra_ln = atan((px(pt)-pos_ln(1))./(py(pt)-pos_ln(2)));
        ra_lf = atan((px(pt)-pos_lf(1))./(py(pt)-pos_lf(2)));
        ra_rn = atan((px(pt)-pos_rn(1))./(py(pt)-pos_rn(2)));
        ra_rf = atan((px(pt)-pos_rf(1))./(py(pt)-pos_rf(2)));
        
        vis_left = (max(ca_ln,ca_lf)<min(ra_ln,ra_lf));
        vis_left2 = (min(ca_ln_2,ca_lf_2)>max(ra_rn,ca_rf));
        vis_left = gpuArray(vis_left(:).*vis_left2(:));
        
        vis_right = (min(ca_rn,ca_rf)>max(ra_rn,ra_rf));
        vis_right2 = (max(ca_rn_2,ca_rf_2)<min(ra_ln,ra_lf));
        vis_right = gpuArray(vis_right(:).*vis_right2(:));
        
        
        % Intensity
        vec_left = cpos_left - [px(pt), py(pt), pz(pt)];
        m_left = plane_area_eff*max(0,sum(vec_left.*normal_surface,2)).*max(0,(-vec_left(:,3)))./sum(vec_left.^2,2).^(3/2);
        
        vec_right = cpos_right - [px(pt), py(pt), pz(pt)];
        m_right = plane_area_eff*max(0,sum(vec_right.*normal_surface,2)).*max(0,(-vec_right(:,3)))./sum(vec_right.^2,2).^(3/2);
        
        % Add contribution to meas from pixel patch
        meas_left = meas_left + m_left.*vis_left;
        meas_right = meas_right + m_right.*vis_right ;
    end
    
    if (sum(meas_left<0)>0 || sum(meas_right<0)>0)
        disp("Negative!")
    end
    
    A_left(:,po)  = reshape(imresize(gather(reshape( meas_left,mp,mp)),1/oversample,'bilinear','Antialiasing',false),[],1);
    A_right(:,po) = reshape(imresize(gather(reshape(meas_right,mp,mp)),1/oversample,'bilinear','Antialiasing',false),[],1);
    
    % Display progress
    if (mod(po,50)==0)
        fprintf("po=%d (%.2f%%)\n",po,po/rec_pix*100)
    end
end
A = [A_left;A_right];
delete(gcp('nocreate'));
end