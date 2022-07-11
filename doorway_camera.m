function doorway_camera(labData,scnNum,samps,facetHeight,fov,res,lambdaGroup,showPlots, bgRemoval)
addpath('Utils','Laboratory Data','Synthetic Data');
close all
%% Initialize test parameters

% Struct p room and measurement 
p.plane_height = facetHeight; % height of facets in meters
p.meas_pix = res; % meas_pix x meas_pix measurement at FOV
p.meas_size = fov; % FOV size in meters (must be square)
p.rec_size = [45,45]; % Number of angles per corner, yields up to N*(N-1)/2 pixels
p.room_dim = [3,3]; % Size of the hidden scene in meters
p.min_angle = pi/8; % 

door_width = 0.508; % in meters
door_center = p.room_dim(2)/2;
pos = [door_center-door_width/2, 0; door_center+door_width/2, 0];

%% Load forward model for each corner

fwd_model_file = strcat('Utils\',num2str(p.meas_pix),'p-',num2str(p.room_dim(1)),'m-',num2str(p.rec_size(1)),'r-',num2str(round(p.meas_size*100)),'cm-min-pi_',num2str(pi/p.min_angle),'.mat');

disp(strcat("Loading Forward Model: ",fwd_model_file));
try
    load(fwd_model_file);
catch
    disp("Forward model not found, generating A matrix");
    tic
    [A, recon_grid] = makeA(pos,1,16,p,1);
    toc
    save(fwd_model_file, 'p','A','recon_grid','pos');
end

A_ext = [A, ones(size(A,1),1)];

%% figure out which pixels can only be positive

fg_range = p.room_dim(1)*.2; % everything inside of this range must be positive
center = [door_center,0];
rangeEdge1 = sqrt((squeeze(recon_grid(1,:,1))-center(1)).^2 + (squeeze(recon_grid(2,:,1))-center(2)).^2);
rangeEdge2 = sqrt((squeeze(recon_grid(1,:,2))-center(1)).^2 + (squeeze(recon_grid(2,:,2))-center(2)).^2);
y_coord = min([squeeze(recon_grid(2,:,1))',squeeze(recon_grid(2,:,2))'],[],2);
x_coord = min([squeeze(abs(recon_grid(1,:,1)-center(1)))',squeeze(abs(recon_grid(1,:,2)-center(1)))'],[],2);
ifiltx = find(x_coord<p.room_dim(1)*1/3);
ifilty = find(y_coord<p.room_dim(1)*15/24);
ifilt = intersect(ifiltx,ifilty);

range = min([rangeEdge1', rangeEdge2'],[],2);
iforeground = find(range<fg_range);

iforeground = union(ifilt,iforeground);

if showPlots
    plot_foreground = -ones(length(recon_grid),1);
    plot_foreground(iforeground) = 1;
    imagesccyl(plot_foreground,recon_grid);
    drawnow;
end
%% Load Measurements, resize to res*res
disp("Loading measurements");
tic

if labData
    % Note: These images are 8-bit JPGs, and must be scaled accordingly
    % for another filetype
    switch(scnNum)
        case 0
            filepath = 'Laboratory Data\One Red Target\';
            newFrame1 = double(imread(strcat(filepath, 'Measurements\scn_0_exp_0_left.jpg')))/256;
            newFrame2 = double(imread(strcat(filepath, 'Measurements\scn_0_exp_0_right.jpg')))/256;
        case 1
            filepath = 'Laboratory Data\One Red Target\';
            newFrame1 = double(imread(strcat(filepath, 'Measurements\scn_1_exp_0_left.jpg')))/256;
            newFrame2 = double(imread(strcat(filepath, 'Measurements\scn_1_exp_0_right.jpg')))/256;
        case 2
            filepath = 'Laboratory Data\One Red Target\';
            newFrame1 = double(imread(strcat(filepath, 'Measurements\scn_2_exp_0_left.jpg')))/256;
            newFrame2 = double(imread(strcat(filepath, 'Measurements\scn_2_exp_0_right.jpg')))/256;
        case 3
            filepath = 'Laboratory Data\One Red Target\';
            newFrame1 = double(imread(strcat(filepath, 'Measurements\scn_3_exp_0_left.jpg')))/256;
            newFrame2 = double(imread(strcat(filepath, 'Measurements\scn_3_exp_0_right.jpg')))/256;
        case 4
            filepath = 'Laboratory Data\Two Targets\';
            newFrame1 = double(imread(strcat(filepath, 'Measurements\scn_1_exp_0_left.jpg')))/256;
            newFrame2 = double(imread(strcat(filepath, 'Measurements\scn_1_exp_0_right.jpg')))/256;
        case 5
            filepath = 'Laboratory Data\Two Targets\';
            newFrame1 = double(imread(strcat(filepath, 'Measurements\scn_5_exp_0_left.jpg')))/256;
            newFrame2 = double(imread(strcat(filepath, 'Measurements\scn_5_exp_0_right.jpg')))/256;
        case 6
            filepath = 'Laboratory Data\Three Targets\';
            newFrame1 = double(imread(strcat(filepath, 'Measurements\scn_1_exp_0_left.jpg')))/256;
            newFrame2 = double(imread(strcat(filepath, 'Measurements\scn_1_exp_0_right.jpg')))/256;
        otherwise
            disp("Data does not exist")
            return
    end
    bg1 = double(imread(strcat(filepath, 'Measurements\scn_bg_exp_0_left.jpg')))/256;
    bg2 = double(imread(strcat(filepath, 'Measurements\scn_bg_exp_0_right.jpg')))/256;
    defSize = 0.15; % size in meters of the original meas
else
    % Note: These images are 16-bit PNGs, and must be scaled accordingly
    % for another filetype
    switch(scnNum)
        case 0
            filepath = 'Synthetic Data\Emissive Example\';
            newFrame1 = double(imread(strcat(filepath, 'Measurements\colorful_object_left_meas.png')))/2^16;
            newFrame2 = double(imread(strcat(filepath, 'Measurements\colorful_object_right_meas.png')))/2^16;
            bg1 = double(imread(strcat(filepath, 'Measurements\colorful_background_left_meas.png')))/2^16;
            bg2 = double(imread(strcat(filepath, 'Measurements\colorful_background_right_meas.png')))/2^16;
        case 1
            filepath = 'Synthetic Data\Two Target Example\';
            newFrame1 = double(imread(strcat(filepath, 'Measurements\frame_0_samps_',num2str(samps),'_left_meas.png')))/2^16;
            newFrame2 = double(imread(strcat(filepath, 'Measurements\frame_0_samps_',num2str(samps),'_right_meas.png')))/2^16;
            bg1 = double(imread(strcat(filepath, 'Measurements\background_samps_',num2str(8000),'_left_meas.png')))/2^16;
            bg2 = double(imread(strcat(filepath, 'Measurements\background_samps_',num2str(8000),'_right_meas.png')))/2^16;
            if (fov == 0.5)
                newFrame1 = double(imread(strcat(filepath, 'Measurements\frame_0_samps_',num2str(samps),'_left_meas_large.png')))/2^16;
                newFrame2 = double(imread(strcat(filepath, 'Measurements\frame_0_samps_',num2str(samps),'_right_meas_large.png')))/2^16;
                bg1 = double(imread(strcat(filepath, 'Measurements\frame_bg_samps_',num2str(8000),'_left_meas_large.png')))/2^16;
                bg2 = double(imread(strcat(filepath, 'Measurements\frame_bg_samps_',num2str(8000),'_right_meas_large.png')))/2^16;
            end
            
        case 2
            filepath = 'Synthetic Data\Four Target Example\';
            newFrame1 = double(imread(strcat(filepath, 'Measurements\frame_0_samps_',num2str(samps),'_left_meas.png')))/2^16;
            newFrame2 = double(imread(strcat(filepath, 'Measurements\frame_0_samps_',num2str(samps),'_right_meas.png')))/2^16;
            bg1 = double(imread(strcat(filepath, 'Measurements\background_samps_',num2str(8000),'_left_meas.png')))/2^16;
            bg2 = double(imread(strcat(filepath, 'Measurements\background_samps_',num2str(8000),'_right_meas.png')))/2^16;
        case 3
            filepath = 'Synthetic Data\Same Color Example\';
            newFrame1 = double(imread(strcat(filepath, 'Measurements\frame_0_samps_',num2str(samps),'_left_meas.png')))/2^16;
            newFrame2 = double(imread(strcat(filepath, 'Measurements\frame_0_samps_',num2str(samps),'_right_meas.png')))/2^16;
            bg1 = double(imread(strcat(filepath, 'Measurements\background_samps_',num2str(8000),'_left_meas.png')))/2^16;
            bg2 = double(imread(strcat(filepath, 'Measurements\background_samps_',num2str(8000),'_right_meas.png')))/2^16;
        case 4
            filepath = 'Synthetic Data\Non-Lambertian Example\';
            newFrame1 = double(imread(strcat(filepath, 'Measurements\frame_0_samps_',num2str(samps),'_left_meas.png')))/2^16;
            newFrame2 = double(imread(strcat(filepath, 'Measurements\frame_0_samps_',num2str(samps),'_right_meas.png')))/2^16;
            bg1 = double(imread(strcat(filepath, 'Measurements\background_samps_',num2str(8000),'_left_meas.png')))/2^16;
            bg2 = double(imread(strcat(filepath, 'Measurements\background_samps_',num2str(8000),'_right_meas.png')))/2^16;
        otherwise
            disp("Data does not exist")
            return
    end
    defSize = 0.3; % size in meters of the original meas
    if (fov == 0.5)
        defSize = 0.5;
    end
end

if ~bgRemoval
    bg1 = zeros(size(bg1));
    bg2 = zeros(size(bg2));
end

defRes1 = size(bg1); % Resolution of original image
defRes2 = size(bg2); % Resolution of original image

cropX1 = 1:round(defRes1(1)*p.meas_size/defSize);
cropY1 = round((defRes1(2)/2)-(defRes1(2)*p.meas_size/defSize/2)+1) : round((defRes1(2)/2)+(defRes1(2)*p.meas_size/defSize/2));
cropX2 = 1:round(defRes1(1)*p.meas_size/defSize);
cropY2 = round((defRes2(2)/2)-(defRes2(2)*p.meas_size/defSize/2)+1) : round((defRes2(2)/2)+(defRes2(2)*p.meas_size/defSize/2));

% bg1 = myimresize(bg1(cropX1, cropY1, :),[p.meas_pix,p.meas_pix]);
bg1 = imresize(bg1(cropX1, cropY1, :),[p.meas_pix,p.meas_pix],"Antialiasing",true);
% bg2 = myimresize(bg2(cropX2, cropY2, :),[p.meas_pix,p.meas_pix]);
bg2 = imresize(bg2(cropX2, cropY2, :),[p.meas_pix,p.meas_pix],"Antialiasing",true);
% newFrame1 = myimresize(newFrame1(cropX1, cropY1, :),[p.meas_pix,p.meas_pix]);
newFrame1 = imresize(newFrame1(cropX1, cropY1, :),[p.meas_pix,p.meas_pix],"Antialiasing",true);
% newFrame2 = myimresize(newFrame2(cropX2, cropY2, :),[p.meas_pix,p.meas_pix]);
newFrame2 = imresize(newFrame2(cropX2, cropY2, :),[p.meas_pix,p.meas_pix],"Antialiasing",true);

meas1 = newFrame1 - bg1; meas2 = newFrame2 - bg2;

toc

%% Save images (optional)
if showPlots

    if labData
        savefilepath = [filepath '/Figures/scn_' num2str(scnNum) '_fov_' num2str(fov) '_res_' num2str(res) '/'];
    else 
        savefilepath = [filepath '/Figures/scn_' num2str(scnNum) '_samps_' num2str(samps) '_fov_' num2str(fov) '_res_' num2str(res) '/'];
    end
    mkdir(savefilepath);

    figure
    image(newFrame1)
    axis square; axis off;
    set(gcf,'Color',[1 1 1]); set(gca,'Color',[0 0 0]); set(gcf,'InvertHardCopy','off'); set(gca,'FontSize',10); 
    fig = gcf; fig.PaperUnits = 'centimeters'; fig.Units = 'centimeters'; 
    fig.PaperPosition = .5*[-1.2 -1.2 10 10]; fig.PaperSize=.5*[8 7.9];
    
    print(fig,[savefilepath 'newFrame_left.pdf'],'-dpdf');
    
    figure
    image(newFrame2)
    axis square; axis off;
    set(gcf,'Color',[1 1 1]); set(gca,'Color',[0 0 0]); set(gcf,'InvertHardCopy','off'); set(gca,'FontSize',10); 
    fig = gcf; fig.PaperUnits = 'centimeters'; fig.Units = 'centimeters'; 
    fig.PaperPosition = .5*[-1.2 -1.2 10 10]; fig.PaperSize=.5*[8 7.9];
    
    print(fig,[savefilepath 'newFrame_right.pdf'],'-dpdf');

    figure
    image(bg1)
    axis square; axis off;
    set(gcf,'Color',[1 1 1]); set(gca,'Color',[0 0 0]); set(gcf,'InvertHardCopy','off'); set(gca,'FontSize',10); 
    fig = gcf; fig.PaperUnits = 'centimeters'; fig.Units = 'centimeters'; 
    fig.PaperPosition = .5*[-1.2 -1.2 10 10]; fig.PaperSize=.5*[8 7.9];
    
    print(fig,[savefilepath 'background_left.pdf'],'-dpdf');

    figure
    image(bg2)
    axis square; axis off;
    set(gcf,'Color',[1 1 1]); set(gca,'Color',[0 0 0]); set(gcf,'InvertHardCopy','off'); set(gca,'FontSize',10); 
    fig = gcf; fig.PaperUnits = 'centimeters'; fig.Units = 'centimeters'; 
    fig.PaperPosition = .5*[-1.2 -1.2 10 10]; fig.PaperSize=.5*[8 7.9];
    
    print(fig,[savefilepath 'background_right.pdf'],'-dpdf');

    dispMeas1 = meas1./max(abs(meas1(:)))/2+0.5;
    dispMeas2 = meas2./max(abs(meas2(:)))/2+0.5;
    figure
    image(dispMeas1)
    axis square; axis off;
    set(gcf,'Color',[1 1 1]); set(gca,'Color',[0 0 0]); set(gcf,'InvertHardCopy','off'); set(gca,'FontSize',10); 
    fig = gcf; fig.PaperUnits = 'centimeters'; fig.Units = 'centimeters'; 
    fig.PaperPosition = .5*[-1.2 -1.2 10 10]; fig.PaperSize=.5*[8 7.9];
    
    print(fig,[savefilepath 'meas_left.pdf'],'-dpdf');

    figure
    image(dispMeas2)
    axis square; axis off;
    set(gcf,'Color',[1 1 1]); set(gca,'Color',[0 0 0]); set(gcf,'InvertHardCopy','off'); set(gca,'FontSize',10); 
    fig = gcf; fig.PaperUnits = 'centimeters'; fig.Units = 'centimeters'; 
    fig.PaperPosition = .5*[-1.2 -1.2 10 10]; fig.PaperSize=.5*[8 7.9];
    
    print(fig,[savefilepath 'meas_right.pdf'],'-dpdf');

end

%% Apply Mask
disp("Discarding corner pixels and first row");
tic
exclude_radius = res/16;
mask = zeros(p.meas_pix);
for i = 1:exclude_radius
    mask(i,:) = sqrt(([1:p.meas_pix]-p.meas_pix/2).^2+i^2) <= exclude_radius;
end
mask = ~mask(:);

meas1=meas1.*reshape(mask,p.meas_pix,p.meas_pix);
meas2=meas2.*reshape(mask,p.meas_pix,p.meas_pix);
m1 = reshape(meas1,p.meas_pix*p.meas_pix,3);
m2 = reshape(meas2,p.meas_pix*p.meas_pix,3);

A_ext = A_ext.*[mask;mask]; % pointwise multiplication of every column
toc
% rows are zeroed, not removed! -sheila

%% Sparse regularized recovery with scaled A
% lambdaGroup = 0 ;
lambdas = [lambdaGroup];

rec_params.nIter = 20000;
rec_params.stepSize = 1e-4;
if (res > 64)
    rec_params.stepSize = 3e-5;
elseif (res < 32)
    rec_params.stepSize = 1e-3;
end
rec_params.plotFreq = 1000;
rec_params.nAmb = 1;
rec_params.renderNum = scnNum;
rec_params.rec_size = p.rec_size(1);
rec_params.meas_size = p.meas_size;
rec_params.meas_pix = p.meas_pix;
rec_params.room_dim = p.room_dim(1);
rec_params.samps = samps;
rec_params.fov = fov;
rec_params.height = facetHeight;

toggles.plotStuff = showPlots;
toggles.saveStuff = 1;
toggles.debias = 0;
toggles.realData = labData;

tic
disp("Getting Sparse Recovery...")
recover(A_ext,[m1;m2],(1/sqrt(3))*ones(size(A_ext,2),3),.1*ones(size(A_ext,2),2),recon_grid,lambdas,rec_params,toggles,iforeground,filepath);
% savefig(2,strcat('scn-',num2str(renderNum),'-',num2str(meas_pix),'p-',num2str(rec_size(1)),'r-',num2str(round(meas_size*100)),'.fig'));
% savefig(2,strcat('grp-',num2str(lambdaGroup),'-slf-', num2str(lambdaSwagSelf), '-cpl-', num2str(lambdaSwagCouple),'-scn-',num2str(renderNum),'.fig'));

toc
end
%% Functions
function recover(A,meas,initV,initC,recon_grid,lambdas,rec_params,toggles,iforeground,filepath)

A_left = A(1:end/2,:);
A_right = A(end/2+1:end,:);
A_left_T = A_left';
A_right_T = A_right';

meas_left = meas(1:end/2,:);
meas_right = meas(end/2+1:end,:);
 
lambdaGroup = lambdas(1);

plotStuff = toggles.plotStuff;
saveStuff = toggles.saveStuff;
realData = toggles.realData;

nIter = rec_params.nIter;
stepSize = rec_params.stepSize;
plotFreq = rec_params.plotFreq;
nAmb = rec_params.nAmb;
renderNum = rec_params.renderNum;
rec_size = rec_params.rec_size;
meas_size = rec_params.meas_size;
meas_pix = rec_params.meas_pix;
room_dim = rec_params.room_dim;
samps = rec_params.samps;
fov = rec_params.fov;
height = rec_params.height;
tic

y = [initV,initC];
x_old = y;
x_new = y;
t_old = 1;

% meas = meas(:);

disp('starting iterations')



% init costs
cost = 1;
cost_data = 0;
cost_groupSparse = 0;

count = 0;

for it = 1:nIter
    

    alpha = (vecnorm(y(:,1:3)'))';
    if (sum(isnan(alpha))>0 || sum(isinf(alpha))>0)
        disp('check issue alpha!')
    end

    if (sum(isnan(y(:)))>0 || sum(isinf(y(:)))>0)
        disp('check issue y!')
    end

    rgb = y(:,1:3)./repmat(alpha,[1,3]);
    c1 = y(:,4).*alpha;
    c2 = y(:,5).*alpha;

    % % % % % % % Differentiable components % % % % % % %
    % L2 component of gradient
    % derivative wrt rgb variables
    res_left = (A_left*(repmat(c1,[1,3]).*rgb) - meas_left); % residual left   
    res_right = (A_right*(repmat(c2,[1,3]).*rgb) - meas_right); % residual left
    gradV = [diag(c1)*A_left_T, diag(c2)*A_right_T]*[res_left; res_right];
    
    % derivative wrt c variables
    gradC1 = [diag(rgb(:,1))*A_left_T, diag(rgb(:,2))*A_left_T, diag(rgb(:,3))*A_left_T]...
        *res_left(:);
    
    gradC2 = [diag(rgb(:,1))*A_right_T, diag(rgb(:,2))*A_right_T, diag(rgb(:,3))*A_right_T]...
        *res_right(:);
    
    grad = [gradV, gradC1, gradC2];
    
    
    % MEGS/SWAGGER component of gradient
    
%     if it>megsStart && lambdaMegs~=0
% 
%         lambdaMegs = lambdaMegs*megsGrow;
%         v_megs = [x_new(:,4);x_new(:,5)];
%         grad_megs = lambdaMegs*2*S*v_megs;
%         grad(:,4:5) = grad(:,4:5) + [grad_megs(1:end/2),grad_megs(end/2+1:end)];
%        
%     end
    
    % Update x
    x_new = y - stepSize * grad;
    
    % keep RGB between 0 and 1
%     x_new(:,1:3) = min(x_new(:,1:3),1);
    x_new(:,1:3) = max(x_new(:,1:3),0);
    
    % % % % % % % Non-differentialble components  % % % % % % %
    % group lasso over the C's 
 
    
    v_group = x_new(1:end-1,4:5);
    v_group_l2 =  sqrt(sum(v_group.^2,2));
    
    i_group_zero = find(v_group_l2<lambdaGroup); % indices to shrink to zero
    v_group_new = v_group;
    v_group_new(i_group_zero,:) = 0;
    v_group_l2(i_group_zero,:) = 1; % avoid divide by zeros
    v_group_new = (1-(lambdaGroup)./v_group_l2 ).*v_group_new;
    v_group_new(iforeground,:) = max(0, v_group_new(iforeground,:));
    x_new(1:end-1,4:5) = v_group_new;
    
%     if it>megsStart && lambdaMegs~=0
%         cm = lambdaMegs*v_megs'*S*v_megs;
%     end
%     if it>(megsStart+2*plotFreq) && lambdaMegs~=0
%        v_megs = [x_new(:,4);x_new(:,5)];
%        
%        if cm<1e-8
%            saveStuff = 0;
%            lambdaMegs = 0;
%            
%        end
%        
%     end
% 
    if sum(isnan(x_new(:)))>0
        stopHere = 1;
    end

    
    t_new = 0.5*(1 + sqrt(1 + 4*t_old^2));

    % update y
    y = x_new + (t_old - 1)/t_new*(x_new - x_old);
    x_old = x_new;
    
    t_old = t_new;
    
    if mod(it,plotFreq) == 0
        temp = x_new(:,4:5);
        
        toc
        tic
        count = count+1;
        
        % % % % % % % Calculate Cost % % % % % % %
        cost_data(count) =  (1/2)*norm( [res_left(:); res_right(:)] )^2;
        v_group = y(1:end-1,4:5);
        
        cost_groupSparse(count) = lambdaGroup*sum( sqrt( sum( v_group.^2,2) ));
        

        cost(count) =   cost_data(count) + cost_groupSparse(count);

        
        % % % % % % % Save progress % % % % % % %
        if saveStuff == 1
            if (realData)
                filename = strcat('recon-grp-',num2str(lambdaGroup), ...
                    '-scn-',num2str(renderNum),...
                    '-ht-',num2str(height),...
                    '-fov-',num2str(fov),...
                    '-res-',num2str(meas_pix),...
                    '.mat');
            else
                filename = strcat('recon-grp-',num2str(lambdaGroup), ...
                    '-scn-',num2str(renderNum),...
                    '-samps-',num2str(samps),...
                    '-ht-',num2str(height),...
                    '-fov-',num2str(fov),...
                    '-res-',num2str(meas_pix),...
                    '.mat');
            end
            savefilepath = [filepath 'data_files/' filename];
            mkdir([filepath 'data_files/']);
            save(savefilepath,'plotFreq','stepSize','lambdaGroup','cost',...
                'cost_groupSparse','cost_data','y','recon_grid','A_left','A_right','meas')
        end
        
        % % % % % % % Plot progress % % % % % % %
        if plotStuff == 1
            c1_zero = zeros(size(c1));
            c1_zero(c1==0) = 1;
            disp(['n zero in c1 = ' num2str(sum(c1_zero))])
            
            it
            if (mod(it,plotFreq*5)==0)
                clf
            end
            c1 = y(:,4);
            c2 = y(:,5);
            [c1_pos, c1_neg] = plotPosNeg(c1);
            [c2_pos, c2_neg] = plotPosNeg(c2);
            
            if norm((c1_pos-c1_neg)-c1)>1e-5
                disp('Code issue. Check!')
            end
            disp_left = diag(c1)*y(:,1:3);
            disp_left_pos = diag(c1_pos)*y(:,1:3);
            disp_left_neg = diag(c1_neg)*y(:,1:3);
            
            disp_right = diag(c2)*y(:,1:3);
            disp_right_pos = diag(c2_pos)*y(:,1:3);
            disp_right_neg = diag(c2_neg)*y(:,1:3);
            
            figure(2)
            
            % Plot costs
            subplot(331)
            semilogy(cost)
            title('Cost')
            xlabel(strcat('Iteration (x',num2str(plotFreq),')'))
            axis square
            grid on
            
            subplot(332)
            semilogy(cost_data)
            title('Cost Data')
            xlabel(strcat('Iteration (x',num2str(plotFreq),')'))
            axis square
            grid on
            
            subplot(333)
            semilogy(cost_groupSparse)
            title('Cost group L1')
            xlabel(strcat('Iteration (x',num2str(plotFreq),')'))
            axis square
            grid on
            

            % Plot scene
            perc_sat = .025; % percent to saturate i.e. greater than 1 values
            scale = max([c1; c2]);
            sat_threshold = prctile_one(reshape([disp_left;disp_right],[],1),1-perc_sat);
%             scaleFact = 1*max(max([disp_left_pos(:); disp_right_pos(:)]));
%             scaleFactNeg = 1*max(max([disp_left_neg(:); disp_right_neg(:)]));
%             scaleFactCombo = max([scaleFact, scaleFactNeg]);

            
            
            
            subplot(345)
            imagecyl(disp_left_pos./max(max(disp_left_pos)),recon_grid);
            title('RGB*C1 pos');
            axis square; set(gca,'color','black'); xlim([0,room_dim]); ylim([0,room_dim]); axis on;
            
            subplot(346)
            imagecyl(disp_left_neg./max(max(disp_left_pos)),recon_grid);
            title('RGB*C1 neg');
            axis square; set(gca,'color','black'); xlim([0,room_dim]); ylim([0,room_dim]); axis on;
            
            subplot(3,4,7)
            imagecyl(c1_pos./scale,recon_grid);
            title('c1_pos');
            axis square; set(gca,'color','black'); xlim([0,room_dim]); ylim([0,room_dim]); axis on;
            
            subplot(3,4,8)
            imagecyl(c1_neg./scale,recon_grid);
            title('c1_neg');
            axis square; set(gca,'color','black'); xlim([0,room_dim]); ylim([0,room_dim]); axis on;
            

            subplot(349)
            imagecyl(disp_right_pos./max(max(disp_right_pos)),recon_grid);
            title('RGB*C2 pos');
            axis square; set(gca,'color','black'); xlim([0,room_dim]); ylim([0,room_dim]); axis on;
            
            subplot(3,4,10)
            imagecyl(disp_right_neg./max(max(disp_right_pos)),recon_grid);
            title('RGB*C2 neg');
            axis square; set(gca,'color','black'); xlim([0,room_dim]); ylim([0,room_dim]); axis on;
            
            subplot(3,4,11)
            imagecyl(c2_pos./scale,recon_grid);
            title('c2_pos');
            axis square; set(gca,'color','black'); xlim([0,room_dim]); ylim([0,room_dim]); axis on;
            
            subplot(3,4,12)
            imagecyl(c2_neg./scale,recon_grid);
            title('c2_neg');
            axis square; set(gca,'color','black'); xlim([0,room_dim]); ylim([0,room_dim]); axis on;
            
            
            
            
            disp_recon = [A_left*disp_left;A_right*disp_right];
            
            
            figure(4)
            subplot(3,3,1)
            imagesc(reshape(meas(:,1),meas_pix,meas_pix*2))
            title('Meas GT Red')
            colorbar; axis off;
            
            subplot(3,3,2)
            imagesc(reshape(meas(:,2),meas_pix,meas_pix*2))
            title('Meas GT Green')
            colorbar; axis off;
            
            subplot(3,3,3)
            imagesc(reshape(meas(:,3),meas_pix,meas_pix*2))
            title('Meas GT Blue')
            colorbar; axis off;
            
            subplot(3,3,4)
            imagesc(reshape(disp_recon(:,1),meas_pix,meas_pix*2))
            title('Recon Red')
            colorbar; axis off;
            
            subplot(3,3,5)
            imagesc(reshape(disp_recon(:,2),meas_pix,meas_pix*2))
            title('Recon Green')
            colorbar; axis off;
            
            subplot(3,3,6)
            imagesc(reshape(disp_recon(:,3),meas_pix,meas_pix*2))
            title('Recon Blue')
            colorbar; axis off;
            
            subplot(3,3,7)
            imagesc(reshape(meas(:,1)-disp_recon(:,1),meas_pix,meas_pix*2))
            title('Meas-recon Red')
            colorbar; axis off;
            
            subplot(3,3,8)
            imagesc(reshape(meas(:,2)-disp_recon(:,2),meas_pix,meas_pix*2))
            title('Meas-recon Green')
            colorbar; axis off;
            
            subplot(3,3,9)
            imagesc(reshape(meas(:,3)-disp_recon(:,3),meas_pix,meas_pix*2))
            title('Meas-recon Blue')
            colorbar; axis off;
            
            drawnow
        end
    end
end

end

function V_x = prctile_one(V,p)
V = sort(V,'ascend');
N = length(V);
x = p*(N-1)+1; % position x
if floor(x) < N
    V_x = V(floor(x)) + mod(x,1)*(V(floor(x)+1) - V(floor(x))); % value
else
    V_x = V(N); % position N
end
end