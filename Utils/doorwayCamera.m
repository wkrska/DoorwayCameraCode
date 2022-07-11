function doorwayCamera(renderNum, bgStr, frameStr, lambdaGroup,dirpath,nIter,stepSize,plotFreq)
%% Init

dirpath = [dirpath 'Figures/Frame' num2str(renderNum) '/'];
close all

global meas_pix

    
meas_pix = 64; % meas_pix x meas_pix image at both corners
meas_size = 0.075; % FOV size in meters
rec_size = [45,45]; % Number of discrete angles per edge
room_dim = [3,3]; % Size of the room in meters

p.rec_size = rec_size;
p.meas_pix = meas_pix;
p.meas_size = meas_size;
p.room_dim = room_dim;
p.min_angle = pi/8; % this was previously 10, which exceeds the max (pi). This gives you a 0 length vector for sweeping through the angles

door_width = .508;
door_center = 1.5;%room_dim(2)/2;
pos = [door_center-door_width/2, -0.002; door_center+door_width/2, -0.002];

%% Load forward model for each corner
tic
filename = strcat(num2str(meas_pix),'p-',num2str(room_dim(1)),'m-',num2str(rec_size(1)),'r-',num2str(round(meas_size*100)),'cm-min-pi_',num2str(pi/p.min_angle),'.mat');

disp(strcat("Loading A mat: ",filename));
try
    load(filename);
catch
    disp("No data file for specified resolution, generating A matrix");
    tic
    [A, recon_grid, S] = makeA3D_cylindrical_plane_modified_sheila(pos,1,16,p,1);
    toc
    save(filename, 'p','meas_pix','meas_size','rec_size','room_dim','A','recon_grid','S');
end



nAmb = 0;
A_ext = [A, ones(size(A,1),1)];


toc
%% figure out which pixels can only be positive (i.e. S_{FG} in paper)

fg_range = 1.1; % everything inside of this range must be positive
center = [door_center,0];
rangeEdge1 = sqrt((squeeze(recon_grid(1,:,1))-center(1)).^2 + (squeeze(recon_grid(2,:,1))-center(2)).^2);
rangeEdge2 = sqrt((squeeze(recon_grid(1,:,2))-center(1)).^2 + (squeeze(recon_grid(2,:,2))-center(2)).^2);
y_coord = min([squeeze(recon_grid(2,:,1))',squeeze(recon_grid(2,:,2))'],[],2);
x_coord = min([squeeze(abs(recon_grid(1,:,1)-center(1)))',squeeze(abs(recon_grid(1,:,2)-center(1)))'],[],2);
ifiltx = find(x_coord<1);
ifilty = find(y_coord<1.85);
ifilt = intersect(ifiltx,ifilty);

range = min([rangeEdge1', rangeEdge2'],[],2);
iforeground = find(range<fg_range);

iforeground = union(ifilt,iforeground);

testScene = zeros(230,3);
testScene(:,1) = .5;
testScene(:,2) = .2;
testScene(:,3) = .5;

testScene(iforeground,1) = 1;
testScene(iforeground,2) = 0;
testScene(iforeground,3) = 0;

figure
hold on
imagecylGraph(testScene./(max(max(testScene))),recon_grid);
axis square; 
axis on;
xlabel('x [m]','Interpreter','latex','FontSize',14)
ylabel('y [m]','Interpreter','latex','FontSize',14)
scatter(pos(:,1),pos(:,2),'filled','b')
xlim([-.1 3.1])
ylim([0 3])

fig = gcf; 
fig.PaperUnits = 'centimeters';  
fig.PaperPosition = [0 0 10 10]; 
fig.Units = 'centimeters'; 
fig.PaperSize=[9.5 9.5]; 
fig.Units = 'centimeters'; 
print(fig,['biangularReconstructionGrid.pdf'],'-dpdf','-r200'); 

%% Load Measurements, resize to 256x256
disp("Loading measurements");
tic

newFrame1 = double(imread([frameStr '_left_meas.png']))/2^16;
newFrame2 = double(imread([frameStr '_right_meas.png']))/2^16;
bg1 = double(imread([bgStr '_left_meas.png']))/2^16;
bg2 = double(imread([bgStr '_right_meas.png']))/2^16;


meas1 = newFrame1 - bg1; meas2 = newFrame2 - bg2;

defRes1 = size(meas1); % Resolution of original image
defRes2 = size(meas2); % Resolution of original image
defSize = 0.15; % size in meters of the original meas
cropX1 = 1:round(defRes1(1)*meas_size/defSize);
cropY1 = round((defRes1(2)/2)-(defRes1(2)*meas_size/defSize/2)+1) : round((defRes1(2)/2)+(defRes1(2)*meas_size/defSize/2));
cropX2 = 1:round(defRes1(1)*meas_size/defSize);
cropY2 = round((defRes2(2)/2)-(defRes2(2)*meas_size/defSize/2)+1) : round((defRes2(2)/2)+(defRes2(2)*meas_size/defSize/2));

meas1 = imresize(meas1(cropX1, cropY1, :),[meas_pix,meas_pix],"Antialiasing",true);
meas2 = imresize(meas2(cropX2, cropY2, :),[meas_pix,meas_pix],"Antialiasing",true);

toc
%% Apply Mask
disp("Discarding corner pixels and first row");
tic
dead_radius = 4;
mask = zeros(p.meas_pix);
for i = 1:dead_radius
    mask(i,:) = sqrt(([1:p.meas_pix]-p.meas_pix/2).^2+i^2) < dead_radius;
end
mask = ~mask(:);

meas1=meas1.*reshape(mask,meas_pix,meas_pix);
meas2=meas2.*reshape(mask,meas_pix,meas_pix);
m1 = reshape(meas1,meas_pix*meas_pix,3);
m2 = reshape(meas2,meas_pix*meas_pix,3);

A_ext = A_ext.*[mask;mask]; % pointwise multiplication of every column
toc
% rows are zeroed, not removed! -sheila

%% Plot Measurments

figure(5)
subplot(3,3,1)
imagesc(bg1(:,:,1))
title('R bg')
colorbar; axis off;


subplot(3,3,2)
imagesc(newFrame1(:,:,1))
title('R new frame')
colorbar; axis off;


subplot(3,3,3)
imagesc(meas1(:,:,1))
title('R diff')
colorbar; axis off;


subplot(3,3,4)
imagesc(bg1(:,:,2))
title('G bg')
colorbar; axis off;


subplot(3,3,5)
imagesc(newFrame1(:,:,2))
title('G new frame')
colorbar; axis off;


subplot(3,3,6)
imagesc(meas1(:,:,2))
title('G diff')
colorbar; axis off;


subplot(3,3,7)
imagesc(bg1(:,:,3))
title('B bg')
colorbar; axis off;


subplot(3,3,8)
imagesc(newFrame1(:,:,3))
title('B new frame')
colorbar; axis off;


subplot(3,3,9)
imagesc(meas1(:,:,3))
title('B diff')
colorbar; axis off;

fig = gcf;
fig.PaperUnits = 'centimeters';
fig.Units = 'centimeters';
fig.PaperPosition = 3*[0 0 8 6];
fig.PaperSize=3*[8 6];

print(fig,[dirpath '/Measurements_' num2str(renderNum) '.pdf'],'-dpdf','-r200');
   



figure(1091)
subplot(231)
imshow(bg1)
title('bg')

subplot(232)
imshow(newFrame1)
title('new frame')

subplot(233)
imshow( (meas1-min(meas1(:))) / (max(meas1(:)) - min(meas1(:)) ))
title('amped dff')

subplot(234)
imshow(bg2)
title('bg')

subplot(235)
imshow(newFrame2)
title('new frame')

subplot(236)
imshow( (meas2-min(meas2(:))) / (max(meas2(:)) - min(meas2(:)) ))
title('amplified dff')

fig = gcf;
fig.PaperUnits = 'centimeters';
fig.Units = 'centimeters';
fig.PaperPosition = 3*[0 0 8 6];
fig.PaperSize=3*[8 6];

print(fig,[dirpath '/MeasurementsRGB_' num2str(renderNum) '.pdf'],'-dpdf','-r200');
   


plotDiff1 = newFrame1 - bg1; 
plotDiff1 = (plotDiff1-min(plotDiff1(:)))/(max(plotDiff1(:))-min(plotDiff1(:)));
plotDiff2 = newFrame2 - bg2;
plotDiff2 = (plotDiff2-min(plotDiff2(:)))/(max(plotDiff2(:))-min(plotDiff2(:)));

close all
figure(1)
imshow(plotDiff1)
axis off
axis square
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.Units = 'centimeters';
fig.PaperPosition = [-.135 -.17 1 1];
fig.PaperSize=[.725 .74];
print(fig,[dirpath '/leftAmpedDiffMeas_' num2str(renderNum) '.pdf'],'-dpdf','-r200');
    
figure(2)
imshow( plotDiff2)
axis off
axis square
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.Units = 'centimeters';
fig.PaperPosition = [-.135 -.17 1 1];
fig.PaperSize=[.725 .74];

print(fig,[dirpath '/rightAmpedDiffMeas_' num2str(renderNum) '.pdf'],'-dpdf','-r200');
    

    
figure(3)
imshow( newFrame1)
axis off
axis square
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.Units = 'centimeters';
fig.PaperPosition = [-.135 -.17 1 1];
fig.PaperSize=[.725 .74];

print(fig,[dirpath '/newFrameLeft_' num2str(renderNum) '.pdf'],'-dpdf','-r200');
    
figure(4)
imshow( newFrame2)
axis off
axis square
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.Units = 'centimeters';
fig.PaperPosition = [-.135 -.17 1 1];
fig.PaperSize=[.725 .74];

print(fig,[dirpath '/newFrameRight_' num2str(renderNum) '.pdf'],'-dpdf','-r200');
    

    
figure(5)
imshow( bg1)
axis off
axis square
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.Units = 'centimeters';
fig.PaperPosition = [-.135 -.17 1 1];
fig.PaperSize=[.725 .74];

print(fig,[dirpath '/bgLeft_' num2str(renderNum) '.pdf'],'-dpdf','-r200');
    
figure(6)
imshow( bg2)
axis off
axis square
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.Units = 'centimeters';
fig.PaperPosition = [-.135 -.17 1 1];
fig.PaperSize=[.725 .74];
print(fig,[dirpath '/bgRight_' num2str(renderNum) '.pdf'],'-dpdf','-r200');
 

%% Sparse regularized recovery with scaled A
% lambdaGroup = 0 ;
lambdas = lambdaGroup;

% nIter = 12000;
% stepSize = 1e-4;
% plotFreq = 1000;
% megsStart = 3000;
% megsGrow = 1.001;
params = [nIter,stepSize,plotFreq,nAmb,renderNum,rec_size(1),meas_size,room_dim(1)];

plotStuff = 1;
saveStuff = 1;
toggles = [plotStuff,saveStuff];
S = sparse(S);
tic
disp("Getting Sparse Recovery...")
recover(A_ext,S,[m1;m2],(1/sqrt(3))*ones(size(A_ext,2),3),.1*ones(size(A_ext,2),2),recon_grid,lambdas,params,toggles,iforeground);
% savefig(2,strcat('scn-',num2str(renderNum),'-',num2str(meas_pix),'p-',num2str(rec_size(1)),'r-',num2str(round(meas_size*100)),'.fig'));
% savefig(2,strcat('grp-',num2str(lambdaGroup),'-slf-', num2str(lambdaSwagSelf), '-cpl-', num2str(lambdaSwagCouple),'-scn-',num2str(renderNum),'.fig'));

toc
end
%% Functions
function recover(A,S,meas,initV,initC,recon_grid,lambdas,params,toggles,iforeground)
global meas_pix

A_left = A(1:end/2,:);
A_right = A(end/2+1:end,:);
A_left_T = A_left';
A_right_T = A_right';




meas_left = meas(1:end/2,:);
meas_right = meas(end/2+1:end,:);
 


lambdaGroup = lambdas(1);

plotStuff = toggles(1);
saveStuff = toggles(2);


nIter = params(1); 
stepSize = params(2);
plotFreq = params(3);
renderNum = params(5);
room_dim=params(8);

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
cost_couple = 0;

count = 0;
initial = 1;

for it = 1:nIter
    

    alpha = (1./vecnorm(y(:,1:3)'))';
    if sum(isnan(alpha))>0
        disp('check issue!')
    end

    rgb = y(:,1:3).*repmat(alpha,[1,3]);
    c1 = y(:,4)./alpha;
    c2 = y(:,5)./alpha;
    

%     rgb = y(:,1:3);
%     c1 = y(:,4);
%     c2 = y(:,5);
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
            filename = strcat('grp-',num2str(lambdaGroup), '-scn-',num2str(renderNum),'.mat');
            save(filename,'plotFreq','stepSize','lambdaGroup','cost',...
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
%             sat_threshold = prctile_one(reshape([disp_left;disp_right],[],1),1-perc_sat);
            

            scaleFact = 1*max(max([disp_left_pos(:); disp_right_pos(:)]));
            scaleFactNeg = 1*max(max([disp_left_neg(:); disp_right_neg(:)]));
            scaleFactCombo = max([scaleFact, scaleFactNeg]);

            subplot(345)
            imagecyl(disp_left_pos./scaleFact,recon_grid);
            title('RGB*C1 pos');
            axis square; set(gca,'color','black'); xlim([0,room_dim]); ylim([0,room_dim]); axis on;
            
            subplot(346)
            imagecyl(disp_left_neg./scaleFactNeg,recon_grid);
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
            imagecyl(disp_right_pos./scaleFact,recon_grid);
            title('RGB*C2 pos');
            axis square; set(gca,'color','black'); xlim([0,room_dim]); ylim([0,room_dim]); axis on;
            
            subplot(3,4,10)
            imagecyl(disp_right_neg./scaleFactNeg,recon_grid);
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
%         
%         % Check if function is runaway
%         if cost(count) > cost(initial)*10 || max(reshape(x_new,[],1)) > 100
%             it = 1;
%             cost = 0;
%             y = [initV,initC];
% 
%             x_old = y;
%             x_new = y;
%             t_old = 1;
% 
%             initial = count + 1;
%             stepSize = stepSize/10;
%             disp(strcat("Step size reduced to: ", num2str(stepSize)));
%         end
    end
end

% debias
% if debias == 1
%     notzero = find(x_new_left(:,1) + x_new_left(:,2) + x_new_left(:,3) ~= 0);
%     Anew = A(:,notzero);
%     debiased = cat(3, lsqnonneg(Anew, meas(:,1)), lsqnonneg(Anew, meas(:,2)) , lsqnonneg(Anew, meas(:,3)));
%     x_new_left(notzero,:) = debiased;
%     notzero = find(x_new_left(:,1) + x_new_left(:,2) + x_new_left(:,3) ~= 0);
%     Anew = A(:,notzero);
%     debiased = cat(3, lsqnonneg(Anew, meas(:,1)), lsqnonneg(Anew, meas(:,2)) , lsqnonneg(Anew, meas(:,3)));
%     x_new_right(notzero,:) = debiased;
% end

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