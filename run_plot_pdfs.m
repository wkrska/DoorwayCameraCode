clc
clear all
close all
addpath("Utils\");

%% --------- User Parameters --------- %%
labData = 1;
scnNum = 6;
samp = 8000;
facetHeight = 1;
fov = .15;
res = 64;
lambdaGroup = 1e-6;
bgRemoval = 1;
t_val = 0.005;

%% --------- Room Parameters --------- %%
door_width = .508;
door_center = 1.5;
room_dim = [3,3]; % Size of the room in meters

%% --------- Load Data File --------- %%
if labData
    switch(scnNum)
        case 0
            dirpath = 'Laboratory Data\One Red Target\';
        case 1
            dirpath = 'Laboratory Data\One Red Target\';
        case 2
            dirpath = 'Laboratory Data\One Red Target\';
        case 3
            dirpath = 'Laboratory Data\One Red Target\';
        case 4
            dirpath = 'Laboratory Data\Two Targets\';
        case 5
            dirpath = 'Laboratory Data\Two Targets\';
        case 6
            dirpath = 'Laboratory Data\Three Targets\';
        otherwise
            disp("Data does not exist")
            return
    end
else
    switch(scnNum)
        case 0
            dirpath = 'Synthetic Data\Emissive Example\';
        case 1
            dirpath = 'Synthetic Data\Two Target Example\';
        case 2
            dirpath = 'Synthetic Data\Four Target Example\';
        case 3
            dirpath = 'Synthetic Data\Same Color Example\';
        case 4
            dirpath = 'Synthetic Data\Non-Lambertian Example\';
        otherwise
            disp("Data does not exist")
            return
    end
end
    
if labData
    load(strcat(dirpath,'data_files/recon-grp-',num2str(lambdaGroup), ...
                    '-scn-',num2str(scnNum),...
                    '-ht-',num2str(facetHeight),...
                    '-fov-',num2str(fov),...
                    '-res-',num2str(res),...
                    '.mat'));
    savefilepath = [dirpath '/Figures/scn_' num2str(scnNum) '_fov_' num2str(fov) '_res_' num2str(res) '/grp_' num2str(lambdaGroup) '_ht_' num2str(facetHeight)];
else
    load(strcat(dirpath,'data_files/recon-grp-',num2str(lambdaGroup), ...
                    '-scn-',num2str(scnNum),...
                    '-samps-',num2str(samp),...
                    '-ht-',num2str(facetHeight),...
                    '-fov-',num2str(fov),...
                    '-res-',num2str(res),...
                    '.mat'));
    savefilepath = [dirpath '/Figures/scn_' num2str(scnNum) '_samps_' num2str(samp) '_fov_' num2str(fov) '_res_' num2str(res) '/grp_' num2str(lambdaGroup) '_ht_' num2str(facetHeight)];
end

c1 = y(:,4);
c2 = y(:,5);
[c1_pos, c1_neg] = plotPosNeg(c1);
[c2_pos, c2_neg] = plotPosNeg(c2);

disp_left = diag(c1)*y(:,1:3);
disp_left_pos = diag(c1_pos)*y(:,1:3);
disp_left_neg = diag(c1_neg)*y(:,1:3);

disp_right = diag(c2)*y(:,1:3);
disp_right_pos = diag(c2_pos)*y(:,1:3);
disp_right_neg = diag(c2_neg)*y(:,1:3);
    
scaleFact = 1*max(max([disp_left_pos(:); disp_right_pos(:)]));
scaleFactNeg = 1*max(max([disp_left_neg(:); disp_right_neg(:)]));
scaleFactCombo = max([scaleFact, scaleFactNeg]);
    
disp_combo = plotCombo(disp_left,disp_right,c1,c2,t_val);
disp_combo_neg = plotCombo(-disp_left,-disp_right,-c1,-c2,t_val);

%% --------- Plot Combo --------- %%
figure
disp("Combo")
imagecyl(disp_combo./scaleFactCombo,recon_grid);
hold on
scatter(door_center-door_width/2, 0,[],[1 1 0],'filled')
scatter(door_center+door_width/2, 0,[],[1 1 0],'filled')
axis square; set(gca,'color','black'); xlim([0,room_dim(1)]); ylim([0,room_dim(2)]); axis on;
set(gcf,'Color',[1 1 1]); set(gca,'Color',[0 0 0]); set(gcf,'InvertHardCopy','off');
set(gca,'FontSize',10)
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.Units = 'centimeters';
fig.PaperPosition = .5*[-1.2 -1.2 10 10]; 
fig.PaperSize=.5*[8 7.9];

print(fig,[savefilepath '_recon_combo.pdf'],'-dpdf');
   
%% --------- Plot Left --------- %%

figure
disp("left pos")
imagecyl(disp_left_pos./scaleFactCombo,recon_grid);
hold on
scatter(door_center-door_width/2, 0,[],[1 1 0],'filled')
scatter(door_center+door_width/2, 0,[],[1 1 0],'filled')
axis square; set(gca,'color','black'); xlim([0,room_dim(1)]); ylim([0,room_dim(2)]); axis on;
set(gcf,'Color',[1 1 1]); set(gca,'Color',[0 0 0]); set(gcf,'InvertHardCopy','off');
set(gca,'FontSize',10)
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.Units = 'centimeters';
fig.PaperPosition = .5*[-1.2 -1.2 10 10]; 
fig.PaperSize=.5*[8 7.9];
    
print(fig,[savefilepath '_recon_left.pdf'],'-dpdf');

%% --------- Plot Right --------- %%

figure
disp("right pos")
imagecyl(disp_right_pos./scaleFactCombo,recon_grid);
hold on
scatter(door_center-door_width/2, 0,[],[1 1 0],'filled')
scatter(door_center+door_width/2, 0,[],[1 1 0],'filled')
axis square; set(gca,'color','black'); xlim([0,room_dim(1)]); ylim([0,room_dim(2)]); axis on;
set(gcf,'Color',[1 1 1]); set(gca,'Color',[0 0 0]); set(gcf,'InvertHardCopy','off');
set(gca,'FontSize',10)

fig = gcf;
fig.PaperUnits = 'centimeters';
fig.Units = 'centimeters';
fig.PaperPosition = .5*[-1.2 -1.2 10 10]; 
fig.PaperSize=.5*[8 7.9];
    
print(fig,[savefilepath '_recon_right.pdf'],'-dpdf');

%% --------- Plot Combo Negative--------- %%
figure
disp("Combo Neg")
imagecyl(disp_combo_neg./scaleFactCombo,recon_grid);
hold on
scatter(door_center-door_width/2, 0,[],[1 1 0],'filled')
scatter(door_center+door_width/2, 0,[],[1 1 0],'filled')
axis square; set(gca,'color','black'); xlim([0,room_dim(1)]); ylim([0,room_dim(2)]); axis on;
set(gcf,'Color',[1 1 1]); set(gca,'Color',[0 0 0]); set(gcf,'InvertHardCopy','off');
set(gca,'FontSize',10)
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.Units = 'centimeters';
fig.PaperPosition = .5*[-1.2 -1.2 10 10]; 
fig.PaperSize=.5*[8 7.9];

print(fig,[savefilepath '_recon_combo_neg.pdf'],'-dpdf');
   
%% --------- Plot Left Negative --------- %%
figure
disp("left neg")
imagecyl(disp_left_neg./scaleFactCombo,recon_grid);
    hold on
scatter(door_center-door_width/2, 0,[],[1 1 0],'filled')
scatter(door_center+door_width/2, 0,[],[1 1 0],'filled')
axis square; set(gca,'color','black'); xlim([0,room_dim(1)]); ylim([0,room_dim(2)]); axis on;
set(gcf,'Color',[1 1 1]); set(gca,'Color',[0 0 0]); set(gcf,'InvertHardCopy','off');
set(gca,'FontSize',10)
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.Units = 'centimeters';
fig.PaperPosition = .5*[-1.2 -1.2 10 10]; 
fig.PaperSize=.5*[8 7.9];
    
print(fig,[savefilepath '_recon_left_neg.pdf'],'-dpdf');

%% --------- Plot Right Negative --------- %%
figure
disp("right neg")
imagecyl(disp_right_neg./scaleFactCombo,recon_grid);
    hold on
scatter(door_center-door_width/2, 0,[],[1 1 0],'filled')
scatter(door_center+door_width/2, 0,[],[1 1 0],'filled')
axis square; set(gca,'color','black'); xlim([0,room_dim(1)]); ylim([0,room_dim(2)]); axis on;
set(gcf,'Color',[1 1 1]); set(gca,'Color',[0 0 0]); set(gcf,'InvertHardCopy','off');
set(gca,'FontSize',10)

fig = gcf;
fig.PaperUnits = 'centimeters';
fig.Units = 'centimeters';
fig.PaperPosition = .5*[-1.2 -1.2 10 10]; 
fig.PaperSize=.5*[8 7.9];

print(fig,[savefilepath '_recon_right_neg.pdf'],'-dpdf');

%%
% end 
% end