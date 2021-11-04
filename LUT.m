% _______________________________________________________________________
% LUT.m
% Look-Up Table
% retrieve chlorophyll by PROSPECT
% _______________________________________________________________________

clear all%清楚工作空间的所有变量、函数和MEX文件
clc%清楚命令窗口的内容，对工作环境中的全部变量无任何影响
addpath PROSPECT-D_Matlab
 
% %%
% N = 1.4;
% Cw = 0.09;
% Anth = 2;
% Cbrown = 0;
% Cab = 0.00001;
% Cm = 0.00001;
% Car = Cab*0.25;
% reflectance = [400:1:2500];
% transmittance = [400:1:2500];
% CT = table(N, Cab, Car, Anth, Cbrown, Cw, Cm, reflectance, transmittance); %创建表
% 
% % Cabs =[1:0.5:100];
% % Cms = [0.002:0.005:0.05];
% % Cws = [0.01:0.005:0.09];
% 
% for Cab = 1:0.1:100
%     for Cm = 0.002:0.005:0.05
%          %for Cw = 0.01:0.005:0.09
%             fprintf('current processing Cab is %d\n',Cab);
%             Car = Cab*0.25;
%             LRT=prospect_DB(N, Cab, Car, Anth, Cbrown, Cw, Cm);     %PROSPECT
%             pig_ref_trans = {N, Cab, Car, Anth, Cbrown, Cw, Cm, LRT(:,2).', LRT(:,3).'};   
%             CT = [CT; pig_ref_trans];  
%          %end
%     end   
% end   
% writetable(CT, 'CHL_LookUpTable2.csv', 'WriteRowNames', true);

%%
RT = shaperead('AlakesR.shp') ; 
CT = readtable('CHL_LookUpTable.csv'); 

%%
chl = [];
tic

for i = 1:1:size(RT, 1)
    
    fprintf('current processing index is %d\n',i);
    
    % image_ref = [RT(i).BR0469 RT{i, 'BR0555'} RT{i, 'BR0645'} RT{i, 'BR0859'} RT{i, 'BR1240'} RT{i, 'BR1640'} RT{i, 'BR2130'}]*0.001;
    image_ref = [RT(i).BR0469 RT(i).BR0555 RT(i).BR0645 RT(i).BR0859 RT(i).BR1240 RT(i).BR1640 RT(i).BR2130]*0.001;
    image_refs = repmat( image_ref, size(CT, 1), 1);
    
    simula_refs = [CT{:, 77} CT{:, 163} CT{:, 253} CT{:, 467} CT{:, 848} CT{:, 1248} CT{:, 1738}];
    
    discre = sum(power(image_refs-simula_refs, 2), 2);%寻找最小光谱差异的数据集

    [m, p] = min(discre);
    RT(i).chl = CT{p, 'Cab'};
    % chl = [chl, CT{p, 'Cab'}];
    disp(m);
    disp(p);
    toc
end

%%
shapewrite(RT,'AlakesR1909_LUT.shp');
