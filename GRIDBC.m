clear;clc;
load('F:\cmip6小论文0915\CMCC\7 特定月份全流域偏差校正\MH.mat');%历史模式 1995-2014
load('F:\cmip6小论文0915\CMCC\7 特定月份全流域偏差校正\OHandF.mat'); %观测 1980-2020
load('F:\cmip6小论文0915\CMCC\7 特定月份全流域偏差校正\MF.mat'); %未来模式 2015-2034
load('F:\cmip6小论文0915\CMCC\7 特定月份全流域偏差校正\Opart.mat');%观测 1995-2014

POBS_H = [1:20]'./21;         
POBS_F = [1:6]'./7;
PSAT_F = [1:20]'./21;            
PSAT_H = [1:20]'./21;
POBS_Hmin = [1:41]'./42;
OBSmont_Hmin=reshape(OHandF,41,4784);
OBSmont_H=reshape(Opart,20,4784);
SATmont_H=reshape(MH,20,4784);
SATmont_F=reshape(MF,20,4784);

%% 订正前原始数据 出图
CMCC20157=SATmont_F(6,:); %提取7月的观测值
CMCCA = CMCC20157';
CMCCB= reshape(CMCCA,92,52);
CMCCC=CMCCB';
[mask,Q]=geotiffread('F:\huaihe\hh_matlab01.tif');%
mask=double(mask);
mask(mask==0)=nan;
JulyCMCC=CMCCC.*mask;
column=NaN(52,2); 
row=NaN(2,96);
A_column=[JulyCMCC column];
B_column=[column A_column];
A_row=[B_column;row];
B_row=[row;A_row];
cmcc20157=B_row;
h=imagesc(cmcc20157);
set(h,'alphadata',~isnan(cmcc20157));
xlabel('Pre(mm/month)','FontName','Times New Roman');
title('SSP5-8.5 CMCC ','fontweight','bold','FontName','Times New Roman','fontsize',24), grid on
text(5,4,'Pre-ave = 194.5','FontSize',18);
clear text
for i =1:56
    f(i)=nanmean(cmcc20157(:));
end
p=mean(f);

xlabel('Pre(mm/month)','FontName','Times New Roman','fontweight','bold','fontsize',20);
set(gca,'yticklabel',35.5:-1:30.0,'fontsize',20);
set(gca,'xticklabel',113.5:2:122,'fontsize',20);
colormap(flipud(parula))% 反转colorbar的颜色
%%仿写流域平均
%     OBS_A=interp1(POBS_Hmin,sort(OBSmont_Hmin),sort(POBS_H),'pchip','extrap');  %未来观测
%     SAT_B= interp1(POBS_H,sort(OBSmont_H),sort(PSAT_F),'pchip','extrap');  %未来观测
%     SAT_C= interp1(PSAT_H,sort(SATmont_H),sort(PSAT_F),'pchip','extrap');  %未来模式
  
%% 插值  
    OBS_A = NaN.*zeros(20,4784);
    for i = 1:4784
        if length(find(isnan(OBSmont_Hmin(:,i)))) >0
            continue
        else
            OBS_A(:,i)=interp1(POBS_Hmin,sort(OBSmont_Hmin(:,i)),sort(POBS_H),'pchip','extrap');
        end
    end
    
    SAT_B = NaN.*zeros(20,4784);
    for i = 1:4784
        if length(find(isnan(OBSmont_H(:,i)))) >0
            continue
        else
            SAT_B(:,i)=interp1(POBS_H,sort(OBSmont_H(:,i)),sort(PSAT_F),'pchip','extrap');
        end
    end
    
     SAT_C = NaN.*zeros(20,4784);
     for i = 1:4784
         if length(find(isnan(SATmont_H(:,i)))) >0
             continue
         else
             SAT_C(:,i)=interp1(PSAT_H,sort(SATmont_H(:,i)),sort(PSAT_F),'pchip','extrap');
         end
     end
%% 检查插值后是否还能出完整的流域图  以找到2014年及以后流域数据是完整的
% ZZZ=OBS_A(20,:);
% zzzz=reshape(ZZZ,92,52);
% zzzzz=zzzz';
% ppp=OBSmont_Hmin(35,:);
% pppp=reshape(ppp,92,52);
% ppppp=pppp';
% %  [mask,Q]=geotiffread('F:\huaihe\hh_matlab01.tif');%
% %   mask=double(mask);
% %   mask(mask==0)=nan;
% %C= ppppp.*mask;
% %h=imagesc(C);  检查至此 也没有问题 模式数据可出图
%% 确定相关系数及拟合
d1 = SAT_B-SAT_C;
coe=[];
load('F:\cmip6小论文0915\CMCC\7 特定月份全流域偏差校正\min_ob.mat');
for i = 1:4784
    if length(find(isnan(d1(:,i)))) == 41
        continue
    else
        COEFF= polyfit(SAT_C(:,i),d1(:,i),2);
        coe=[coe COEFF];
    end
end
%grid_coe = reshape(coe,3,4784);
load('F:\cmip6小论文0915\CMCC\7 特定月份全流域偏差校正\ssp370\grid_coe.mat');
grid_coe=grid_coe;
MO=[];
RESULT=[];

for i=1:4784
    SATEDCDF_ini= polyval(grid_coe(:,i),SATmont_F(:,i))+SATmont_F(:,i);
    %MO=[MO sort(OBS_A(:,i))];
    min_ob=min_ob;
    for j=1:20
    if SATEDCDF_ini(j)<0
       SATEDCDF_ini(j)=min_ob(i); %min_ob到第59列也是有值的
    end
       SATEDCDF = SATEDCDF_ini;
    end
    RESULT=[RESULT SATEDCDF];
end

%%  格点订正完成，将RESULT（20*4784）转为矩阵
A=RESULT';
B=A(:,6);%列数为年份数 2015年对应第一列1
C=reshape(B,92,52);
D=C'; %改年为2015年7月的淮河流域订正降水数据
clear A B C POBS_F POBS_H POBS_Hmin PSAT_F PSAT_H MF MH d1;
%% 导入淮河TIFF 出图
[mask,Q]=geotiffread('F:\huaihe\hh_matlab01.tif');%
mask=double(mask);
mask(mask==0)=nan;
p=imagesc(mask);    %mask大小没有问题
JulyMOD= D.*mask;
column=NaN(52,2); 
row=NaN(2,96);
A_column=[JulyMOD column];
B_column=[column A_column];
A_row=[B_column;row];
B_row=[row;A_row];
ssp585=B_row;
h=imagesc(ssp585);
set(h,'alphadata',~isnan(ssp585))
xlabel('Pre(mm/month)','FontName','Times New Roman','fontsize',20);
title('(f) SSP5-8.5','fontweight','bold','FontName','Times New Roman','fontsize',24), grid on
text(7,4,'Pre-ave = 172.3','FontSize',12);
clear text
%% 计算流域平均降水
for i =1:56
    f(i)=nanmean(ssp585(:));
end
p=mean(f);
%% 导出tiff 可至arcgis编辑
[mask,Q]=geotiffread('F:\huaihe\hh_matlab01.tif');%
mask=double(mask);
mask(mask==0)=nan;
Lon = min(Q.LongitudeLimits)+0.1/2:0.1:max(Q.LongitudeLimits)-0.1/2;
Lat = min(Q.LatitudeLimits)+0.1/2:0.1:max(Q.LatitudeLimits)-0.1/2;
filepath = 'F:\cmip6小论文0915\CMCC\7 特定月份全流域偏差校正\';     
R = georefcells(Q.LatitudeLimits,Q.LongitudeLimits,size(mask)); % 设置一个地理坐标
file_name = 'JulyOBS';
geotiffwrite('JulyOBS',JulyOBS, R);
%% 观测值出图 2017 2018年缺数据
OBS20167=OBSmont_Hmin(41,:); %提取7月的观测值
OBSA = OBS20167';
OBSB= reshape(OBSA,92,52);
OBSC=OBSB';
[mask,Q]=geotiffread('F:\huaihe\hh_matlab01.tif');%
mask=double(mask);
mask(mask==0)=nan;
JulyOBS=OBSC.*mask;
column=NaN(52,2); 
row=NaN(2,96);
A_column=[JulyOBS column];
B_column=[column A_column];
A_row=[B_column;row];
B_row=[row;A_row];
obs20167=B_row;
h=imagesc(obs20167);
set(h,'alphadata',~isnan(obs20167));
xlabel('Pre(mm/month)','FontName','Times New Roman','fontweight','bold','fontsize',20);
set(gca,'yticklabel',35.5:-1:30.0,'fontsize',20);
set(gca,'xticklabel',113.5:2:122,'fontsize',20);
title('Observed ','fontweight','bold','FontName','Times New Roman','fontsize',24), grid on
text(5,4,'Pre-ave = 296.7','FontSize',18);
clear text
for i =1:56
    f(i)=nanmean(obs20167(:));
end
p=mean(f);
colormap(flipud(parula))% 反转colorbar的颜色