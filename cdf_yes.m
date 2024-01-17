clear;clc;
OBS=xlsread('F:\cmip6小论文0915\CMCC\0 情景数据\CMCC585data.xlsx','历史观测值'); %1995-2014
SAT=xlsread('F:\cmip6小论文0915\CMCC\0 情景数据\CMCC585data.xlsx','历史模式数据');%1995-2014
SAT_F=xlsread('F:\cmip6小论文0915\CMCC\0 情景数据\CMCC585data.xlsx','未来模式数据','A1:B72'); %2015-2020
OBS_F=xlsread('F:\cmip6小论文0915\CMCC\0 情景数据\CMCC585data.xlsx','未来观测值');   %2015-2020

D_OBS = OBS(:,1);  %读取第一列全部  观测值对应的时间
OBSdata=(OBS(:,2));%读取第二列全部  观测值对应的数值
D_SAT = SAT(:,1);   %读取第一列全部 模拟值对应的时间
SATdata=(SAT(:,2)); %读取第二列全部 模拟值对应的数值
D_SAT_F = SAT_F(:,1);   %读取第一列全部 未来模拟值对应的时间
SATdata_F=(SAT_F(:,2)); %读取第二列全部 未来模拟值对应的数值
SATCDF_ok = NaN(length(SATdata_F),1);  %创建一个length*1的二维矩阵
D_OBS_F = OBS_F(:,1);  %读取第一列全部  观测值对应的时间
OBSdata_F=(OBS_F(:,2));%读取第二列全部  观测值对应的数值
%COE_cdf=[]

for i=1:12     % for each month of the year
    OBSmont = OBSdata(find(month(datetime(datevec(D_OBS)))==i));  %读出不同年份在某一月的观测值
    OBSmont_F = OBSdata_F(find(month(datetime(datevec(D_OBS_F)))==i));
    ID_SATmont = find(month(datetime(datevec(D_SAT)))==i);        
    SATmont = SATdata(ID_SATmont);                                %读出不同年份在某一月的模拟值
    ID_SATmont_F = find(month(datetime(datevec(D_SAT_F)))==i); 
    SATmont_F = SATdata_F(ID_SATmont_F);                                %读出不同年份在某一月的模拟值
    
    POBS = [1:length(OBSmont)]'./(length(OBSmont)+1);             %计算观测值的累积分布函数
    PSAT = [1:length(SATmont)]'./(length(SATmont)+1);             %计算模拟值的累积分布函数
    
    % Interpolation of input data (SAT) to the same percentiles of
    % 将模拟数据内插到同观测数据相同百分位数
    % benchmark data (OBS)
    SATint= interp1(PSAT,sort(SATmont),sort(POBS),'linear','extrap');  %将模式值从小到大排序
    %D=sort(OBSmont);
    % Computation of the differences between the CDFs of OBS and SAT data
    DIFF=sort(OBSmont)-SATint;  %两组数据排序后，观测-模拟 作差
    
    % Fitting of a polynomial curve to DIFF
    COEFF= polyfit(SATint,DIFF, 2);   %p = polyfit(x,y,n) 返回次数为 n 的多项式 p(x) 的系数
    %COE_cdf=[COE_cdf COEFF]
    %coe_monthly=reshape(COE_cdf,6,12)
    
    % Evaluation of the polynomial curve to SAT data 
    SATCDF= polyval(COEFF,SATmont_F)+SATmont_F;  %y = polyval(p,x) 计算多项式 p 在 x 的每个点处的值
    SATCDF_ok(ID_SATmont_F) = SATCDF;
    
    % Comparison of cdf curves estimated for: benchmark data (OBSmont), 
    % data to modify (SATmont), corrected data (SATCDF)
    set(gcf,'position',[ 530, 190, 1111, 794])
    subplot(3,4,i)
    plot( sort(OBSmont_F),(1:length(OBSmont_F))/(length(OBSmont_F)+1),'Color',0.7*[1,1,1], 'linewidth',7)
    hold on
    plot(sort(SATmont_F),(1:length(SATmont_F))/(length(SATmont_F)+1), 'b-','linewidth',4)  %b-表示蓝色实线
    plot( sort(SATCDF),(1:length(SATCDF))/(length(SATCDF)+1), 'r--', 'linewidth',2)  %r--表示红色虚线
    xlabel('data'), ylabel('Cumulative Density Function')  %x轴是数据值，y轴是累积分布
    monthTitle =['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
    title(monthTitle(i,:),'fontweight','bold','fontsize',10), grid on
    if i==12, legend ('Reference data','Original biased data','Corrected data',4), end
    
    
    M_STAT_OBS(i,1)= nanmean(OBSmont); V_STAT_OBS(i,1) = nanvar(OBSmont);
    M_STAT_ST(i,1)= nanmean(SATCDF);   V_STAT_ST(i,1) = nanvar(SATCDF);   %剔除Nan值后求均值和方差
    
end
