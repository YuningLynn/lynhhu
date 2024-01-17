clear;clc;
OBS_H=xlsread('F:\cmip6小论文0915\CMCC\0 情景数据\CMCC585data.xlsx','历史观测值');  %1995-2014
OBS_Hmin=xlsread('F:\cmip6小论文0915\论文计算过程\5.EDCDF校正\SSP126\BCC\history_pre.xlsx','1980-2020');%此步骤是为了率定某小概率下的降水值 1980=2020
SAT_H=xlsread('F:\cmip6小论文0915\CMCC\0 情景数据\CMCC585data.xlsx','历史模式数据'); %1995-2014
SAT_F=xlsread('F:\cmip6小论文0915\CMCC\0 情景数据\CMCC585data.xlsx','未来模式数据'); %2015-2034
OBS_F=xlsread('F:\cmip6小论文0915\CMCC\0 情景数据\CMCC585data.xlsx','未来观测值');   %2015-2020

D_OBS_H = OBS_H(:,1);  %读取第一列全部  历史观测值对应的时间
OBSdata_H=(OBS_H(:,2));%读取第二列全部  历史观测值对应的数值
D_OBS_Hmin = OBS_Hmin(:,1);  %读取第一列全部  历史观测值对应的时间
OBSdata_Hmin=(OBS_Hmin(:,2));%读取第二列全部  历史观测值对应的数值

% loading data to modify
D_SAT_F = SAT_F(:,1);   %读取第一列全部 未来模拟值对应的时间
SATdata_F=(SAT_F(:,2)); %读取第二列全部 未来模拟值对应的数值
SATEDCDF_ok = NaN(length(SATdata_F),1);  %创建一个length*1的二维矩阵

%载入历史的模式模拟数据
D_SAT_H = SAT_H(:,1);   %读取第一列全部 历史模拟值对应的时间
SATdata_H=(SAT_H(:,2)); %读取第二列全部 历史模拟值对应的数值

D_OBS_F = OBS_F(:,1);  %读取第一列全部  观测值对应的时间
OBSdata_F=(OBS_F(:,2));%读取第二列全部  观测值对应的数值

%ID_SATmont_F=[];
%ID_SATmont_H=[];
MF=[];%用来存各年的校正值
MO=[];%用来存各年的最小值
Mobs=[];%用来存放1995-2014各月的观测值
Mmod=[];%用来存放1995-2014各月的观测值
% coe=[];%用来存每个月的系数
sat_monthlyF=[]
Obshis=[];
for i=1:12
    OBSmont_H = OBSdata_H(find(month(datetime(datevec(D_OBS_H)))==i));  %读出不同年份在某一月的观测值
    OBSmont_F = OBSdata_F(find(month(datetime(datevec(D_OBS_F)))==i));
    OBSmont_Hmin = OBSdata_Hmin(find(month(datetime(datevec(D_OBS_Hmin)))==i));  %读出不同年份在某一月的观测值
    
    ID_SATmont_F = find(month(datetime(datevec(D_SAT_F)))==i); 
    SATmont_F = SATdata_F(ID_SATmont_F);                                %读出不同年份在某一月的模拟值
    ID_SATmont_H = find(month(datetime(datevec(D_SAT_H)))==i);        
    SATmont_H = SATdata_H(ID_SATmont_H);                                %读出不同年份在某一月的模拟值
    
    
    POBS_H = [1:length(OBSmont_H)]'./(length(OBSmont_H)+1);           %计算观测值的累积分布函数
    POBS_F = [1:length(OBSmont_F)]'./(length(OBSmont_F)+1);
    PSAT_F = [1:length(SATmont_F)]'./(length(SATmont_F)+1);             %计算模拟值的累积分布函数
    PSAT_H = [1:length(SATmont_H)]'./(length(SATmont_H)+1);
    POBS_Hmin = [1:length(OBSmont_Hmin)]'./(length(OBSmont_Hmin)+1);           %计算观测值的累积分布函数
    
    OBS_A=interp1(POBS_Hmin,sort(OBSmont_Hmin),sort(POBS_H),'pchip','extrap');  %未来观测
    SAT_B= interp1(POBS_H,sort(OBSmont_H),sort(PSAT_F),'pchip','extrap');  %未来观测
    SAT_C= interp1(PSAT_H,sort(SATmont_H),sort(PSAT_F),'pchip','extrap');  %未来模式
    d1 = SAT_B-SAT_C;
    COEFF= polyfit(SAT_C,d1,2);   %p = polyfit(x,y,n) 返回次数为 n 的多项式 p(x) 的系数
%   coe = [coe COEFF]
    %monthly_coe = reshape(coe,3,12)
    SATEDCDF_ini= polyval(COEFF,SATmont_F)+SATmont_F;  %y = polyval(p,x) 计算多项式 p 在 x 的每个点处的值
    
    %................此段是将初步订正的值挑出负值，然后替换为长系列里面的最小值
    MO=[MO sort(OBS_A)];
    min_ob=MO(1,:)';
    for j=1:20
    if SATEDCDF_ini(j)<0
       SATEDCDF_ini(j)=min_ob(i)
    end
       SATEDCDF = SATEDCDF_ini
    end
%     %......................................................................
     Mobs=[Mobs OBSmont_H];%用来存放1995-2014各月的观测值
     Mmod=[Mmod SATmont_H];%用来存放1995-2014各月的观测值       这两行是分月的输入值，用于分月份的分位数映射
%   %......................................................................
    SATEDCDF_ok(ID_SATmont_F) = SATEDCDF;
    MF = [MF SATEDCDF]
    Obshis=[Obshis OBSmont_H];
    sat_monthlyF=[sat_monthlyF SATmont_F] 
    set(gcf,'position',[ 530, 190, 1111, 794])
    subplot(3,4,i)
    plot( sort(OBSmont_F),(1:length(OBSmont_F))/(length(OBSmont_F)+1),'Color',0.7*[1,1,1], 'linewidth',7) %已知的未来观测数据
    hold on
    plot(sort(SATmont_F),(1:length(SATmont_F))/(length(SATmont_F)+1), 'b-','linewidth',4)  %b-表示蓝色实线，未来模式数据
    plot( sort(SATEDCDF),(1:length(SATEDCDF))/(length(SATEDCDF)+1), 'r--', 'linewidth',2)  %r--表示红色虚线，订正后的模式数据
    %plot( sort(SAT_A),(1:length(SAT_A))/(length(SAT_A)+1), 'o-', 'linewidth',2)  %r--表示红色虚线
    xlabel('Pre(mm)'), ylabel('Cumulative Density Function')  %x轴是数据值，y轴是累积分布
    monthTitle =['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
    title(monthTitle(i,:),'fontweight','bold','fontsize',10), grid on
    if i==12, legend ('Reference data','Original biased data','Corrected data',4), end
end