clear;clc;
OBS=xlsread('F:\cmip6С����0915\CMCC\0 �龰����\CMCC585data.xlsx','��ʷ�۲�ֵ'); %1995-2014
SAT=xlsread('F:\cmip6С����0915\CMCC\0 �龰����\CMCC585data.xlsx','��ʷģʽ����');%1995-2014
SAT_F=xlsread('F:\cmip6С����0915\CMCC\0 �龰����\CMCC585data.xlsx','δ��ģʽ����','A1:B72'); %2015-2020
OBS_F=xlsread('F:\cmip6С����0915\CMCC\0 �龰����\CMCC585data.xlsx','δ���۲�ֵ');   %2015-2020

D_OBS = OBS(:,1);  %��ȡ��һ��ȫ��  �۲�ֵ��Ӧ��ʱ��
OBSdata=(OBS(:,2));%��ȡ�ڶ���ȫ��  �۲�ֵ��Ӧ����ֵ
D_SAT = SAT(:,1);   %��ȡ��һ��ȫ�� ģ��ֵ��Ӧ��ʱ��
SATdata=(SAT(:,2)); %��ȡ�ڶ���ȫ�� ģ��ֵ��Ӧ����ֵ
D_SAT_F = SAT_F(:,1);   %��ȡ��һ��ȫ�� δ��ģ��ֵ��Ӧ��ʱ��
SATdata_F=(SAT_F(:,2)); %��ȡ�ڶ���ȫ�� δ��ģ��ֵ��Ӧ����ֵ
SATCDF_ok = NaN(length(SATdata_F),1);  %����һ��length*1�Ķ�ά����
D_OBS_F = OBS_F(:,1);  %��ȡ��һ��ȫ��  �۲�ֵ��Ӧ��ʱ��
OBSdata_F=(OBS_F(:,2));%��ȡ�ڶ���ȫ��  �۲�ֵ��Ӧ����ֵ
%COE_cdf=[]

for i=1:12     % for each month of the year
    OBSmont = OBSdata(find(month(datetime(datevec(D_OBS)))==i));  %������ͬ�����ĳһ�µĹ۲�ֵ
    OBSmont_F = OBSdata_F(find(month(datetime(datevec(D_OBS_F)))==i));
    ID_SATmont = find(month(datetime(datevec(D_SAT)))==i);        
    SATmont = SATdata(ID_SATmont);                                %������ͬ�����ĳһ�µ�ģ��ֵ
    ID_SATmont_F = find(month(datetime(datevec(D_SAT_F)))==i); 
    SATmont_F = SATdata_F(ID_SATmont_F);                                %������ͬ�����ĳһ�µ�ģ��ֵ
    
    POBS = [1:length(OBSmont)]'./(length(OBSmont)+1);             %����۲�ֵ���ۻ��ֲ�����
    PSAT = [1:length(SATmont)]'./(length(SATmont)+1);             %����ģ��ֵ���ۻ��ֲ�����
    
    % Interpolation of input data (SAT) to the same percentiles of
    % ��ģ�������ڲ嵽ͬ�۲�������ͬ�ٷ�λ��
    % benchmark data (OBS)
    SATint= interp1(PSAT,sort(SATmont),sort(POBS),'linear','extrap');  %��ģʽֵ��С��������
    %D=sort(OBSmont);
    % Computation of the differences between the CDFs of OBS and SAT data
    DIFF=sort(OBSmont)-SATint;  %������������󣬹۲�-ģ�� ����
    
    % Fitting of a polynomial curve to DIFF
    COEFF= polyfit(SATint,DIFF, 2);   %p = polyfit(x,y,n) ���ش���Ϊ n �Ķ���ʽ p(x) ��ϵ��
    %COE_cdf=[COE_cdf COEFF]
    %coe_monthly=reshape(COE_cdf,6,12)
    
    % Evaluation of the polynomial curve to SAT data 
    SATCDF= polyval(COEFF,SATmont_F)+SATmont_F;  %y = polyval(p,x) �������ʽ p �� x ��ÿ���㴦��ֵ
    SATCDF_ok(ID_SATmont_F) = SATCDF;
    
    % Comparison of cdf curves estimated for: benchmark data (OBSmont), 
    % data to modify (SATmont), corrected data (SATCDF)
    set(gcf,'position',[ 530, 190, 1111, 794])
    subplot(3,4,i)
    plot( sort(OBSmont_F),(1:length(OBSmont_F))/(length(OBSmont_F)+1),'Color',0.7*[1,1,1], 'linewidth',7)
    hold on
    plot(sort(SATmont_F),(1:length(SATmont_F))/(length(SATmont_F)+1), 'b-','linewidth',4)  %b-��ʾ��ɫʵ��
    plot( sort(SATCDF),(1:length(SATCDF))/(length(SATCDF)+1), 'r--', 'linewidth',2)  %r--��ʾ��ɫ����
    xlabel('data'), ylabel('Cumulative Density Function')  %x��������ֵ��y�����ۻ��ֲ�
    monthTitle =['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
    title(monthTitle(i,:),'fontweight','bold','fontsize',10), grid on
    if i==12, legend ('Reference data','Original biased data','Corrected data',4), end
    
    
    M_STAT_OBS(i,1)= nanmean(OBSmont); V_STAT_OBS(i,1) = nanvar(OBSmont);
    M_STAT_ST(i,1)= nanmean(SATCDF);   V_STAT_ST(i,1) = nanvar(SATCDF);   %�޳�Nanֵ�����ֵ�ͷ���
    
end
