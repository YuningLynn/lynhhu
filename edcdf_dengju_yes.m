clear;clc;
OBS_H=xlsread('F:\cmip6С����0915\CMCC\0 �龰����\CMCC585data.xlsx','��ʷ�۲�ֵ');  %1995-2014
OBS_Hmin=xlsread('F:\cmip6С����0915\���ļ������\5.EDCDFУ��\SSP126\BCC\history_pre.xlsx','1980-2020');%�˲�����Ϊ���ʶ�ĳС�����µĽ�ˮֵ 1980-2020
SAT_H=xlsread('F:\cmip6С����0915\CMCC\0 �龰����\CMCC585data.xlsx','��ʷģʽ����'); %1995-2014
SAT_F=xlsread('F:\cmip6С����0915\CMCC\0 �龰����\CMCC585data.xlsx','δ��ģʽ����'); %2015-2034
OBS_F=xlsread('F:\cmip6С����0915\CMCC\0 �龰����\CMCC585data.xlsx','δ���۲�ֵ');   %2015-2020

D_OBS_H = OBS_H(:,1);  %��ȡ��һ��ȫ��  ��ʷ�۲�ֵ��Ӧ��ʱ��
OBSdata_H=(OBS_H(:,2));%��ȡ�ڶ���ȫ��  ��ʷ�۲�ֵ��Ӧ����ֵ
D_OBS_Hmin = OBS_Hmin(:,1);  %��ȡ��һ��ȫ��  ��ʷ�۲�ֵ��Ӧ��ʱ��
OBSdata_Hmin=(OBS_Hmin(:,2));%��ȡ�ڶ���ȫ��  ��ʷ�۲�ֵ��Ӧ����ֵ

% loading data to modify
D_SAT_F = SAT_F(:,1);   %��ȡ��һ��ȫ�� δ��ģ��ֵ��Ӧ��ʱ��
SATdata_F=(SAT_F(:,2)); %��ȡ�ڶ���ȫ�� δ��ģ��ֵ��Ӧ����ֵ
SATEDCDF_ok = NaN(length(SATdata_F),1);  %����һ��length*1�Ķ�ά����

%������ʷ��ģʽģ������
D_SAT_H = SAT_H(:,1);   %��ȡ��һ��ȫ�� ��ʷģ��ֵ��Ӧ��ʱ��
SATdata_H=(SAT_H(:,2)); %��ȡ�ڶ���ȫ�� ��ʷģ��ֵ��Ӧ����ֵ

D_OBS_F = OBS_F(:,1);  %��ȡ��һ��ȫ��  �۲�ֵ��Ӧ��ʱ��
OBSdata_F=(OBS_F(:,2));%��ȡ�ڶ���ȫ��  �۲�ֵ��Ӧ����ֵ

D=[];%�������d1
MF=[];%����������У��ֵ
MO=[];%������������Сֵ
% Mobs=[];%�������1995-2014���µĹ۲�ֵ
% Mmod=[];%�������1995-2014���µĹ۲�ֵ
sat_monthlyF=[];
for i=1:12
    OBSmont_H = OBSdata_H(find(month(datetime(datevec(D_OBS_H)))==i));  %������ͬ�����ĳһ�µĹ۲�ֵ
    OBSmont_F = OBSdata_F(find(month(datetime(datevec(D_OBS_F)))==i));
    OBSmont_Hmin = OBSdata_Hmin(find(month(datetime(datevec(D_OBS_Hmin)))==i));  %������ͬ�����ĳһ�µĹ۲�ֵ
    
    ID_SATmont_F = find(month(datetime(datevec(D_SAT_F)))==i); 
    SATmont_F = SATdata_F(ID_SATmont_F);                                %������ͬ�����ĳһ�µ�ģ��ֵ
    ID_SATmont_H = find(month(datetime(datevec(D_SAT_H)))==i);        
    SATmont_H = SATdata_H(ID_SATmont_H);                                %������ͬ�����ĳһ�µ�ģ��ֵ
    
    
    POBS_H = [1:length(OBSmont_H)]'./(length(OBSmont_H)+1);           %����۲�ֵ���ۻ��ֲ�����
    POBS_F = [1:length(OBSmont_F)]'./(length(OBSmont_F)+1);
    PSAT_F = [1:length(SATmont_F)]'./(length(SATmont_F)+1);             %����ģ��ֵ���ۻ��ֲ�����
    PSAT_H = [1:length(SATmont_H)]'./(length(SATmont_H)+1);
    POBS_Hmin = [1:length(OBSmont_Hmin)]'./(length(OBSmont_Hmin)+1);           %����۲�ֵ���ۻ��ֲ�����
    
    OBS_A=interp1(POBS_Hmin,sort(OBSmont_Hmin),sort(POBS_H),'pchip','extrap');  %��ʷ�۲�
    SAT_B= interp1(POBS_H,sort(OBSmont_H),sort(PSAT_F),'pchip','extrap');  %δ���۲�
    SAT_C= interp1(PSAT_H,sort(SATmont_H),sort(PSAT_F),'pchip','extrap');  %δ��ģʽ\
    
    x = [1:1:length(SATmont_F)];
    y = SATmont_F;
    z = SAT_C;
    [y, ID_SATmont_F] =  sort (y);
    x = x(ID_SATmont_F);
    zz= z';
    zzz=zz(:,x);
    zzzz=zzz';
    
    r=SAT_B;
    rr = r';
    rrr=rr(:,x);
    rrrr=rrr';
    d1= rrrr-zzzz;
    D=[D d1];
    SAT_A=SATmont_F+d1;
    SATEDCDF_ini=SAT_A;   
    %................�˶��ǽ�����������ֵ������ֵ��Ȼ���滻Ϊ��ϵ���������Сֵ
    MO=[MO sort(OBS_A)];
    min_ob=MO(1,:)';
    for j=1:20
    if SATEDCDF_ini(j)<0
       SATEDCDF_ini(j)=min_ob(i);
    end
       SATEDCDF = SATEDCDF_ini;
    end
    sat_monthlyF = [sat_monthlyF SATmont_F];
%     %......................................................................
%     Mobs=[Mobs OBSmont_H];%�������1995-2014���µĹ۲�ֵ
%     Mmod=[Mmod SATmont_H];%�������1995-2014���µĹ۲�ֵ       �������Ƿ��µ�����ֵ�����ڷ��·ݵķ�λ��ӳ��
%   %......................................................................
    SATEDCDF_ok(ID_SATmont_F) = SATEDCDF;
    MF = [MF SATEDCDF]
    %SAT_COR=reshape(MF',1,240) %����Ҫ��δ��������
     
    set(gcf,'position',[ 530, 190, 1111, 794])
    subplot(3,4,i)
    plot( sort(OBSmont_F),(1:length(OBSmont_F))/(length(OBSmont_F)+1),'Color',0.7*[1,1,1], 'linewidth',7) %��֪��δ���۲�����
    hold on
    plot(sort(SATmont_F),(1:length(SATmont_F))/(length(SATmont_F)+1), 'b-','linewidth',4)  %b-��ʾ��ɫʵ�ߣ�δ��ģʽ����
    plot( sort(SATEDCDF),(1:length(SATEDCDF))/(length(SATEDCDF)+1), 'r--', 'linewidth',2)  %r--��ʾ��ɫ���ߣ��������ģʽ����
    %plot( sort(SAT_A),(1:length(SAT_A))/(length(SAT_A)+1), 'o-', 'linewidth',2)  %r--��ʾ��ɫ����
    xlabel('Pre(mm)'), ylabel('Cumulative Density Function')  %x��������ֵ��y�����ۻ��ֲ�
    monthTitle =['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
    title(monthTitle(i,:),'fontweight','bold','fontsize',10), grid on
    if i==12, legend ('Reference data','Original biased data','Corrected data',4), end
end