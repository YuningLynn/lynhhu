clear;clc;
%%
ncols=128; %����
nrows=72;%����
%%������������%%
readFilePath='F:\1961-2020��0.5x0.5���� ����������\����4\*07.txt';  %1980-1994
readPathStr='F:\1961-2020��0.5x0.5���� ����������\����4\';
fileList=dir(readFilePath);
fileNum=length(fileList);%�ļ�����
for i=1:fileNum
    name=fileList(i).name;
    splitName=strsplit(name,'.txt');  %��.����ȡ.ǰ����ַ���
    varStr = splitName{1};
    fileName=strcat(readPathStr,name);%������ ���ǻ��������ļ�������·��
    fid=fopen(fileName,'r');
    FormatString=repmat('%f',1,ncols);
    Pre(:,:,i)=cell2mat(textscan(fid,FormatString,nrows,'HeaderLines',6));%��ȡ�й���Χ��05�������ˮ����  ����ǰ����
    Pre_new=Pre(35:47,81:100,i);%��������Ӧ������ ��ά����
    Pre_new(Pre_new<0)=nan;
    lon=[112:0.5:121.5]; %����
    lat=[30.5:0.5:36.5] ;%γ��
    lon1=[112:0.1:121.5];
    lat1=[30.5:0.1:36.5];
    [x,y]=meshgrid(lon,lat);
    [X,Y]=meshgrid(lon1,lat1);
    Pre_NEW(:,:) = interp2(x,y,Pre_new(:,:),X,Y,'bilinear');
    Pre_NEW00(:,:)= Pre_NEW(6:57,3:94);
    [mask,Q]=geotiffread('F:\huaihe\hh_matlab01.tif');%
    mask=double(mask);
    mask(mask==0)=nan;
    Pre_end(:,:,i)= Pre_NEW00.*mask;
    %h=imagesc(Pre_end(:,:,i))
    clear Pre_NEW Pre_new
end

%%
%%%%�������뽫���������۲����ݵ���Pre_end��ģʽ���ݵĵ���ͬ��

date=ncinfo('F:\cmip6\CMCC-CM2-SR5\pr_day_CMCC-CM2-SR5_ssp585_r1i1p1f1_gn_20150101-20391231.nc');
lon=ncread('F:\cmip6\CMCC-CM2-SR5\pr_day_CMCC-CM2-SR5_ssp585_r1i1p1f1_gn_20150101-20391231.nc','lon');%����
lat=ncread('F:\cmip6\CMCC-CM2-SR5\pr_day_CMCC-CM2-SR5_ssp585_r1i1p1f1_gn_20150101-20391231.nc','lat');%γ��
time=ncread('F:\cmip6\CMCC-CM2-SR5\pr_day_CMCC-CM2-SR5_ssp585_r1i1p1f1_gn_20150101-20391231.nc','time');
pr=ncread('F:\cmip6\CMCC-CM2-SR5\pr_day_CMCC-CM2-SR5_ssp585_r1i1p1f1_gn_20150101-20391231.nc','pr');
%epr=nan(120,180,18250); %��γ�ȣ����ȣ����ݼ���������
for i=1:7300 %������ĿǰҪ��2015-2034�������
    epr(:,:,i)=flipud ( rot90 ( pr(:,:,i) ) ); %��ʱ����ת90��
end
study_area = [112.141102328 121.301102328 30.9867570496 36.1367570496];%Ϊͳһ����γ�� �޸ĵ��о�����Χ
%%study_area = [112.141102328 121.341102328 30.9867570496 36.1867570496];%ʵ�ʵ��о�����Χ
 [mask,Q]=geotiffread('F:\huaihe\hh_matlab01.tif');%
  mask=double(mask);
  mask(mask==0)=nan;
scale_want=0.1;
[lon_want,lat_want]=meshgrid( study_area(1) : scale_want : study_area(2), study_area(3): scale_want : study_area(4));
[lon_original,lat_original]=meshgrid(lon',lat');
ch_pr=nan(52,92,7300);%��Ӧ��С����γ��Ҫ�ٷ�һ��
for i=1:7300%����
    ch_pr(:,:,i)=interp2(lon_original,lat_original,epr(:,:,i),lon_want,lat_want);%˫���Բ�ֵ
end
clear pr lon_want lat_want ;%����������ʾ����
modelpre= ch_pr*86400;
%%�˲������2015-2034ģʽ��ˮ����һ��Ҫ��������µĽ�ˮ
JUN=[];
for n=1:20 %20�� 2015-2034
Julpre_n=modelpre(:,:,182+365*(n-1):212+365*(n-1));%���7��ÿ��Ľ�ˮ
Julpremon_n=sum(Julpre_n,3);
JUN=[JUN Julpremon_n];%������20���7�½�ˮ
end
Julpre_end=reshape(JUN,52,92,20);%����20���7�½�ˮ�����ʾ

%C= Julpre_end(:,:,1).*mask;
%h=imagesc(C);  ������� Ҳû������ ģʽ���ݿɳ�ͼ
clear Julpremon_n Julpre_n epr ch_pr;
%%���� 2015-2034��20���7��ģʽ��ˮ���������
%%
load('F:\cmip6С����0915\CMCC\7 �ض��·�ȫ����ƫ��У��\hismodp.mat'); %�����ch_pr��1850-2014������½�ˮ
hismodp=ch_pr;
clear ch_pr;
Jhismodp=[];
for m=1:165
    Jhismodp_m=hismodp(:,:,7+12*(m-1));
    Jhismodp=[Jhismodp Jhismodp_m];
end
Jhismodp=Jhismodp([1:52],[13341: 15180]) ;
Jhismodp=Jhismodp*2678400; %�ɴ˵õ�ģʽ������1995-2014���ֵ 52*1840=52*92*20��

%%
%������������˹۲⡢ģʽ������7�µ�ֵ������
%1.��������Pre_end 52*92*41 (1980-2020)
%2.δ��ģʽ����Julpre_end 52*92*20(2015-2034)
%3.��ʷģʽ����Jhismodp 52*1840 (1995-2014)
%������Ҫ�������������ݰ�ʱ�����к�դ��������У�
% OBS=reshape(Pre_end,[52,3772]); %�����������ŵ���״ͼ 1980-2020 41*92
% MODfur=reshape(Julpre_end,[52,1840]); %�����������ŵ���״ͼ 2015-2034 20*92
% MODhis=Jhismodp; %�����������ŵ���״ͼ 1995-2014 20*92
OBS=Pre_end;
MODfur=Julpre_end;
Jhismodp=reshape(Jhismodp,52,92,20); %ת����ά
%D= Jhismodp(:,:,1).*mask;
%h=imagesc(D);  ������� ��ʷģʽ����Ҳ���Գ�����������ͼ
MODhis=Jhismodp;
MF=[];
MH=[];
OHandF=[];
for oo=1:size(MODfur,1)%�У�92
    for jj=1:size(MODfur,2)%�У�52
        for pp=1:20%ҳ��20
            J=MODfur(oo,jj,pp);
            MF=[MF J]; 
        end
    end
end
%MFΪ2015-2034���ģʽ���� ������ٰ����������ȡԪ�� ��1-20Ϊ�������Ͻǵ�һ��Ԫ������20��
for oo=1:size(OBS,1)%�У�92
    for jj=1:size(OBS,2)%�У�52
        for pp=1:41%ҳ��41
            O=OBS(oo,jj,pp);
            OHandF=[OHandF O];
        end
    end
end
%OHandFΪ1980-2020��Ĺ۲����� ������ٰ����������ȡԪ�� ��1-40Ϊ�������Ͻǵ�һ��Ԫ������40��
for oo=1:size(MODhis,1)
    for jj=1:size(MODhis,2)
        for pp=1:size(MODhis,3)
            JH=MODhis(oo,jj,pp);
            MH=[MH JH];
        end
    end
end
%MHΪ1995-2014���ģʽ���� ������ٰ����������ȡԪ�� ��1-20Ϊ�������Ͻǵ�һ��Ԫ������20��
clear J O JH time study_are scale_want;

%%��������ȡ1995-2014��Ĺ۲�����
OBS_part=OBS(:,:,16:35);
Opart=[];
for oo=1:size(OBS_part,1)
    for jj=1:size(OBS_part,2)
        for pp=1:size(OBS_part,3)
            OOO=OBS_part(oo,jj,pp);
            Opart=[Opart OOO];
        end
    end
end
clear OOO;

