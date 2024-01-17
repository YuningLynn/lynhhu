clear;clc;
%%
ncols=128; %列数
nrows=72;%行数
%%读入气象数据%%
readFilePath='F:\1961-2020年0.5x0.5气象 网格月数据\算例4\*07.txt';  %1980-1994
readPathStr='F:\1961-2020年0.5x0.5气象 网格月数据\算例4\';
fileList=dir(readFilePath);
fileNum=length(fileList);%文件个数
for i=1:fileNum
    name=fileList(i).name;
    splitName=strsplit(name,'.txt');  %在.处截取.前面的字符串
    varStr = splitName{1};
    fileName=strcat(readPathStr,name);%这个语句 就是获得了这个文件的完整路径
    fid=fopen(fileName,'r');
    FormatString=repmat('%f',1,ncols);
    Pre(:,:,i)=cell2mat(textscan(fid,FormatString,nrows,'HeaderLines',6));%读取中国范围内05°格网降水数据  不读前六行
    Pre_new=Pre(35:47,81:100,i);%淮河区对应的行列 二维数组
    Pre_new(Pre_new<0)=nan;
    lon=[112:0.5:121.5]; %经度
    lat=[30.5:0.5:36.5] ;%纬度
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
%%%%上述代码将多年的气象观测数据导入Pre_end，模式数据的导入同理

date=ncinfo('F:\cmip6\CMCC-CM2-SR5\pr_day_CMCC-CM2-SR5_ssp585_r1i1p1f1_gn_20150101-20391231.nc');
lon=ncread('F:\cmip6\CMCC-CM2-SR5\pr_day_CMCC-CM2-SR5_ssp585_r1i1p1f1_gn_20150101-20391231.nc','lon');%经度
lat=ncread('F:\cmip6\CMCC-CM2-SR5\pr_day_CMCC-CM2-SR5_ssp585_r1i1p1f1_gn_20150101-20391231.nc','lat');%纬度
time=ncread('F:\cmip6\CMCC-CM2-SR5\pr_day_CMCC-CM2-SR5_ssp585_r1i1p1f1_gn_20150101-20391231.nc','time');
pr=ncread('F:\cmip6\CMCC-CM2-SR5\pr_day_CMCC-CM2-SR5_ssp585_r1i1p1f1_gn_20150101-20391231.nc','pr');
%epr=nan(120,180,18250); %（纬度，经度，数据集的天数）
for i=1:7300 %天数，目前要做2015-2034年的天数
    epr(:,:,i)=flipud ( rot90 ( pr(:,:,i) ) ); %逆时针旋转90°
end
study_area = [112.141102328 121.301102328 30.9867570496 36.1367570496];%为统一数组纬度 修改的研究区范围
%%study_area = [112.141102328 121.341102328 30.9867570496 36.1867570496];%实际的研究区范围
 [mask,Q]=geotiffread('F:\huaihe\hh_matlab01.tif');%
  mask=double(mask);
  mask(mask==0)=nan;
scale_want=0.1;
[lon_want,lat_want]=meshgrid( study_area(1) : scale_want : study_area(2), study_area(3): scale_want : study_area(4));
[lon_original,lat_original]=meshgrid(lon',lat');
ch_pr=nan(52,92,7300);%对应大小，经纬度要再反一下
for i=1:7300%天数
    ch_pr(:,:,i)=interp2(lon_original,lat_original,epr(:,:,i),lon_want,lat_want);%双线性插值
end
clear pr lon_want lat_want ;%工作区不显示变量
modelpre= ch_pr*86400;
%%此步骤完成2015-2034模式降水，下一步要算各年七月的降水
JUN=[];
for n=1:20 %20年 2015-2034
Julpre_n=modelpre(:,:,182+365*(n-1):212+365*(n-1));%求出7月每天的降水
Julpremon_n=sum(Julpre_n,3);
JUN=[JUN Julpremon_n];%是连续20年的7月降水
end
Julpre_end=reshape(JUN,52,92,20);%连续20年的7月降水逐个显示

%C= Julpre_end(:,:,1).*mask;
%h=imagesc(C);  检查至此 也没有问题 模式数据可出图
clear Julpremon_n Julpre_n epr ch_pr;
%%到此 2015-2034这20年的7月模式降水已输入完成
%%
load('F:\cmip6小论文0915\CMCC\7 特定月份全流域偏差校正\hismodp.mat'); %载入的ch_pr是1850-2014年的逐月降水
hismodp=ch_pr;
clear ch_pr;
Jhismodp=[];
for m=1:165
    Jhismodp_m=hismodp(:,:,7+12*(m-1));
    Jhismodp=[Jhismodp Jhismodp_m];
end
Jhismodp=Jhismodp([1:52],[13341: 15180]) ;
Jhismodp=Jhismodp*2678400; %由此得到模式数据在1995-2014年的值 52*1840=52*92*20年

%%
%上述内容完成了观测、模式数据在7月的值，包括
%1.气象数据Pre_end 52*92*41 (1980-2020)
%2.未来模式数据Julpre_end 52*92*20(2015-2034)
%3.历史模式数据Jhismodp 52*1840 (1995-2014)
%接下来要把上述三个数据按时间序列和栅格次序排列：
% OBS=reshape(Pre_end,[52,3772]); %淮河流域连着的形状图 1980-2020 41*92
% MODfur=reshape(Julpre_end,[52,1840]); %淮河流域连着的形状图 2015-2034 20*92
% MODhis=Jhismodp; %淮河流域连着的形状图 1995-2014 20*92
OBS=Pre_end;
MODfur=Julpre_end;
Jhismodp=reshape(Jhismodp,52,92,20); %转成三维
%D= Jhismodp(:,:,1).*mask;
%h=imagesc(D);  检查至此 历史模式数据也可以出完整的流域图
MODhis=Jhismodp;
MF=[];
MH=[];
OHandF=[];
for oo=1:size(MODfur,1)%行：92
    for jj=1:size(MODfur,2)%列：52
        for pp=1:20%页：20
            J=MODfur(oo,jj,pp);
            MF=[MF J]; 
        end
    end
end
%MF为2015-2034年的模式数据 按年份再按矩阵横着提取元素 如1-20为矩阵左上角第一个元素连续20年
for oo=1:size(OBS,1)%列：92
    for jj=1:size(OBS,2)%行：52
        for pp=1:41%页：41
            O=OBS(oo,jj,pp);
            OHandF=[OHandF O];
        end
    end
end
%OHandF为1980-2020年的观测数据 按年份再按矩阵横着提取元素 如1-40为矩阵左上角第一个元素连续40年
for oo=1:size(MODhis,1)
    for jj=1:size(MODhis,2)
        for pp=1:size(MODhis,3)
            JH=MODhis(oo,jj,pp);
            MH=[MH JH];
        end
    end
end
%MH为1995-2014年的模式数据 按年份再按矩阵横着提取元素 如1-20为矩阵左上角第一个元素连续20年
clear J O JH time study_are scale_want;

%%以下是提取1995-2014年的观测数据
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

