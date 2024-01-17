getwd() 
library(readxl)
MON=read_excel("F:/Cmip6小论文0915/CMCC/5 QQM-QUANT校正/1995-2014CMCC与观测数据.xlsx",sheet=1) #月份
FUTUREmod=read_excel("F:/Cmip6小论文0915/CMCC/5 QQM-QUANT校正/ssp126/ssp126未来模式数据.xlsx",sheet=1) #未来模式
FUTUREobs=read_excel("F:/Cmip6小论文0915/CMCC/5 QQM-QUANT校正/ssp126/ssp126未来模式数据.xlsx",sheet=3) #未来观测
pa = unlist(c(MON[1])) #1995-2014obs
pm = unlist(c(MON[2])) #1995-2014 CMCC model
pmf = unlist(c(FUTUREmod[1])) #2015-2020 CMCCmodel
paf = unlist(c(FUTUREobs[1])) #2015-2020 obs
op <- par(mfrow=c(2,3))
#qqplot(pa,pm)
#abline(a=0,b=1,col="red")#畫一條紅色的對角線y=x，a是intercept，b是slope
ss1= seq(0.00000001,1, length = 20) 
ss2= seq(0.00000001,1, length = 6) 
pa.quan= quantile(pa, ss1) #ss个分位点及其对应的值  原序列重新排序
pm.quan= quantile(pm, ss1) #ss个分位点及其对应的值
paf.quan=quantile (paf,ss2)
pmf.quan=quantile(pmf,ss2)
new=pmf.quan

par(mfrow=c(2,2),mai=c(0.5,0.5,0.2,0.2))

#for (TF in 1:2)
#{
  
  #if(TF==1) pm.tf = lm(pa.quan~pm.quan+0)  #一次函数 
  #if(TF==1) pm.tf = lm(pa.quan~pm.quan)    #正比例函数 y=ax+b  pm.tf是传递函数的系数
  #if(TF==2) pm.tf = nls(pa.quan ~ b*pm.quan^c, start=list(b=0.5,c=0.5))
  #if(TF==4) 
   # pm.tf = nls(pa.quan ~ b*(pm.quan-x0)^c, start=list(b=0.5,c=0.6,x0=0))
  
  
  #pm.quan.tf =predict(pm.tf) #订正后的排序模式值 1995-2014
   model<-nls(pa.quan ~ b*pm.quan^c, start=list(b=0.5,c=0.5))
   pmf.quan.tf =predict(model,newdata=new) #订正后的排序模式值 2015-2020
  #pm.quan.tf = predict(pm.tf, data.frame(pm.quan=pm.quan))
  
  #plot(paf.quan, pmf.quan.tf) #观测值与校正后的模式值
  #abline(a=0,b=1,col="red")
#}

inverse_ecdf <- function(x, prob) { 
  if (is.unsorted(x)) x <- sort(x)
  n <- length(x)
  approx(seq(1/(n+1), n/(n+1), length = n), x, prob)$y #cdf的逆函数
}

Fn.pmf <- ecdf(pmf)

pmf.tf = inverse_ecdf(pmf.quan.tf, Fn.pmf(pmf))#導出補點後的y 根据历史模式的概率反推对应频率下的订正模式值
plot(pmf, pmf.tf)  #(x为订正前 y为订正后)
abline(0,1)


par(mfrow=c(2,1),mai=c(0.8,0.8,0.6,0.7))
qqplot(paf, pmf,main='校正前观测降雨与模式降雨对比图',xlab="obs_pre(mm)",ylab="mod_pre(mm)"); abline(a=0,b=1,col="red")     #订正前
qqplot(paf, pmf.tf,main='校正后观测降雨与模式降雨对比图',xlab="obs_pre(mm)",ylab="mod_pre(mm)"); abline(a=0,b=1,col="red")  #订正后

df<-data.frame(pmf.tf)
write.csv(df,"F:/Cmip6小论文0915/CMCC/4 QQM-TF校正/ssp126/QUANT15-20 JAN订正",row.names= FALSE)