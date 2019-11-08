nc64startup

f=netcdf('sss_subp_921.nc')
sss_921=f{'SSS'}(:);
f=netcdf('sss_subp_969.nc')
sss_969=f{'SSS'}(:);
f=netcdf('sss_subp_970.nc')
sss_970=f{'SSS'}(:);
f=netcdf('sss_subp_ave.nc')
sss_ave=f{'SSS_AVE'}(:);

f=netcdf('SSS_60W0E_50N65N_2015.nc')
sss_obs_5yr=f{'SSS'}(1:56);

for i=1:56
tyear_obs(i)=1956+i;
end

for i=1:165
tyear(i)=1850+i-1;
end

for i=1:161
sss_921_5yr(i)=mean(sss_921(i:i+4));
sss_969_5yr(i)=mean(sss_969(i:i+4));
sss_970_5yr(i)=mean(sss_970(i:i+4));
sss_ave_5yr(i)=mean(sss_ave(i:i+4));
end

dsss_ave=detrend(sss_ave_5yr);
trend_sss_ave=sss_ave_5yr-dsss_ave;

figure
plot(tyear(3:163),sss_ave_5yr-mean(sss_ave_5yr(106:161)),'b','LineWidth',2)
hold
plot(tyear(3:163),trend_sss_ave-mean(trend_sss_ave(106:161)),'b-','LineWidth',2)
plot(tyear(3:163),sss_921_5yr-mean(sss_921_5yr(106:161)),'b-')
plot(tyear(3:163),sss_969_5yr-mean(sss_969_5yr(106:161)),'b-')
plot(tyear(3:163),sss_970_5yr-mean(sss_970_5yr(106:161)),'b-')
plot(tyear_obs,detrend(sss_obs_5yr,'constant'),'r','LineWidth',2)
axis([1852 2012 -0.12 0.1])
title('Subpolar North Atlantic SSS Anomalies (PSU)')
xlabel('Year')
grid



