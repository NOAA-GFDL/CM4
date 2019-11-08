nc64startup

f=netcdf('moc_atl_921_26_5N_max.nc')
moc_921=f{'MOC'}(:);
f=netcdf('moc_atl_969_26_5N_max.nc')
moc_969=f{'MOC'}(:);
f=netcdf('moc_atl_970_26_5N_max.nc')
moc_970=f{'MOC'}(:);
f=netcdf('moc_atl_ave_26_5N.nc')
moc_ave=f{'MOC_AVE'}(:);

for i=1:165
tyear(i)=1850+i-1;
end

for i=1:161
moc_921_5yr(i)=mean(moc_921(i:i+4));
moc_969_5yr(i)=mean(moc_969(i:i+4));
moc_970_5yr(i)=mean(moc_970(i:i+4));
moc_ave_5yr(i)=mean(moc_ave(i:i+4));
end

dmoc_ave=detrend(moc_ave_5yr);
trend_moc_ave=moc_ave_5yr-dmoc_ave;

figure
plot(tyear(3:163),detrend(moc_ave_5yr,'constant'),'b','LineWidth',2)
hold
plot(tyear(3:163),detrend(trend_moc_ave,'constant'),'b-','LineWidth',2)
plot(tyear(3:163),detrend(moc_921_5yr,'constant'),'b-')
plot(tyear(3:163),detrend(moc_969_5yr,'constant'),'b-')
plot(tyear(3:163),detrend(moc_970_5yr,'constant'),'b-')
axis([1852 2012 -3 3])
title('AMOC Anomalies at 26N (Sv)')
xlabel('Year')
grid

