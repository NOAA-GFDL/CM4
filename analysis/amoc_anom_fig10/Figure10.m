nc64startup

f=netcdf('amoc_26.5N_z_151_650.nc')
moc_control=f{'AMOC'}(:);

for i=1:500
tyear_con(i)=i;
end

for i=1:496
moc_control_5yr(i)=mean(moc_control(i:i+4));
end

figure
plot(tyear_con(3:498),detrend(moc_control_5yr,'constant'),'k','LineWidth',2)
axis([0 500 -3 3])
title('AMOC Anomalies at 26N (Sv)')
xlabel('Year')
grid
