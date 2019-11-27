nc64startup

f=netcdf('amoc_151_650_bottom_ave_26_5N.nc')
moc_980=f{'AMOC_AVE'}(:);
moc_980=[1.2,moc_980'];

f=netcdf('amoc_sm_151_650_bottom_ave_26_5N.nc')
moc_sm_980=f{'AMOC_SM_AVE'}(:);
moc_sm_980=[0,moc_sm_980'];

f=netcdf('moc_vertical_ave_2004_2015.nc')
moc_rapid=f{'MOC_AVE'}(:);
depth_rapid=-f{'DEPTH'}(:);

depth=-[0 5 15 25 40 62.5 87.5 112.5 137.5 175 225 275 350 450 550 650 750 850 950 1050 1150 1250 1350 1450 1625 1875 2250 2750 3250 3750 4250 4750 5250 5750 6250 6750];

figure
plot(moc_980-moc_sm_980,depth,'r','LineWidth',2)
hold
plot(moc_rapid,depth_rapid,'k','LineWidth',3)
legend ('CM4.0','RAPID')
axis([-5 20 -6000 0])
xlabel('AMOC at 26N (Sv)')
ylabel('Depth (m)')
grid





