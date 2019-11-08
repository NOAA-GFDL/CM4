
clear all
nc64startup;

nc_96_static  = netcdf('/archive/oar.gfdl.cmip6/CM4/warsaw_201710_om4_v1.0.1/CM4_piControl_C/gfdl.ncrc4-intel16-prod-openmp/pp/ocean_monthly_1x1deg/ocean_monthly_1x1deg.static.nc', 'nowrite');

%nc_96_static  = netcdf('/data_cmip6/CMIP6/CMIP/NOAA-GFDL/GFDL-CM4/piControl/r1i1p1f1/Ofx/sftof/gr/v20180701/sftof_Ofx_GFDL-CM4_piControl_r1i1p1f1_gr.nc', 'nowrite');

nc_PI_ts = netcdf('/net2/h1g/CM4_paper/cm4_data/PI_C/TOS/atmos.000101-065012.tos.nc', 'nowrite')
nc_AMIP_ts = netcdf('/archive/h1g/Had_sst/20190708/fms/hadisst_sst.data_187001_201904.nc', 'nowrite')
                                                       

rad = pi/180;
v_lat   = ncvar('lat'    ,   nc_96_static );
v_lon   = ncvar('lon'    ,   nc_96_static );
v_lat_bnds  = ncvar('lat_bnds'    ,  nc_96_static );
v_lon_bnds  = ncvar('lon_bnds'    ,  nc_96_static );

v_ocean  = ncvar('sftof' ,   nc_96_static);
lat   = v_lat(:);
lon   = v_lon(:);
lat_bnds   = v_lat_bnds(:,:);
lon_bnds   = v_lon_bnds(:,:);

ocean = v_ocean;
ocean = ocean* 0.01;
ocean  = floor(ocean);
for j = 1:size(lat,1),    
  for i = 1:size(lon,1),
   if ocean(j,i) > 1.e10
      ocean(j,i)= 0.0;
   end
  end
end 

nlat = size(lat,1);
nlon  = size(lon,1);
for j = 1:nlat
  wt(j) = 0.5*(sin(rad*lat_bnds(j,2)) - sin(rad*lat_bnds(j,1)));
end


nfilt = 1;
seas = 'ann';
lim_globe = [-90, 90, -180, 180];
lim_nh    = [ 0,  90, -180, 180];
lim_sh    = [-90,  0, -180, 180];
lim_sh_ex = [-90,-45, -180, 180];
lim_nh_ex = [ 30, 90, -180, 180];
lim_trop  = [-30, 30, -180, 180];
lim_rest  = [-45, 90, -180, 180];
lim_ant  = [-90, -75, -180, 180];


PIyears = 650;
PIyear = [1:PIyears];
PIyear5 = [3:PIyears-2];


% using HadISST from 187001-201412
AMIPyears = 145;

% year 1850 corresponds to CM4.0 piControl 251
% year 1869 corresponds to CM4.0 piControl 270
% year 1870 corresponds to CM4.0 piControl 271
% year 1880 corresponds to CM4.0 piControl 281
AMIPyear = [1:AMIPyears] + 270;
AMIPyear5 = [3:AMIPyears-2] + 270;

mask = ocean;

%%% ==========================================================
% tos unit: Celcius
ts_PI_ann_ocean = series_mask(PIyears, 0, lat, lon, wt, nc_PI_ts, 'tos', lim_globe, 'ann', nfilt, mask);
ts_PI_ann_ocean = ts_PI_ann_ocean + 273.15   % convert from C to K
ts_PI_ann_ocean5 = filt5(ts_PI_ann_ocean);

% sst unit: Kelvin
ts_AMIP_ann_ocean = series_mask(AMIPyears, 1869, lat, lon, wt, nc_AMIP_ts, 'sst', lim_globe, 'ann', nfilt, mask);
ts_AMIP_ann_ocean5 = filt5(ts_AMIP_ann_ocean);

%%% ==========================================================
PImean = mean(ts_PI_ann_ocean(281:301))
AMIPmean = mean(ts_AMIP_ann_ocean(11:31))
PIbias = PImean - AMIPmean

figure
vaxis = [-20, 680, 289.5, 292.5];
plot(PIyear5,  ts_PI_ann_ocean5, '-r','LineWidth', 3) 
hold
plot(AMIPyear5, ts_AMIP_ann_ocean5, '-k','LineWidth', 3) 
plot([151,650], [PImean,PImean],'-r','LineWidth', 1)
plot([151,650], [AMIPmean,AMIPmean],'-k','LineWidth', 1)
%plot(PIyear(1:3),  ts_PI_ann_ocean(1:3), '-r','LineWidth', 3) 
plot([151,151],[200,300], '-k','LineWidth', 2) 
plot([650,650],[200,300], '-k','LineWidth', 2) 
plot([251,251],[200,300], '--k','LineWidth', 0.8) 
plot([415,415],[200,300], '--k','LineWidth', 0.8) 
axis(vaxis);
xlabel('Year','FontWeight','bold','FontSize',18);
ylabel('Gobal Mean SST (K)','FontWeight','bold','FontSize',18);
set(gca, 'FontWeight','bold','FontSize',12)

% Arrow with two head at both end and text between
x = [0.769 0.869];    % adjust length and location of arrow 
y = [0.85 0.85];    
Xadj = 1.19;      % adjust location of left arrow starting point (the sum of this with 'x' should not be negative)
annotation('textarrow',x,y,'String',' CM4.0 piControl years 1-500 ','FontWeight','bold','FontSize',14,'Linewidth',2)
annotation('textarrow',-x+Xadj,y,'String','','FontSize',14,'Linewidth',2)

% Arrow with two head at both end and text between
x1 = [0.575 0.61];    % adjust length and location of arrow 
y1 = [0.78 0.78];    
Xadj1 = 1.04;      % adjust location of left arrow starting point (the sum of this with 'x' should not be negative)
annotation('textarrow',x1,y1,'String','historical', 'FontWeight', 'bold','FontSize',14,'Linewidth',2)
annotation('textarrow',-x1+Xadj1,y1,'String','','FontSize',14,'Linewidth',2)

legend({'CM4.0 piControl','HadISST'},'Position', [0.69 0.2 0.1 0.05],'FontWeight','bold','FontSize',13);
legend boxoff

print -depsc   tos_PI_HadISST.eps

