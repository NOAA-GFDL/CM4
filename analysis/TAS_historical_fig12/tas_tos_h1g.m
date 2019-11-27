
clear all;

nc64startup;

%nc_96_static  = netcdf('/net2/h1g/CM4_paper/cm4_data/c96/atmos.static.nc', 'nowrite');

nc_96_static  = netcdf('/archive/oar.gfdl.cmip6/CM4/warsaw_201710_om4_v1.0.1/CM4_piControl_C/gfdl.ncrc4-intel16-prod-openmp/pp/ocean_monthly_1x1deg/ocean_monthly_1x1deg.static.nc', 'nowrite');

nc_PI_tas    = netcdf('/net2/h1g/CM4_paper/cm4_data/PI_C/TAS/tas_360x180.nc', 'nowrite');
nc_PI_tos    = netcdf('/net2/h1g/CM4_paper/cm4_data/PI_C/TOS/atmos.000101-065012.tos.nc', 'nowrite');
nc_PI_seaice = netcdf('/net2/h1g/CM4_paper/cm4_data/PI_C/Sea_Ice_Frac/ice_1x1deg.000101-065012.siconc.nc', 'nowrite');

nc_H1_tas    = netcdf('/net2/h1g/CM4_paper/cm4_data/H1/TAS/tas_360x180.nc', 'nowrite');
nc_H1_tos    = netcdf('/net2/h1g/CM4_paper/cm4_data/H1/TOS/ocean_monthly_1x1deg.185001-201412.tos.nc', 'nowrite');
nc_H1_seaice = netcdf('/net2/h1g/CM4_paper/cm4_data/H1/Sea_Ice_Frac/ice_1x1deg.185001-201412.siconc.nc', 'nowrite');

nc_H2_tas    = netcdf('/net2/h1g/CM4_paper/cm4_data/H2/TAS/tas_360x180.nc', 'nowrite');
nc_H2_tos    = netcdf('/net2/h1g/CM4_paper/cm4_data/H2/TOS/ocean_monthly_1x1deg.185001-201412.tos.nc', 'nowrite');
nc_H2_seaice = netcdf('/net2/h1g/CM4_paper/cm4_data/H2/Sea_Ice_Frac/ice_1x1deg.185001-201412.siconc.nc', 'nowrite');

nc_H3_tas    = netcdf('/net2/h1g/CM4_paper/cm4_data/H3/TAS/tas_360x180.nc', 'nowrite');
nc_H3_tos    = netcdf('/net2/h1g/CM4_paper/cm4_data/H3/TOS/ocean_monthly_1x1deg.185001-201412.tos.nc', 'nowrite');
nc_H3_seaice = netcdf('/net2/h1g/CM4_paper/cm4_data/H3/Sea_Ice_Frac/ice_1x1deg.185001-201412.siconc.nc', 'nowrite');

nc_gistemp = netcdf('/net2/h1g/CM4_paper/obsdata/gistemp1200_ERSSTv4.nc', 'nowrite');

v_lat_giss   = ncvar('lat'    ,   nc_gistemp);
v_lon_giss   = ncvar('lon'    ,   nc_gistemp);
load /net2/h1g/CM4_paper/obsdata/giss_nh
load /net2/h1g/CM4_paper/obsdata/giss_sh

giss_year = giss_nh(:,1);
giss_nh_ann = giss_nh(:,14);
giss_sh_ann = giss_sh(:,14);

giss_year5 = [1883:2013];

giss_nh_ann5 = filt5(giss_nh(:,14));
giss_sh_ann5 = filt5(giss_sh(:,14));

giss_nh_ann = giss_nh_ann/100;
giss_sh_ann = giss_sh_ann/100;
giss_globe_ann = (giss_nh_ann + giss_sh_ann)/2;

giss_nh_ann5 = giss_nh_ann5/100;
giss_sh_ann5 = giss_sh_ann5/100;
giss_globe_ann5 = (giss_nh_ann5 + giss_sh_ann5)/2;

gy1 = 1;
gy2 = 21;
n_giss_nh_ann = giss_nh_ann - mean(giss_nh_ann(gy1:gy2));
n_giss_sh_ann = giss_sh_ann - mean(giss_sh_ann(gy1:gy2));
n_giss_globe_ann = giss_globe_ann - mean(giss_globe_ann(gy1:gy2));

n_giss_nh_ann5 = giss_nh_ann5 - mean(giss_nh_ann(gy1:gy2));
n_giss_sh_ann5 = giss_sh_ann5 - mean(giss_sh_ann(gy1:gy2));
n_giss_globe_ann5 = giss_globe_ann5 - mean(giss_globe_ann(gy1:gy2));

rad = pi/180;
v_lat   = ncvar('lat'    ,   nc_96_static);
v_lon   = ncvar('lon'    ,   nc_96_static);
v_lat_bnds  = ncvar('lat_bnds'    ,   nc_96_static);
v_lon_bnds  = ncvar('lon_bnds'    ,   nc_96_static);

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
lim_nh_ex = [ 45, 90, -180, 180];
lim_trop  = [-45, 45, -180, 180];
lim_rest  = [-45, 90, -180, 180];
lim_ant  = [-90, -75, -180, 180];

PIyears = 650;
PIyear = [1:PIyears] + 1850 - 251;
PIyear5 = [3:PIyears-2] + 1850 - 251;
PIyear11 = [6:PIyears-5] + 1850 - 251;

Hyears = 165;
yearH = [1850:2014];
yearH5 = yearH(3:163);

%%% ==========================================================

%tas_H1_ann_globe = series_globe(165, 1850, lat, lon, wt, nc_H1_tas, 'tas', lim_globe, 'ann', nfilt);
%tas_H2_ann_globe = series_globe(165, 1850, lat, lon, wt, nc_H2_tas, 'tas', lim_globe, 'ann', nfilt);
%tas_H3_ann_globe = series_globe(165, 1850, lat, lon, wt, nc_H3_tas, 'tas', lim_globe, 'ann', nfilt);
%tas_PI_ann_globe = series_globe(PIyears, 0, lat, lon, wt, nc_PI_tas, 'tas', lim_globe, 'ann', nfilt);

tas_H1_ann_globe = series_globe_tas_tos(165, 1850, lat, lon, wt, nc_H1_tas, 'tas', nc_H1_tos, 'tos', nc_H1_seaice, 'siconc', ocean,lim_globe, 'ann', nfilt);
tas_H2_ann_globe = series_globe_tas_tos(165, 1850, lat, lon, wt, nc_H2_tas, 'tas', nc_H2_tos, 'tos', nc_H2_seaice, 'siconc', ocean,lim_globe, 'ann', nfilt);
tas_H3_ann_globe = series_globe_tas_tos(165, 1850, lat, lon, wt, nc_H3_tas, 'tas', nc_H3_tos, 'tos', nc_H3_seaice, 'siconc', ocean,lim_globe, 'ann', nfilt);
tas_PI_ann_globe = series_globe_tas_tos(PIyears, 0, lat, lon, wt,nc_PI_tas, 'tas', nc_PI_tos, 'tos', nc_PI_seaice, 'siconc', ocean,lim_globe, 'ann', nfilt);

tas_H1_ann_globe5 = filt5(tas_H1_ann_globe);
tas_H2_ann_globe5 = filt5(tas_H2_ann_globe);
tas_H3_ann_globe5 = filt5(tas_H3_ann_globe);
tas_PI_ann_globe5 = filt5(tas_PI_ann_globe);

tas_PI_ann_globe11 = filt11(tas_PI_ann_globe);

ya = 31;
yb = 51;

tas_H_ann_globe = [tas_H1_ann_globe, tas_H2_ann_globe ,tas_H3_ann_globe]; 
tas_HM_ann_globe = mean(tas_H_ann_globe,2);
tas_HMIN_ann_globe = min(tas_H_ann_globe, [],2);
tas_HMAX_ann_globe = max(tas_H_ann_globe, [],2);

tas_H_ann_globe5 = [tas_H1_ann_globe5, tas_H2_ann_globe5 ,tas_H3_ann_globe5]; 
tas_HM_ann_globe5 = mean(tas_H_ann_globe5,2);
tas_HMIN_ann_globe5 = min(tas_H_ann_globe5, [],2);
tas_HMAX_ann_globe5 = max(tas_H_ann_globe5, [],2);

n_tas_HM_ann_globe =  tas_HM_ann_globe - mean(tas_HM_ann_globe(ya:yb));
n_tas_H1_ann_globe = tas_H1_ann_globe - mean(tas_HM_ann_globe(ya:yb));
n_tas_H2_ann_globe = tas_H2_ann_globe - mean(tas_HM_ann_globe(ya:yb));
n_tas_H3_ann_globe = tas_H3_ann_globe - mean(tas_HM_ann_globe(ya:yb));
n_tas_HMIN_ann_globe = tas_HMIN_ann_globe - mean(tas_HM_ann_globe(ya:yb));
n_tas_HMAX_ann_globe = tas_HMAX_ann_globe - mean(tas_HM_ann_globe(ya:yb));

n_tas_HM_ann_globe5 =  tas_HM_ann_globe5 - mean(tas_HM_ann_globe(ya:yb));
n_tas_H1_ann_globe5 = tas_H1_ann_globe5 - mean(tas_HM_ann_globe(ya:yb));
n_tas_H2_ann_globe5 = tas_H2_ann_globe5 - mean(tas_HM_ann_globe(ya:yb));
n_tas_H3_ann_globe5 = tas_H3_ann_globe5 - mean(tas_HM_ann_globe(ya:yb));
n_tas_HMIN_ann_globe5 = tas_HMIN_ann_globe5 - mean(tas_HM_ann_globe(ya:yb));
n_tas_HMAX_ann_globe5 = tas_HMAX_ann_globe5 - mean(tas_HM_ann_globe(ya:yb));

nnn = size(yearH,2);
for k = 1:nnn
  kk = nnn+1-k;
  yy_globe(k) = n_tas_HMAX_ann_globe(k);
  yy_globe(k+nnn) = n_tas_HMIN_ann_globe(kk);
  tt_globe(k) = yearH(k);
  tt_globe(k+nnn) = yearH(kk);
end
yy_globe(nnn+nnn+1) = n_tas_HMAX_ann_globe(1);
tt_globe(nnn+nnn+1) = yearH(1);

nnn5 = size(yearH5,2);
for k = 1:nnn5
  kk = nnn5+1-k;
  yy_globe5(k) = n_tas_HMAX_ann_globe5(k);
  yy_globe5(k+nnn5) = n_tas_HMIN_ann_globe5(kk);
  tt_globe5(k) = yearH5(k);
  tt_globe5(k+nnn5) = yearH5(kk);
end
yy_globe5(nnn5+nnn5+1) = n_tas_HMAX_ann_globe5(1);
tt_globe5(nnn5+nnn5+1) = yearH5(1);

% average piControl year 281-301, corresponding to first CM4 historical year 1880-1900
%    note piControl year 251 --> historical 1850 
ya_pi = 281;
yb_pi = 301;
n_tas_PI_ann_globe = tas_PI_ann_globe - mean(tas_PI_ann_globe(ya_pi:yb_pi));
n_tas_PI_ann_globe5 = tas_PI_ann_globe5 - mean(tas_PI_ann_globe(ya_pi:yb_pi));
n_tas_PI_ann_globe11 = tas_PI_ann_globe11 - mean(tas_PI_ann_globe(ya_pi:yb_pi));

%%% ==========================================================

%tas_H1_ann_nh = series_globe(165, 1850, lat, lon, wt, nc_H1_tas, 'tas', lim_nh, 'ann', nfilt);
%tas_H2_ann_nh = series_globe(165, 1850, lat, lon, wt, nc_H2_tas, 'tas', lim_nh, 'ann', nfilt);
%tas_H3_ann_nh = series_globe(165, 1850, lat, lon, wt, nc_H3_tas, 'tas', lim_nh, 'ann', nfilt);
%tas_PI_ann_nh = series_globe(PIyears, 0, lat, lon, wt, nc_PI_tas, 'tas', lim_nh, 'ann', nfilt);

tas_H1_ann_nh = series_globe_tas_tos(165, 1850, lat, lon, wt, nc_H1_tas, 'tas', nc_H1_tos, 'tos', nc_H1_seaice, 'siconc', ocean, lim_nh, 'ann', nfilt);
tas_H2_ann_nh = series_globe_tas_tos(165, 1850, lat, lon, wt, nc_H2_tas, 'tas', nc_H2_tos, 'tos', nc_H2_seaice, 'siconc', ocean, lim_nh, 'ann', nfilt);
tas_H3_ann_nh = series_globe_tas_tos(165, 1850, lat, lon, wt, nc_H3_tas, 'tas', nc_H3_tos, 'tos', nc_H3_seaice, 'siconc', ocean, lim_nh, 'ann', nfilt);
tas_PI_ann_nh = series_globe_tas_tos(PIyears, 0,lat, lon, wt, nc_PI_tas, 'tas', nc_PI_tos, 'tos', nc_PI_seaice, 'siconc', ocean, lim_nh, 'ann', nfilt);
 

tas_H1_ann_nh5 = filt5(tas_H1_ann_nh);
tas_H2_ann_nh5 = filt5(tas_H2_ann_nh);
tas_H3_ann_nh5 = filt5(tas_H3_ann_nh);
tas_PI_ann_nh5 = filt5(tas_PI_ann_nh);
tas_PI_ann_nh11 = filt11(tas_PI_ann_nh);

ya = 31;
yb = 51;

tas_H_ann_nh = [tas_H1_ann_nh, tas_H2_ann_nh ,tas_H3_ann_nh]; 
tas_HM_ann_nh = mean(tas_H_ann_nh,2);
tas_HMIN_ann_nh = min(tas_H_ann_nh, [],2);
tas_HMAX_ann_nh = max(tas_H_ann_nh, [],2);

tas_H_ann_nh5 = [tas_H1_ann_nh5, tas_H2_ann_nh5 ,tas_H3_ann_nh5]; 
tas_HM_ann_nh5 = mean(tas_H_ann_nh5,2);
tas_HMIN_ann_nh5 = min(tas_H_ann_nh5, [],2);
tas_HMAX_ann_nh5 = max(tas_H_ann_nh5, [],2);

n_tas_HM_ann_nh =  tas_HM_ann_nh - mean(tas_HM_ann_nh(ya:yb));
n_tas_H1_ann_nh = tas_H1_ann_nh - mean(tas_HM_ann_nh(ya:yb));
n_tas_H2_ann_nh = tas_H2_ann_nh - mean(tas_HM_ann_nh(ya:yb));
n_tas_H3_ann_nh = tas_H3_ann_nh - mean(tas_HM_ann_nh(ya:yb));
n_tas_HMIN_ann_nh = tas_HMIN_ann_nh - mean(tas_HM_ann_nh(ya:yb));
n_tas_HMAX_ann_nh = tas_HMAX_ann_nh - mean(tas_HM_ann_nh(ya:yb));


n_tas_HM_ann_nh5 =  tas_HM_ann_nh5 - mean(tas_HM_ann_nh(ya:yb));
n_tas_H1_ann_nh5 = tas_H1_ann_nh5 - mean(tas_HM_ann_nh(ya:yb));
n_tas_H2_ann_nh5 = tas_H2_ann_nh5 - mean(tas_HM_ann_nh(ya:yb));
n_tas_H3_ann_nh5 = tas_H3_ann_nh5 - mean(tas_HM_ann_nh(ya:yb));
n_tas_HMIN_ann_nh5 = tas_HMIN_ann_nh5 - mean(tas_HM_ann_nh(ya:yb));
n_tas_HMAX_ann_nh5 = tas_HMAX_ann_nh5 - mean(tas_HM_ann_nh(ya:yb));

nnn = size(yearH,2);
for k = 1:nnn
  kk = nnn+1-k;
  yy_nh(k) = n_tas_HMAX_ann_nh(k);
  yy_nh(k+nnn) = n_tas_HMIN_ann_nh(kk);
  tt_nh(k) = yearH(k);
  tt_nh(k+nnn) = yearH(kk);
end
yy_nh(nnn+nnn+1) = n_tas_HMAX_ann_nh(1);
tt_nh(nnn+nnn+1) = yearH(1);

nnn5 = size(yearH5,2);
for k = 1:nnn5
  kk = nnn5+1-k;
  yy_nh5(k) = n_tas_HMAX_ann_nh5(k);
  yy_nh5(k+nnn5) = n_tas_HMIN_ann_nh5(kk);
  tt_nh5(k) = yearH5(k);
  tt_nh5(k+nnn5) = yearH5(kk);
end
yy_nh5(nnn5+nnn5+1) = n_tas_HMAX_ann_nh5(1);
tt_nh5(nnn5+nnn5+1) = yearH5(1);


ya_pi = 281;
yb_pi = 301;
n_tas_PI_ann_nh = tas_PI_ann_nh - mean(tas_PI_ann_nh(ya_pi:yb_pi));
n_tas_PI_ann_nh5 = tas_PI_ann_nh5 - mean(tas_PI_ann_nh(ya_pi:yb_pi));
n_tas_PI_ann_nh11 = tas_PI_ann_nh11 - mean(tas_PI_ann_nh(ya_pi:yb_pi));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%tas_H1_ann_sh = series_globe(165, 1850, lat, lon, wt, nc_H1_tas, 'tas', lim_sh, 'ann', nfilt);
%tas_H2_ann_sh = series_globe(165, 1850, lat, lon, wt, nc_H2_tas, 'tas', lim_sh, 'ann', nfilt);
%tas_H3_ann_sh = series_globe(165, 1850, lat, lon, wt, nc_H3_tas, 'tas', lim_sh, 'ann', nfilt);
%tas_PI_ann_sh = series_globe(PIyears, 0, lat, lon, wt, nc_PI_tas, 'tas', lim_sh, 'ann', nfilt);

tas_H1_ann_sh = series_globe_tas_tos(165, 1850, lat, lon, wt, nc_H1_tas, 'tas', nc_H1_tos, 'tos', nc_H1_seaice, 'siconc', ocean, lim_sh, 'ann', nfilt);
tas_H2_ann_sh = series_globe_tas_tos(165, 1850, lat, lon, wt, nc_H2_tas, 'tas', nc_H2_tos, 'tos', nc_H2_seaice, 'siconc', ocean, lim_sh, 'ann', nfilt);
tas_H3_ann_sh = series_globe_tas_tos(165, 1850, lat, lon, wt, nc_H3_tas, 'tas', nc_H3_tos, 'tos', nc_H3_seaice, 'siconc', ocean, lim_sh, 'ann', nfilt);
tas_PI_ann_sh = series_globe_tas_tos(PIyears, 0,lat, lon, wt, nc_PI_tas, 'tas', nc_PI_tos, 'tos', nc_PI_seaice, 'siconc', ocean, lim_sh, 'ann', nfilt);

tas_H1_ann_sh5 = filt5(tas_H1_ann_sh);
tas_H2_ann_sh5 = filt5(tas_H2_ann_sh);
tas_H3_ann_sh5 = filt5(tas_H3_ann_sh);
tas_PI_ann_sh5 = filt5(tas_PI_ann_sh);
tas_PI_ann_sh11 = filt11(tas_PI_ann_sh);

ya = 31;
yb = 51;

tas_H_ann_sh = [tas_H1_ann_sh, tas_H2_ann_sh ,tas_H3_ann_sh]; 
tas_HM_ann_sh = mean(tas_H_ann_sh,2);
tas_HMIN_ann_sh = min(tas_H_ann_sh, [],2);
tas_HMAX_ann_sh = max(tas_H_ann_sh, [],2);

tas_H_ann_sh5 = [tas_H1_ann_sh5, tas_H2_ann_sh5 ,tas_H3_ann_sh5]; 
tas_HM_ann_sh5 = mean(tas_H_ann_sh5,2);
tas_HMIN_ann_sh5 = min(tas_H_ann_sh5, [],2);
tas_HMAX_ann_sh5 = max(tas_H_ann_sh5, [],2);

n_tas_HM_ann_sh =  tas_HM_ann_sh - mean(tas_HM_ann_sh(ya:yb));
n_tas_H1_ann_sh = tas_H1_ann_sh - mean(tas_HM_ann_sh(ya:yb));
n_tas_H2_ann_sh = tas_H2_ann_sh - mean(tas_HM_ann_sh(ya:yb));
n_tas_H3_ann_sh = tas_H3_ann_sh - mean(tas_HM_ann_sh(ya:yb));
n_tas_HMIN_ann_sh = tas_HMIN_ann_sh - mean(tas_HM_ann_sh(ya:yb));
n_tas_HMAX_ann_sh = tas_HMAX_ann_sh - mean(tas_HM_ann_sh(ya:yb));


n_tas_HM_ann_sh5 =  tas_HM_ann_sh5 - mean(tas_HM_ann_sh(ya:yb));
n_tas_H1_ann_sh5 = tas_H1_ann_sh5 - mean(tas_HM_ann_sh(ya:yb));
n_tas_H2_ann_sh5 = tas_H2_ann_sh5 - mean(tas_HM_ann_sh(ya:yb));
n_tas_H3_ann_sh5 = tas_H3_ann_sh5 - mean(tas_HM_ann_sh(ya:yb));
n_tas_HMIN_ann_sh5 = tas_HMIN_ann_sh5 - mean(tas_HM_ann_sh(ya:yb));
n_tas_HMAX_ann_sh5 = tas_HMAX_ann_sh5 - mean(tas_HM_ann_sh(ya:yb));

nnn = size(yearH,2);
for k = 1:nnn
  kk = nnn+1-k;
  yy_sh(k) = n_tas_HMAX_ann_sh(k);
  yy_sh(k+nnn) = n_tas_HMIN_ann_sh(kk);
  tt_sh(k) = yearH(k);
  tt_sh(k+nnn) = yearH(kk);
end
yy_sh(nnn+nnn+1) = n_tas_HMAX_ann_sh(1);
tt_sh(nnn+nnn+1) = yearH(1);

nnn5 = size(yearH5,2);
for k = 1:nnn5
  kk = nnn5+1-k;
  yy_sh5(k) = n_tas_HMAX_ann_sh5(k);
  yy_sh5(k+nnn5) = n_tas_HMIN_ann_sh5(kk);
  tt_sh5(k) = yearH5(k);
  tt_sh5(k+nnn5) = yearH5(kk);
end
yy_sh5(nnn5+nnn5+1) = n_tas_HMAX_ann_sh5(1);
tt_sh5(nnn5+nnn5+1) = yearH5(1);


ya_pi = 281;
yb_pi = 301;
n_tas_PI_ann_sh = tas_PI_ann_sh - mean(tas_PI_ann_sh(ya_pi:yb_pi));
n_tas_PI_ann_sh5 = tas_PI_ann_sh5 - mean(tas_PI_ann_sh(ya_pi:yb_pi));
n_tas_PI_ann_sh11 = tas_PI_ann_sh11 - mean(tas_PI_ann_sh(ya_pi:yb_pi));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ==== began of h1g added, 2019-07-02 ======
% ==== script for Fig.12 in CM4 paper
figure
vaxis = [1830, 2020, -0.55, 1.3];
%set(gca, 'FontWeight','bold','FontSize',20)
pos1 = [0.24 0.55 0.45 0.45];
subplot('Position',pos1)
patch('Xdata', tt_globe5, 'Ydata', yy_globe5, 'FaceColor', 'r');
hold
plot(PIyear5,  n_tas_PI_ann_globe5, '-b','LineWidth', 2) 
plot(giss_year5, n_giss_globe_ann5, '-k','LineWidth', 3) 
plot([0,3000], [0,0],'-k','LineWidth', 1)
text( 1840,1.1, '(a)  Globe','FontWeight','bold','FontSize',14)

axis(vaxis);
%xlabel('Year','FontWeight','bold','FontSize',20);
ylabel('Temperature Anomalies (K)','FontWeight','bold','FontSize',12);
set(gca,'XTick',[1840 1860 1880 1900 1920 1940 1960 1980 2000 2020]);
%set(gca,'XMinorTick','on')
%set(gca,'YMinorTick','on')

legend( {'spread of three CM4.0 historical','CM4.0 piControl','GISTEMP v4'}, 'Position',[0.795 0.75 0.1 0.2],'FontWeight','bold','FontSize',10) 
legend boxoff

pos2 = [0.069 0.065 0.45 0.45];
subplot('Position',pos2)
patch('Xdata', tt_sh5, 'Ydata', yy_sh5, 'FaceColor', 'r');
hold
plot(PIyear5,  n_tas_PI_ann_sh5, '-b','LineWidth', 2) 
plot(giss_year5, n_giss_sh_ann5, '-k','LineWidth', 3) 
plot([0,3000], [0,0],'-k','LineWidth', 1)
axis(vaxis);
xlabel('Year','FontWeight','bold','FontSize',13);
ylabel('Temperature Anomalies (K)','FontWeight','bold','FontSize',12);
text( 1840,1.1, '(b)  Southern Hemisphere','FontWeight','bold','FontSize',14)
set(gca,'XTick',[1840 1860 1880 1900 1920 1940 1960 1980 2000 2020]);

pos3 = [0.539 0.065 0.45 0.45];
subplot('Position',pos3)
patch('Xdata', tt_nh5, 'Ydata', yy_nh5, 'FaceColor', 'r');
hold
plot(PIyear5,  n_tas_PI_ann_nh5, '-b','LineWidth', 2) 
plot(giss_year5, n_giss_nh_ann5, '-k','LineWidth', 3) 
plot([0,3000], [0,0],'-k','LineWidth', 1)
axis(vaxis);
xlabel('Year','FontWeight','bold','FontSize',13);
%ylabel('Temperature Anomalies (K)','FontWeight','bold','FontSize',16);
text( 1840,1.1, '(c)  Northern Hemisphere','FontWeight','bold','FontSize',14)
set(gca,'XTick',[1840 1860 1880 1900 1920 1940 1960 1980 2000 2020]);
print('tas_tos_global_nh_sh','-depsc') 


% ==== end of h1g added, 2019-07-02 ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




