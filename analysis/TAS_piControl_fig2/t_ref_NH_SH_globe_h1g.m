clear all
nc64startup;

nc_96_static  = netcdf('/net2/h1g/CM4_paper/cm4_data/c96/atmos.static.nc', 'nowrite');
nc_PI_tas = netcdf('/net2/h1g/CM4_paper/cm4_data/PI_C/TAS/tas.nc', 'nowrite');

rad = pi/180;

v_lat   = ncvar('lat'    ,   nc_96_static);
v_lon   = ncvar('lon'    ,   nc_96_static);
v_lat_bnds  = ncvar('lat_bnds'    ,   nc_96_static);
v_lon_bnds  = ncvar('lon_bnds'    ,   nc_96_static);
v_land  = ncvar('land_mask' ,   nc_96_static);
lat   = v_lat(:);
lon   = v_lon(:);
lat_bnds   = v_lat_bnds(:,:);
lon_bnds   = v_lon_bnds(:,:);
land  = v_land(:,:);
land = ceil(land);
ocean = floor(1.-land);
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

% === change by h1g, 2019-07-03===
%lim_nh_ex = [ 35, 90, -180, 180];
lim_nh_ex = [ 45, 90, -180, 180];
% === end change by h1g, 2019-07-03===

lim_trop  = [-30, 30, -180, 180];
lim_rest  = [-45, 90, -180, 180];
lim_rest1 = [-55, 90, -180, 180];
lim_sp    = [-90,-55, -180, 180];
lim_sh_md = [-55,-30, -180, 180];
lim_ant   = [-90,-65, -180, 180];

PIyears = 650;
PIyear = [1:PIyears] + 1850 - 250;


% === change by h1g, 2019-07-03===
%ya = 250;
%yb = 405;
ya = 151;
yb = 650;
% === end change by h1g, 2019-07-03===

ya2 = 250+40;
yb2 = 405+40;

ya3 = 250+82;
yb3 = 405+82;

tas_PI_ann_globe = series_globe(PIyears, 0, lat, lon, wt, nc_PI_tas, 'tas', lim_globe, 'ann', nfilt);
tas_PI_ann_nh = series_globe(PIyears, 0, lat, lon, wt, nc_PI_tas, 'tas', lim_nh, 'ann', nfilt);
tas_PI_ann_sh = series_globe(PIyears, 0, lat, lon, wt, nc_PI_tas, 'tas', lim_sh, 'ann', nfilt);
tas_PI_ann_nh_ex = series_globe(PIyears, 0, lat, lon, wt, nc_PI_tas, 'tas', lim_nh_ex, 'ann', nfilt);
tas_PI_ann_sh_ex = series_globe(PIyears, 0, lat, lon, wt, nc_PI_tas, 'tas', lim_sh_ex, 'ann', nfilt);
tas_PI_ann_trop = series_globe(PIyears, 0, lat, lon, wt, nc_PI_tas, 'tas', lim_trop, 'ann', nfilt);


n_tas_PI_ann_globe = tas_PI_ann_globe - mean(tas_PI_ann_globe(ya:yb));
n_tas_PI_ann_nh = tas_PI_ann_nh - mean(tas_PI_ann_nh(ya:yb));
n_tas_PI_ann_sh = tas_PI_ann_sh - mean(tas_PI_ann_sh(ya:yb));
n_tas_PI_ann_nh_ex = tas_PI_ann_nh_ex - mean(tas_PI_ann_nh_ex(ya:yb));
n_tas_PI_ann_sh_ex = tas_PI_ann_sh_ex - mean(tas_PI_ann_sh_ex(ya:yb));
n_tas_PI_ann_trop = tas_PI_ann_trop - mean(tas_PI_ann_trop(ya:yb));

start_year = 150;
start_year1 = start_year + 1;
end_year = PIyears;

kk=10;
kk1 = kk - 1;
for k = 1:PIyears-kk
  n_tas_PI_ann_globe10(k) = mean(n_tas_PI_ann_globe(k:k+kk));
  n_tas_PI_ann_nh_ex10(k) = mean(n_tas_PI_ann_nh_ex(k:k+kk));
  n_tas_PI_ann_sh_ex10(k) = mean(n_tas_PI_ann_sh_ex(k:k+kk));
  n_tas_PI_ann_trop10(k) = mean(n_tas_PI_ann_trop(k:k+kk));
  PIyear10(k) = mean(PIyear(k:k+kk)); 
end

figure
set(gca, 'FontWeight','bold','FontSize',18)
set(gca,'XMinorTick','on','YMinorTick','on')

vaxis = [-20, 690, -0.8, 1.1];
plot(PIyear10 - 1600, n_tas_PI_ann_globe10, '-k','LineWidth', 2)
hold
plot(PIyear10 - 1600, n_tas_PI_ann_nh_ex10, '-b','LineWidth', 1)
plot(PIyear10 - 1600, n_tas_PI_ann_sh_ex10, '-r','LineWidth', 1)
plot([-3000, 3000],[0,0],  '-k','LineWidth', 1)
plot([151, 151],[-2,2],  '-k','LineWidth', 1)
plot([650, 650],[-2,2],  '-k','LineWidth', 1)
axis(vaxis);
xlabel('Year','FontWeight','bold','FontSize',18);
ylabel('Temperature Anomalies (K)','FontWeight','bold','FontSize',18);
legend( {'Globe', '45N-90N','45S-90S'} , 'Position', [0.68 0.18 0.10 0.07],'FontWeight','bold','FontSize',13) 
legend boxoff

% Arrow with two head at both end and text between
x = [0.760 0.860];    % adjust length and location of arrow 
y = [0.85 0.85];    
Xadj = 1.177;      % adjust location of left arrow starting point (the sum of this with 'x' should not be negative)
annotation('textarrow',x,y,'String',' CM4.0 piControl years 1-500 ','FontWeight','bold','FontSize',14,'Linewidth',2)
annotation('textarrow',-x+Xadj,y,'String','','FontSize',14,'Linewidth',2)

     
print('tas_piControl_global_nh45_sh45','-depsc') 


