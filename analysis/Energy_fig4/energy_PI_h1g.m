
nc64startup;

nc_static  = netcdf('/net2/h1g/CM4_paper/cm4_data/c96/atmos.static.nc', 'nowrite');

nc_PI_TAS = netcdf('/net2/h1g/CM4_paper/cm4_data/PI_C/TAS/tas.nc', 'nowrite');
nc_PI_RLUT = netcdf('/net2/h1g/CM4_paper/cm4_data/PI_C/RLUT/rlut.nc', 'nowrite');
nc_PI_RSUT = netcdf('/net2/h1g/CM4_paper/cm4_data/PI_C/RSUT/rsut.nc', 'nowrite');
nc_PI_RSDT = netcdf('/net2/h1g/CM4_paper/cm4_data/PI_C/RSDT/rsdt.nc', 'nowrite');
nc_PI_RLUS = netcdf('/net2/h1g/CM4_paper/cm4_data/PI_C/RLUS/rlus.nc', 'nowrite');
nc_PI_RSUS = netcdf('/net2/h1g/CM4_paper/cm4_data/PI_C/RSUS/rsus.nc', 'nowrite');
nc_PI_RLDS = netcdf('/net2/h1g/CM4_paper/cm4_data/PI_C/RLDS/rlds.nc', 'nowrite');
nc_PI_RSDS = netcdf('/net2/h1g/CM4_paper/cm4_data/PI_C/RSDS/rsds.nc', 'nowrite');
nc_PI_HFLS = netcdf('/net2/h1g/CM4_paper/cm4_data/PI_C/HFLS/hfls.nc', 'nowrite');
nc_PI_HFSS = netcdf('/net2/h1g/CM4_paper/cm4_data/PI_C/HFSS/hfss.nc', 'nowrite');
nc_PI_PRSN = netcdf('/net2/h1g/CM4_paper/cm4_data/PI_C/PRSN/prsn.nc', 'nowrite');
nc_PI_PR = netcdf('/net2/h1g/CM4_paper/cm4_data/PI_C/PR/pr.nc', 'nowrite');


LV = 2.5e06;
LF = 3.34e05;


rad = pi/180;

v_lat    = ncvar('lat'    ,   nc_static);
v_lon    = ncvar('lon'    ,   nc_static);
v_lat_bnds   = ncvar('lat_bnds'    ,   nc_static);
v_lon_bnds   = ncvar('lon_bnds'    ,   nc_static);
v_land   = ncvar('land_mask' ,   nc_static);
lat    = v_lat (:);
lon    = v_lon (:);
lat_bnds    = v_lat_bnds (:,:);
lon_bnds    = v_lon_bnds (:,:);
land  = v_land (:,:);
land  = ceil(land );
ocean  = floor(1.-land );
nlat  = size(lat ,1);
nlon   = size(lon ,1);
for j = 1:nlat 
  wt (j) = 0.5*(sin(rad*lat_bnds (j,2)) - sin(rad*lat_bnds (j,1)));
end

nfilt = 1;
lim = [-90, 90, -180, 180];
seas = 'ann';



PIyears = 650;
PIyear = [1:PIyears] - 250 + 1850;

ya = 250;
yb = 400;

tas_PI_ann_globe = series_globe(PIyears, 0, lat, lon, wt, nc_PI_TAS, 'tas', lim, 'ann', nfilt);
rlut_PI_ann_globe = series_globe(PIyears, 0, lat, lon, wt, nc_PI_RLUT, 'rlut', lim, 'ann', nfilt);
rsut_PI_ann_globe = series_globe(PIyears, 0, lat, lon, wt, nc_PI_RSUT, 'rsut', lim, 'ann', nfilt);
rsdt_PI_ann_globe = series_globe(PIyears, 0, lat, lon, wt, nc_PI_RSDT, 'rsdt', lim, 'ann', nfilt);
rlus_PI_ann_globe = series_globe(PIyears, 0, lat, lon, wt, nc_PI_RLUS, 'rlus', lim, 'ann', nfilt);
rsus_PI_ann_globe = series_globe(PIyears, 0, lat, lon, wt, nc_PI_RSUS, 'rsus', lim, 'ann', nfilt);
rlds_PI_ann_globe = series_globe(PIyears, 0, lat, lon, wt, nc_PI_RLDS, 'rlds', lim, 'ann', nfilt);
rsds_PI_ann_globe = series_globe(PIyears, 0, lat, lon, wt, nc_PI_RSDS, 'rsds', lim, 'ann', nfilt);
hfls_PI_ann_globe = series_globe(PIyears, 0, lat, lon, wt, nc_PI_HFLS, 'hfls', lim, 'ann', nfilt);
hfss_PI_ann_globe = series_globe(PIyears, 0, lat, lon, wt, nc_PI_HFSS, 'hfss', lim, 'ann', nfilt);
prsn_PI_ann_globe = series_globe(PIyears, 0, lat, lon, wt, nc_PI_PRSN, 'prsn', lim, 'ann', nfilt);
pr_PI_ann_globe = series_globe(PIyears, 0, lat, lon, wt, nc_PI_PR, 'pr', lim, 'ann', nfilt);


control_tas = mean(tas_PI_ann_globe(ya:yb));

net_surf_PI_ann_globe = rsds_PI_ann_globe + rlds_PI_ann_globe - rsus_PI_ann_globe - rlus_PI_ann_globe ...
    - hfls_PI_ann_globe - hfss_PI_ann_globe - LF*prsn_PI_ann_globe;
net_toa_PI_ann_globe = - rlut_PI_ann_globe - rsut_PI_ann_globe + rsdt_PI_ann_globe;
net_atm_PI_ann_globe = net_toa_PI_ann_globe - net_surf_PI_ann_globe;

mean_net_toa_PI_ann_globe = mean(net_toa_PI_ann_globe(ya:yb))
mean_net_surf_PI_ann_globe = mean(net_surf_PI_ann_globe(ya:yb))
mean_net_atm_PI_ann_globe = mean(net_atm_PI_ann_globe(ya:yb))

meantot_net_toa_PI_ann_globe = mean(net_toa_PI_ann_globe(151:PIyears))
meantot_net_surf_PI_ann_globe = mean(net_surf_PI_ann_globe(151:PIyears))
meantot_net_atm_PI_ann_globe = mean(net_atm_PI_ann_globe(151:PIyears))




for y = 1:PIyears-19
  net_toa_PI_ann_globe_20(y) = mean(net_toa_PI_ann_globe(y:y+19));
  net_surf_PI_ann_globe_20(y) = mean(net_surf_PI_ann_globe(y:y+19));
  net_atm_PI_ann_globe_20(y) = mean(net_atm_PI_ann_globe(y:y+19));
  tas_PI_ann_globe_20(y) = mean(tas_PI_ann_globe(y:y+19));
  PIyear_20(y) = mean(PIyear(y:y+19));
end


load thetao_PImat.dat

c = 3.925e3;
masso = 1.38179e+21;
rr = 6.4e+06;
area = 4*pi*rr*rr;
sec_year = 86400*365.14;


tt_PI = thetao_PImat;
th_PI = tt_PI(:,2);
OPIyears = size(th_PI,1);
OPIyear = [1:OPIyears] - 250 + 1850;
heat_PI = th_PI*c*masso/area;
for j = 1:OPIyears-1
  orate_PI(j) = (heat_PI(j+1)-heat_PI(j))/sec_year;
end
mean_orate_PI = mean(orate_PI)
OPIyears_20 = OPIyears - 1 - 19;
for j = 1:OPIyears_20
  orate_PI_20(j) = mean(orate_PI(j:j+19));
  OPIyear_20(j) = mean(OPIyear(j:j+19));
end

meantot_net_ocean_PI_ann_globe = mean(orate_PI(151:size(orate_PI,2)))
meantot_net_surf_ocean_PI_ann_globe = meantot_net_surf_PI_ann_globe - meantot_net_ocean_PI_ann_globe



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

% Plot the first two sets of data
ph1(1) = plot(PIyear_20 - 1750, net_toa_PI_ann_globe_20, '-r','LineWidth',1); 
hold on
ph1(2) = plot(PIyear_20 - 1750, net_surf_PI_ann_globe_20, '-k','LineWidth',1);
xlabel('piControl Year','FontSize',18);
ylabel('Net Downward Flux (W/m^2)','FontSize',18);

% Create the first legend
lh1 = legend({'TOA','SURF'}, 'Position', [0.2 0.75 0.10 0.2],'FontWeight','bold','FontSize',16);
set(lh1,'box','off')
lh1_position = get(lh1,'Position');
% Now set any axis properties that you want (Font, Ticks, etc)
fontSize = 16;
xlimits = [0, 500];
ylimits = [0, 0.50];
xticks  = [0, 100, 200, 300, 400, 500];
yticks  = [0, 0.1, 0.2, 0.3, 0.4, 0.5];
set(gca,'XLim',xlimits,'YLim',ylimits,'FontSize',fontSize,'XTick',xticks,'YTick',yticks,'YMinorTick','on');


% Name the first axis ax1 and create a second axis on top of the first
ax1 = gca;
pos= get(ax1,'Position')
ax2 = axes('Position',pos) ;
% Plot the second two sets of data on ax2
ph2(1) = plot(OPIyear_20 - 1750, orate_PI_20, '-b','LineWidth',1);
hold on
ph2(2) = plot(PIyear_20 - 1750, net_atm_PI_ann_globe_20, '-g','LineWidth',1); 
% Again set any axis properties that you want (Font, Ticks, etc) to make sure that the Legend fonts are the same
set(ax2,'XLim',xlimits,'YLim',ylimits,'FontSize',fontSize,'XTick',xticks,'YTick',yticks);
% Now, link the first axis to the second and make the second invisible
linkaxes([ax1 ax2],'xy');
set(ax2,'Color','none','Position',pos,'XTick',[],'YTick',[],'Box','off');
ax2.XAxis.Visible = 'off';


% Now make the second legend just above the first
lh2 = legend({'OCEAN', 'TOA minus SURF (0.08)'}, 'Position', [0.55 0.75 0.10 0.2],'FontWeight','bold','FontSize',16,'boxoff') ;
set(lh2,'box','off')

print -depsc atmosphere_ocean_balance_final_PI.eps
print -dpng  atmosphere_ocean_balance_final_PI.png



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure

% Plot the first two sets of data
 plot(PIyear_20(151:PIyears-19) - 1750, net_toa_PI_ann_globe_20(151:PIyears-19), '-r','LineWidth',1); 
hold on
 plot(PIyear_20(151:PIyears-19) - 1750, net_surf_PI_ann_globe_20(151:PIyears-19), '-k','LineWidth',1);
 plot(OPIyear_20(151:OPIyears_20) - 1750, orate_PI_20(151:OPIyears_20), '-b','LineWidth',1);
 plot(PIyear_20(151:PIyears-19) - 1750, net_atm_PI_ann_globe_20(151:PIyears-19), '-g','LineWidth',1); 
 fontSize = 16;
xlimits = [0, 500];
ylimits = [0, 0.50];
xticks  = [0, 100, 200, 300, 400, 500];
yticks  = [0, 0.1, 0.2, 0.3, 0.4, 0.5];
set(gca,'XLim',xlimits,'YLim',ylimits,'FontSize',fontSize,'XTick',xticks,'YTick',yticks,'YMinorTick','on');
xlabel('piControl Year','FontWeight','bold','FontSize',18);
ylabel('Net Downward Flux (W/m^2)','FontWeight','bold','FontSize',18);

text(50, 0.481,'_____','Color','red','FontWeight','bold','FontSize',14 );
text(105,0.47,'TOA', 'Color','k',    'FontWeight','bold','FontSize',14 );
 
text(50, 0.451,'_____','Color','k', 'FontWeight','bold','FontSize',14 );
text(105,0.44,'SURF', 'Color','k', 'FontWeight','bold','FontSize',14 ) ;

text(200,0.481,'_____','Color','b','FontWeight','bold','FontSize',14 );
text(255,0.47,'OCN',  'Color','k',    'FontWeight','bold','FontSize',14 ) ;
 
text(200, 0.451,'_____','Color','g', 'FontWeight','bold','FontSize',14 )
text(255,0.44,'TOA minus SURF (0.08)', 'Color','k', 'FontWeight','bold','FontSize',14 ) ;

print -depsc atmosphere_ocean_balance_piControl.eps

 
