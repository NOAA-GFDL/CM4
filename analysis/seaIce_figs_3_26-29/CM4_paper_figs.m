function CM4_paper_figs(fig_dir,work_dir)

close all

addpath('./utils/export_fig/')
addpath('./utils/')
addpath('./utils/m_map/')

month_str={'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'};

%%%%%%LOAD MODEL DATA

%%CM4 Hist Ensemble Member 1
model='CM4';
run='CM4_historical';

dataDir=[work_dir,'/raw_data/',model,'/',run];

var=['siconc'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
SIC_nh=data;
w_nh=w;
x_nh=x;
y_nh=y;
ifXY_nh=ifXY;

var=['siconc_sh'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
SIC_sh=data;
w_sh=w;
x_sh=x;
y_sh=y;
ifXY_sh=ifXY;

var=['sithick'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
SIT_nh=data;

var=['sithick_sh'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
SIT_sh=data;

[SIC_clim_nh_E1 SIC_std_nh_E1 SIT_clim_nh_E1 SIT_std_nh_E1 extent_nh_E1 extent_clim_nh_E1 extent_anom_nh_E1 volume_nh_E1 volume_clim_nh_E1 volume_anom_nh_E1]=computeSeaIceDiagnostics(SIC_nh,SIT_nh,w_nh);
[SIC_clim_sh_E1 SIC_std_sh_E1 SIT_clim_sh_E1 SIT_std_sh_E1 extent_sh_E1 extent_clim_sh_E1 extent_anom_sh_E1 volume_sh_E1 volume_clim_sh_E1 volume_anom_sh_E1]=computeSeaIceDiagnostics(SIC_sh,SIT_sh,w_sh);

time_hist=1950+1/24:1/12:2015-1/24;
hist_years=1950:1:2014;
nT_hist=length(time_hist);

%define common period of 1979-2015 for computing climatology
ind1=(1979-1950)*12+1;
ind2=(2015-1950)*12;
[SIC_clim_common_nh_E1 SIC_std_common_nh_E1 SIT_clim_common_nh_E1 SIT_std_common_nh_E1 extent_common_nh_E1 extent_clim_common_nh_E1 extent_anom_common_nh_E1 volume_common_nh_E1 volume_clim_common_nh_E1 volume_anom_common_nh_E1]=computeSeaIceDiagnostics(SIC_nh(:,ind1:ind2),SIT_nh(:,ind1:ind2),w_nh);
[SIC_clim_common_sh_E1 SIC_std_common_sh_E1 SIT_clim_common_sh_E1 SIT_std_common_sh_E1 extent_common_sh_E1 extent_clim_common_sh_E1 extent_anom_common_sh_E1 volume_common_sh_E1 volume_clim_common_sh_E1 volume_anom_common_sh_E1]=computeSeaIceDiagnostics(SIC_sh(:,ind1:ind2),SIT_sh(:,ind1:ind2),w_sh);

%%CM4 Hist Ensemble Member 2
model='CM4';
run='CM4_historical_noBling_Run895Rst290';

dataDir=[work_dir,'/raw_data/',model,'/',run];

var=['siconc'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
SIC_nh=data;
w_nh=w;
x_nh=x;
y_nh=y;
ifXY_nh=ifXY;

var=['siconc_sh'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
SIC_sh=data;
w_sh=w;
x_sh=x;
y_sh=y;
ifXY_sh=ifXY;

var=['sithick'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
SIT_nh=data;

var=['sithick_sh'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
SIT_sh=data;

[SIC_clim_nh_E2 SIC_std_nh_E2 SIT_clim_nh_E2 SIT_std_nh_E2 extent_nh_E2 extent_clim_nh_E2 extent_anom_nh_E2 volume_nh_E2 volume_clim_nh_E2 volume_anom_nh_E2]=computeSeaIceDiagnostics(SIC_nh,SIT_nh,w_nh);
[SIC_clim_sh_E2 SIC_std_sh_E2 SIT_clim_sh_E2 SIT_std_sh_E2 extent_sh_E2 extent_clim_sh_E2 extent_anom_sh_E2 volume_sh_E2 volume_clim_sh_E2 volume_anom_sh_E2]=computeSeaIceDiagnostics(SIC_sh,SIT_sh,w_sh);

time_hist=1950+1/24:1/12:2015-1/24;
hist_years=1950:1:2014;
nT_hist=length(time_hist);

%define common period of 1979-2015 for computing climatology
ind1=(1979-1950)*12+1;
ind2=(2015-1950)*12;
[SIC_clim_common_nh_E2 SIC_std_common_nh_E2 SIT_clim_common_nh_E2 SIT_std_common_nh_E2 extent_common_nh_E2 extent_clim_common_nh_E2 extent_anom_common_nh_E2 volume_common_nh_E2 volume_clim_common_nh_E2 volume_anom_common_nh_E2]=computeSeaIceDiagnostics(SIC_nh(:,ind1:ind2),SIT_nh(:,ind1:ind2),w_nh);
[SIC_clim_common_sh_E2 SIC_std_common_sh_E2 SIT_clim_common_sh_E2 SIT_std_common_sh_E2 extent_common_sh_E2 extent_clim_common_sh_E2 extent_anom_common_sh_E2 volume_common_sh_E2 volume_clim_common_sh_E2 volume_anom_common_sh_E2]=computeSeaIceDiagnostics(SIC_sh(:,ind1:ind2),SIT_sh(:,ind1:ind2),w_sh);

%%CM4 Hist Ensemble Member 3
model='CM4';
run='CM4_historical_noBling_Run895Rst332';

dataDir=[work_dir,'/raw_data/',model,'/',run];

var=['siconc'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
SIC_nh=data;
w_nh=w;
x_nh=x;
y_nh=y;
ifXY_nh=ifXY;

var=['siconc_sh'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
SIC_sh=data;
w_sh=w;
x_sh=x;
y_sh=y;
ifXY_sh=ifXY;

var=['sithick'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
SIT_nh=data;

var=['sithick_sh'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
SIT_sh=data;

[SIC_clim_nh_E3 SIC_std_nh_E3 SIT_clim_nh_E3 SIT_std_nh_E3 extent_nh_E3 extent_clim_nh_E3 extent_anom_nh_E3 volume_nh_E3 volume_clim_nh_E3 volume_anom_nh_E3]=computeSeaIceDiagnostics(SIC_nh,SIT_nh,w_nh);
[SIC_clim_sh_E3 SIC_std_sh_E3 SIT_clim_sh_E3 SIT_std_sh_E3 extent_sh_E3 extent_clim_sh_E3 extent_anom_sh_E3 volume_sh_E3 volume_clim_sh_E3 volume_anom_sh_E3]=computeSeaIceDiagnostics(SIC_sh,SIT_sh,w_sh);

time_hist=1950+1/24:1/12:2015-1/24;
hist_years=1950:1:2014;
nT_hist=length(time_hist);

%define common period of 1979-2015 for computing climatology
ind1=(1979-1950)*12+1;
ind2=(2015-1950)*12;
[SIC_clim_common_nh_E3 SIC_std_common_nh_E3 SIT_clim_common_nh_E3 SIT_std_common_nh_E3 extent_common_nh_E3 extent_clim_common_nh_E3 extent_anom_common_nh_E3 volume_common_nh_E3 volume_clim_common_nh_E3 volume_anom_common_nh_E3]=computeSeaIceDiagnostics(SIC_nh(:,ind1:ind2),SIT_nh(:,ind1:ind2),w_nh);
[SIC_clim_common_sh_E3 SIC_std_common_sh_E3 SIT_clim_common_sh_E3 SIT_std_common_sh_E3 extent_common_sh_E3 extent_clim_common_sh_E3 extent_anom_common_sh_E3 volume_common_sh_E3 volume_clim_common_sh_E3 volume_anom_common_sh_E3]=computeSeaIceDiagnostics(SIC_sh(:,ind1:ind2),SIT_sh(:,ind1:ind2),w_sh);

%%%%%%%LOAD OBS DATA

%p25 grid

%%load obs data

obsdir=['/work/mib/raw_data/NSIDC/NASATEAM_merged_v1.1_nrt'];
load([obsdir,'/sic_nh_regrid_GFDL_CM4.mat'],'data','x','y','ifXY','w','area_weights')
SIC_obs_nh=data/100;
w_obs_nh=w;

SIT_obs_nh=zeros(size(SIC_obs_nh));

obsdir=['/work/mib/raw_data/NSIDC/NASATEAM_merged_v1.1_nrt_sh'];
load([obsdir,'/sic_sh_regrid_GFDL_CM4.mat'],'data','x','y','ifXY','w','area_weights')
SIC_obs_sh=data/100;
w_obs_sh=w;

SIT_obs_sh=zeros(size(SIC_obs_sh));

[SIC_clim_obs_nh SIC_std_obs_nh SIT_clim_obs_nh SIT_std_obs_nh extent_obs_nh extent_clim_obs_nh extent_anom_obs_nh volume_obs_nh volume_clim_obs_nh volume_anom_obs_nh]=computeSeaIceDiagnostics(SIC_obs_nh,SIT_obs_nh,w_obs_nh);
[SIC_clim_obs_sh SIC_std_obs_sh SIT_clim_obs_sh SIT_std_obs_sh extent_obs_sh extent_clim_obs_sh extent_anom_obs_sh volume_obs_sh volume_clim_obs_sh volume_anom_obs_sh]=computeSeaIceDiagnostics(SIC_obs_sh,SIT_obs_sh,w_obs_sh);

%define common period of 1979-2014 for computing climatology
ind1=(1979-1979)*12+1;
ind2=(2015-1979)*12;

[SIC_clim_obs_common_nh SIC_std_obs_common_nh SIT_clim_obs_common_nh SIT_std_obs_common_nh extent_obs_common_nh extent_clim_obs_common_nh extent_anom_obs_common_nh volume_obs_common_nh volume_clim_obs_common_nh volume_anom_obs_common_nh]=computeSeaIceDiagnostics(SIC_obs_nh(:,ind1:ind2),SIT_obs_nh(:,ind1:ind2),w_obs_nh);
[SIC_clim_obs_common_sh SIC_std_obs_common_sh SIT_clim_obs_common_sh SIT_std_obs_common_sh extent_obs_common_sh extent_clim_obs_common_sh extent_anom_obs_common_sh volume_obs_common_sh volume_clim_obs_common_sh volume_anom_obs_common_sh]=computeSeaIceDiagnostics(SIC_obs_sh(:,ind1:ind2),SIT_obs_sh(:,ind1:ind2),w_obs_sh);

%%nsidc grid
time_obs=1979+1/24:1/12:2018-1/24;
obs_years=1979:2017;

obsdir=['/work/mib/raw_data/NSIDC/NASATEAM_merged_v1.1_nrt_2017'];
load([obsdir,'/sic.mat'],'data','x','y','ifXY','w')
SIC_NSIDCgrid_nh=data/100;
w_NSIDC_nh=w;
x_NSIDC=x;
y_NSIDC=y;
ifXY_NSIDC=ifXY;

clearvars w

SIT_NSIDCgrid_nh=zeros(size(SIC_NSIDCgrid_nh));

[SIC_clim_NSIDCgrid_nh SIC_std_NSIDCgrid_nh SIT_clim_NSIDCgrid_nh SIT_std_NSIDCgrid_nh extent_NSIDCgrid_nh extent_clim_NSIDCgrid_nh extent_anom_NSIDCgrid_nh volume_NSIDCgrid_nh volume_clim_NSIDCgrid_nh volume_anom_NSIDCgrid_nh]=computeSeaIceDiagnostics(SIC_NSIDCgrid_nh(:,:),SIT_NSIDCgrid_nh(:,:),w_NSIDC_nh);

%define common period of 1979-2014 for computing climatology
ind1=(1979-1979)*12+1;
ind2=(2015-1979)*12;

[SIC_clim_NSIDCgrid_common_nh SIC_std_NSIDCgrid_common_nh SIT_clim_NSIDCgrid_common_nh SIT_std_NSIDCgrid_common_nh extent_NSIDCgrid_common_nh extent_clim_NSIDCgrid_common_nh extent_anom_NSIDCgrid_common_nh volume_NSIDCgrid_common_nh volume_clim_NSIDCgrid_common_nh volume_anom_NSIDCgrid_common_nh]=computeSeaIceDiagnostics(SIC_NSIDCgrid_nh(:,ind1:ind2),SIT_NSIDCgrid_nh(:,ind1:ind2),w_NSIDC_nh);

obsdir=['/work/mib/raw_data/NSIDC/NASATEAM_merged_v1.1_nrt_2017_sh'];
load([obsdir,'/sic.mat'],'data','x','y','ifXY','w')
SIC_NSIDCgrid_sh=data/100;
w_NSIDC_sh=w;
x_NSIDC=x;
y_NSIDC=y;
ifXY_NSIDC_sh=ifXY;

SIT_NSIDCgrid_sh=zeros(size(SIC_NSIDCgrid_sh));

[SIC_clim_NSIDCgrid_sh SIC_std_NSIDCgrid_sh SIT_clim_NSIDCgrid_sh SIT_std_NSIDCgrid_sh extent_NSIDCgrid_sh extent_clim_NSIDCgrid_sh extent_anom_NSIDCgrid_sh volume_NSIDCgrid_sh volume_clim_NSIDCgrid_sh volume_anom_NSIDCgrid_sh]=computeSeaIceDiagnostics(SIC_NSIDCgrid_sh(:,:),SIT_NSIDCgrid_sh(:,:),w_NSIDC_sh);

%define common period of 1979-2014 for computing climatology
ind1=(1979-1979)*12+1;
ind2=(2015-1979)*12;

[SIC_clim_NSIDCgrid_common_sh SIC_std_NSIDCgrid_common_sh SIT_clim_NSIDCgrid_common_sh SIT_std_NSIDCgrid_common_sh extent_NSIDCgrid_common_sh extent_clim_NSIDCgrid_common_sh extent_anom_NSIDCgrid_common_sh volume_NSIDCgrid_common_sh volume_clim_NSIDCgrid_common_sh volume_anom_NSIDCgrid_common_sh]=computeSeaIceDiagnostics(SIC_NSIDCgrid_sh(:,ind1:ind2),SIT_NSIDCgrid_sh(:,ind1:ind2),w_NSIDC_sh);

%%%%%%%%MAKE FIGURES

%%%FIGURE 1: NH/SH SIE

    nTileX=2;
    nTileY = 1;
        
    figWidth = 9.5; %6.5 inches = 39 picas, for a two column figure
    deltaX   = 0.6;
    deltaX2  = 0.2;
    deltaY   = 0.5;
    deltaY2  = 0.3;
    gapX     = 0.4;
    gapY     = 0.6;
    panel_scale=0.8;

ax=make_axes(nTileX,nTileY,figWidth,deltaX,deltaX2,deltaY,deltaY2,gapX,gapY,panel_scale);

axes(ax(1,1))

fact=1e-12;
L=1;

hold on
plot(1:12,fact*extent_clim_NSIDCgrid_common_nh,'k','LineWidth',2)
plot(1:12,fact*(extent_clim_common_nh_E1+extent_clim_common_nh_E2+extent_clim_common_nh_E3)/3,'r','LineWidth',2)
plot(1:12,fact*extent_clim_common_nh_E1,'r','LineWidth',L)
plot(1:12,fact*extent_clim_common_nh_E2,'r','LineWidth',L)
plot(1:12,fact*extent_clim_common_nh_E3,'r','LineWidth',L)
plot(1:12,fact*(extent_clim_common_nh_E1+extent_clim_common_nh_E2+extent_clim_common_nh_E3)/3,'r','LineWidth',2)
xlabel('Month')
ylabel('SIE (Million km^2)')
xlim([0.5 12.5])
ylim([6 17])
grid on
title('Arctic SIE Climatology (1979-2014)')
set(gca,'xTick',1:12,'xTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'} )
legend('NSIDC Obs','GFDL CM4.0','location','SouthWest')
text(-0.5,17.6,'A','FontWeight','bold')

EnsMean=(extent_clim_common_nh_E1+extent_clim_common_nh_E2+extent_clim_common_nh_E3)/3;
disp('Arctic RMSE SIE')
RMS_nh=sqrt(sum((EnsMean-extent_clim_NSIDCgrid_common_nh).^2)/12)

axes(ax(2,1))

hold on
plot(1:12,fact*extent_clim_NSIDCgrid_common_sh,'k','LineWidth',2)
plot(1:12,fact*(extent_clim_common_sh_E1+extent_clim_common_sh_E2+extent_clim_common_sh_E3)/3,'r','LineWidth',2)
plot(1:12,fact*extent_clim_common_sh_E1,'r','LineWidth',L)
plot(1:12,fact*extent_clim_common_sh_E2,'r','LineWidth',L)
plot(1:12,fact*extent_clim_common_sh_E3,'r','LineWidth',L)
plot(1:12,fact*(extent_clim_common_sh_E1+extent_clim_common_sh_E2+extent_clim_common_sh_E3)/3,'r','LineWidth',2)
xlabel('Month')
xlim([0.5 12.5])
ylim([0 23])
grid on
title('Antarctic SIE Climatology (1979-2014)')
set(gca,'xTick',1:12,'xTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'} )
legend('NSIDC Obs','GFDL CM4.0','location','SouthEast')
text(-0.2,24,'B','FontWeight','bold')

EnsMean=(extent_clim_common_sh_E1+extent_clim_common_sh_E2+extent_clim_common_sh_E3)/3;
disp('Antarctic RMSE SIE')
RMS_sh=sqrt(sum((EnsMean-extent_clim_NSIDCgrid_common_sh).^2)/12)

set(gcf,'Color','w')
filename=[fig_dir,'/SIE_clim.pdf'];
export_fig(gcf,filename,'-a1','-pdf');

%%%FIGURE 2: NH/SH SIE timeseries

    nTileX=2;
    nTileY = 2;
        
    figWidth = 9.5; %6.5 inches = 39 picas, for a two column figure
    deltaX   = 0.7;
    deltaX2  = 0.2;
    deltaY   = 0.5;
    deltaY2  = 0.3;
    gapX     = 0.4;
    gapY     = 0.8;
    panel_scale=0.8;

ax=make_axes(nTileX,nTileY,figWidth,deltaX,deltaX2,deltaY,deltaY2,gapX,gapY,panel_scale);


fact=1e-12;
L=1;

size(extent_NSIDCgrid_nh)
size(extent_NSIDCgrid_sh)
size(obs_years)

axes(ax(1,1))
month=3;
hold on

plot(obs_years,fact*extent_NSIDCgrid_nh(month:12:end),'k','LineWidth',2)
plot(hist_years,fact*(extent_nh_E1(month:12:end)+extent_nh_E2(month:12:end)+extent_nh_E3(month:12:end))/3,'r','LineWidth',2)
plot(hist_years,fact*extent_nh_E1(month:12:end),'r','LineWidth',L)
plot(hist_years,fact*extent_nh_E2(month:12:end),'r','LineWidth',L)
plot(hist_years,fact*extent_nh_E3(month:12:end),'r','LineWidth',L)
grid on
legend('NSIDC Obs','GFDL CM4.0','location','SouthWest')
xlim([1950 2018])
ylim([14.3 17.5])
xlabel('Time (years)')
ylabel('SIE (Million km^2)')
title('March Arctic SIE')
text(1943,17.8,'A','FontWeight','bold')

EnsMean=(extent_nh_E1(month:12:end)+extent_nh_E2(month:12:end)+extent_nh_E3(month:12:end))/3;


axes(ax(2,1))
month=9;
hold on
plot(obs_years,fact*extent_NSIDCgrid_nh(month:12:end),'k','LineWidth',2)
plot(hist_years,fact*(extent_nh_E1(month:12:end)+extent_nh_E2(month:12:end)+extent_nh_E3(month:12:end))/3,'r','LineWidth',2)
plot(hist_years,fact*extent_nh_E1(month:12:end),'r','LineWidth',L)
plot(hist_years,fact*extent_nh_E2(month:12:end),'r','LineWidth',L)
plot(hist_years,fact*extent_nh_E3(month:12:end),'r','LineWidth',L)
grid on
xlim([1950 2018])
ylim([3.5 9.5])
xlabel('Time (years)')
title('September Arctic SIE')
text(1947,10,'B','FontWeight','bold')

EnsMean=(extent_nh_E1(month:12:end)+extent_nh_E2(month:12:end)+extent_nh_E3(month:12:end))/3;
disp('Arctic Sept Ens Mean')

axes(ax(1,2))
month=3;
hold on
plot(obs_years,fact*extent_NSIDCgrid_sh(month:12:end),'k','LineWidth',2)
plot(hist_years,fact*(extent_sh_E1(month:12:end)+extent_sh_E2(month:12:end)+extent_sh_E3(month:12:end))/3,'r','LineWidth',2)
plot(hist_years,fact*extent_sh_E1(month:12:end),'r','LineWidth',L)
plot(hist_years,fact*extent_sh_E2(month:12:end),'r','LineWidth',L)
plot(hist_years,fact*extent_sh_E3(month:12:end),'r','LineWidth',L)
plot(hist_years,fact*(extent_sh_E1(month:12:end)+extent_sh_E2(month:12:end)+extent_sh_E3(month:12:end))/3,'r','LineWidth',2)
grid on
xlim([1950 2018])
xlabel('Time (years)')
ylabel('SIE (Million km^2)')
title('March Antarctic SIE')
text(1943,6.4,'C','FontWeight','bold')

axes(ax(2,2))
month=9;
hold on
plot(obs_years,fact*extent_NSIDCgrid_sh(month:12:end),'k','LineWidth',2)
plot(hist_years,fact*(extent_sh_E1(month:12:end)+extent_sh_E2(month:12:end)+extent_sh_E3(month:12:end))/3,'r','LineWidth',2)
plot(hist_years,fact*extent_sh_E1(month:12:end),'r','LineWidth',L)
plot(hist_years,fact*extent_sh_E2(month:12:end),'r','LineWidth',L)
plot(hist_years,fact*extent_sh_E3(month:12:end),'r','LineWidth',L)
plot(hist_years,fact*(extent_sh_E1(month:12:end)+extent_sh_E2(month:12:end)+extent_sh_E3(month:12:end))/3,'r','LineWidth',2)
grid on
xlim([1950 2018])
xlabel('Time (years)')
title('September Antarctic SIE')
text(1947,24.6,'D','FontWeight','bold')

set(gcf,'Color','w')
filename=[fig_dir,'/SIE_interannnual.pdf'];
export_fig(gcf,filename,'-a1','-pdf');

%%%FIGURE 3: SIC biases

    nTileX=2;
    nTileY = 2;
        
    figWidth = 6.5; %6.5 inches = 39 picas, for a two column figure
    deltaX   = 0.4;
    deltaX2  = 0.02;
    deltaY   = 0.5;
    deltaY2  = 0.3;
    gapX     = 0.02;
    gapY     = 0.1;
    panel_scale=1;

ax=make_axes(nTileX,nTileY,figWidth,deltaX,deltaX2,deltaY,deltaY2,gapX,gapY,panel_scale);

SIC_clim_common_nh_EnsMean=(SIC_clim_common_nh_E1+SIC_clim_common_nh_E2+SIC_clim_common_nh_E3)/3;
SIC_clim_common_sh_EnsMean=(SIC_clim_common_sh_E1+SIC_clim_common_sh_E2+SIC_clim_common_sh_E3)/3;

axes(ax(1,1))
tmp_obs_p25=SIC_clim_obs_common_nh(:,3);

tmp=SIC_clim_common_nh_EnsMean(:,3);
make_polar_plot_pcolor(x_nh,y_nh,ifXY_nh,tmp-tmp_obs_p25,45,'',0,-1,1);
make_polar_contour_plot(x_nh,y_nh,ifXY_nh,tmp_obs_p25,45,'','k',1)
text(-0.75,0.85,'A','FontWeight','bold')
text(-0.45,0.91,'March SIC Bias','FontWeight','bold','FontSize',16)

axes(ax(2,1))
tmp_obs_p25=SIC_clim_obs_common_nh(:,9);

tmp=SIC_clim_common_nh_EnsMean(:,9);
make_polar_plot_pcolor(x_nh,y_nh,ifXY_nh,tmp-tmp_obs_p25,30,'',0,-1,1);
make_polar_contour_plot(x_nh,y_nh,ifXY_nh,tmp_obs_p25,30,'','k',1)
text(-0.47,0.58,'B','FontWeight','bold')
text(-0.25,0.6,'Sept SIC Bias','FontWeight','bold','FontSize',16)

axes(ax(1,2))
tmp_obs_p25=SIC_clim_obs_common_sh(:,3);

tmp=SIC_clim_common_sh_EnsMean(:,3);
make_polar_plot_pcolor_sh(x_sh,y_sh,ifXY_sh,tmp-tmp_obs_p25,30,'',0,-1,1);
make_polar_contour_plot_sh(x_sh,y_sh,ifXY_sh,tmp_obs_p25,30,'','k',1)
text(-0.47,0.58,'C','FontWeight','bold')

                 pos=get(gca,'position');

                 c=colorbar('horizontal');

                 set(gca,'position',pos)

                 x_pos=get(c, 'Position');
                 y_pos=[x_pos(1)+0.02 x_pos(2)-0.07 2.1*x_pos(3) 0.7*x_pos(4)];
                 set(c,'Position',y_pos)


axes(ax(2,2))
tmp_obs_p25=SIC_clim_obs_common_sh(:,9);

tmp=SIC_clim_common_sh_EnsMean(:,9);
make_polar_plot_pcolor_sh(x_sh,y_sh,ifXY_sh,tmp-tmp_obs_p25,40,'',0,-1,1);
make_polar_contour_plot_sh(x_sh,y_sh,ifXY_sh,tmp_obs_p25,40,'','k',1)
text(-0.75,0.85,'D','FontWeight','bold')

set(gcf,'Color','w')
filename=[fig_dir,'/SIC_bias'];
export_fig(gcf,filename,'-a1','-png','-r600','-painters');


