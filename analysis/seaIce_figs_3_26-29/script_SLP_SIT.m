function script_SLP_SIT(fig_dir,work_dir)

close all

addpath('./utils/export_fig/')
addpath('./utils/')
addpath('./utils/m_map/')

load('cmap_red_seq.mat')
cmapSIT=cmap;

load('cmap_jet3.mat','cmap')
cmapSLP=cmap;

SLPyr1=2005;
SLPyr2=2014+1;

SITyr1=1998;
SITyr2=2007+1;

CM4SITyr1=2005;
CM4SITyr2=2014+1;

%%load SLP model data

%%CM4 Hist Ensemble Member 1
model='CM4';
run='CM4_historical';

dataDir=[work_dir,'/raw_data/',model,'/',run];
var=['slp'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
SLP_E1=data;
w_nh_1deg=w;
x_nh_1deg=x;
y_nh_1deg=y;
ifXY_nh_1deg=ifXY;

%%CM4 Hist Ensemble Member 2
model='CM4';
run='CM4_historical_noBling_Run895Rst290';

dataDir=[work_dir,'/raw_data/',model,'/',run];
var=['slp'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
SLP_E2=data;

%%CM4 Hist Ensemble Member 3
model='CM4';
run='CM4_historical_noBling_Run895Rst332';

dataDir=[work_dir,'/raw_data/',model,'/',run];
var=['slp'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
SLP_E3=data;

%%CM4 AMIP 
model='CM4';
run='CM4_amip';

dataDir=[work_dir,'/raw_data/',model,'/',run];
var=['slp'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
SLP_amip=data;

%compute clim over 1979-2014
ind1=(SLPyr1-1950)*12+1;
ind2=(SLPyr2-1950)*12;

[SLP_E1_clim SLP_E1_anom]=computeclim(SLP_E1(:,ind1:ind2));
[SLP_E2_clim SLP_E2_anom]=computeclim(SLP_E2(:,ind1:ind2));
[SLP_E3_clim SLP_E3_anom]=computeclim(SLP_E3(:,ind1:ind2));
[SLP_amip_clim SLP_amip_anom]=computeclim(SLP_amip(:,ind1:ind2));

SLP_EnsMean_clim=(SLP_E1_clim+SLP_E2_clim+SLP_E3_clim)/3;

%%%ERA-interim 1979-2014
 
model='era_interim';
run='slp';
dataDir=['/work/mib/raw_data/',model,'/',run];
var='slp_regrid_CM4';
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
SLP_era=data/100;%convert to hPa
w_era=w;
x_era=x;
y_era=y;
ifXY_era=ifXY;

%compute clim over 1979-2014
ind1=(SLPyr1-1979)*12+1;
ind2=(SLPyr2-1979)*12;
[SLP_era_clim SLP_era_anom]=computeclim(SLP_era(:,ind1:ind2));

%load sea ice thickness

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

var=['sithick'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
SIT_nh=data;

[SIC_clim_nh_E1 SIC_std_nh_E1 SIT_clim_nh_E1 SIT_std_nh_E1 extent_nh_E1 extent_clim_nh_E1 extent_anom_nh_E1 volume_nh_E1 volume_clim_nh_E1 volume_anom_nh_E1]=computeSeaIceDiagnostics(SIC_nh,SIT_nh,w_nh);

time_hist=1950+1/24:1/12:2015-1/24;
hist_years=1950:1:2014;
nT_hist=length(time_hist);

%define common period of 2010-2015 for computing climatology
ind1=(CM4SITyr1-1950)*12+1;
ind2=(CM4SITyr2-1950)*12;
[SIC_clim_common_nh_E1 SIC_std_common_nh_E1 SIT_clim_common_nh_E1 SIT_std_common_nh_E1 extent_common_nh_E1 extent_clim_common_nh_E1 extent_anom_common_nh_E1 volume_common_nh_E1 volume_clim_common_nh_E1 volume_anom_common_nh_E1]=computeSeaIceDiagnostics(SIC_nh(:,ind1:ind2),SIT_nh(:,ind1:ind2),w_nh);

%%CM4 Hist Ensemble Member 2
model='CM4';
run='CM4_historical_noBling_Run895Rst290';

dataDir=[work_dir,'/raw_data/',model,'/',run];

var=['siconc'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
SIC_nh=data;

var=['sithick'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
SIT_nh=data;

[SIC_clim_nh_E2 SIC_std_nh_E2 SIT_clim_nh_E2 SIT_std_nh_E2 extent_nh_E2 extent_clim_nh_E2 extent_anom_nh_E2 volume_nh_E2 volume_clim_nh_E2 volume_anom_nh_E2]=computeSeaIceDiagnostics(SIC_nh,SIT_nh,w_nh);

%define common period of 2010-2015 for computing climatology
ind1=(CM4SITyr1-1950)*12+1;
ind2=(CM4SITyr2-1950)*12;
[SIC_clim_common_nh_E2 SIC_std_common_nh_E2 SIT_clim_common_nh_E2 SIT_std_common_nh_E2 extent_common_nh_E2 extent_clim_common_nh_E2 extent_anom_common_nh_E2 volume_common_nh_E2 volume_clim_common_nh_E2 volume_anom_common_nh_E2]=computeSeaIceDiagnostics(SIC_nh(:,ind1:ind2),SIT_nh(:,ind1:ind2),w_nh);

%%CM4 Hist Ensemble Member 3
model='CM4';
run='CM4_historical_noBling_Run895Rst332';

dataDir=[work_dir,'/raw_data/',model,'/',run];

var=['siconc'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
SIC_nh=data;

var=['sithick'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
SIT_nh=data;

[SIC_clim_nh_E3 SIC_std_nh_E3 SIT_clim_nh_E3 SIT_std_nh_E3 extent_nh_E3 extent_clim_nh_E3 extent_anom_nh_E3 volume_nh_E3 volume_clim_nh_E3 volume_anom_nh_E3]=computeSeaIceDiagnostics(SIC_nh,SIT_nh,w_nh);

time_hist=1950+1/24:1/12:2015-1/24;
hist_years=1950:1:2014;
nT_hist=length(time_hist);

%define common period of 2010-2015 for computing climatology
ind1=(CM4SITyr1-1950)*12+1;
ind2=(CM4SITyr2-1950)*12;

[SIC_clim_common_nh_E3 SIC_std_common_nh_E3 SIT_clim_common_nh_E3 SIT_std_common_nh_E3 extent_common_nh_E3 extent_clim_common_nh_E3 extent_anom_common_nh_E3 volume_common_nh_E3 volume_clim_common_nh_E3 volume_anom_common_nh_E3]=computeSeaIceDiagnostics(SIC_nh(:,ind1:ind2),SIT_nh(:,ind1:ind2),w_nh);

SIT_clim_common_nh_EnsMean=(SIT_clim_common_nh_E1+SIT_clim_common_nh_E2+SIT_clim_common_nh_E3)/3;

%%%LOAD OM4 CORE2 data

%p25
model='OM4p25';
run='OM4p25_IAF_BLING_CFC_csf';

dataDir=[work_dir,'/raw_data/',model,'/',run];

var=['siconc_nh'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
SIC_nh=data;

var=['sithick_nh'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
SIT_nh=data;

[SIC_clim_nh_p25 SIC_std_nh_p25 SIT_clim_nh_p25 SIT_std_nh_p25 extent_nh_p25 extent_clim_nh_p25 extent_anom_nh_p25 volume_nh_p25 volume_clim_nh_p25 volume_anom_nh_p25]=computeSeaIceDiagnostics(SIC_nh,SIT_nh,w_nh);

time_CORE2=1948+1/24:1/12:2008-1/24;
CORE2_years=1948:1:2007;
nT_CORE2=length(time_CORE2);

%define common period of 1979-2007 for computing climatology
ind1=(SITyr1-1948)*12+1;
ind2=(SITyr2-1948)*12;
[SIC_clim_common_nh_p25 SIC_std_common_nh_p25 SIT_clim_common_nh_p25 SIT_std_common_nh_p25 extent_common_nh_p25 extent_clim_common_nh_p25 extent_anom_common_nh_p25 volume_common_nh_p25 volume_clim_common_nh_p25 volume_anom_common_nh_p25]=computeSeaIceDiagnostics(SIC_nh(:,ind1:ind2),SIT_nh(:,ind1:ind2),w_nh);

%AWI SIT (spans 2011-2018)

dataDir=['/work/mib/raw_data/CRYOSAT_AWI'];
filename=[dataDir,'/sit.mat'];
load(filename,'SIT','SIT_clim','SIC','volume','volume_clim','time','x','y','ifXY','w','area_weights')

AWI_SIT=SIT;
AWI_SIT_clim=SIT_clim;
AWI_volume=volume;
AWI_volume_clim=volume_clim;
AWI_time=time;
AWI_x=x;
AWI_y=y;
AWI_ifXY=ifXY;

%%add velocity vectors

K=3;
head_length=2.5;
shaft_width=0.8;
scale=50;

yrLim=[2010 2018];
dataDir=['/work/mib/raw_data/OSISAF/drift_lr_',num2str(yrLim(1)),'-',num2str(yrLim(2))];

filename=[dataDir,'/siu_trueE.mat'];
load(filename,'data','time','x','y','ifXY','w','area_weights')
x_OS=x;
y_OS=y;
ifXY_OS=ifXY;
siu=100*data;

filename=[dataDir,'/siv_trueN.mat'];
load(filename,'data','time','x','y','ifXY','w','area_weights')
siv=100*data;

ind1=(2011-2010)*12+1;
ind2=(2019-2010)*12;

[siu_clim siu_anom]=computeclim(siu(:,ind1:ind2));
[siv_clim siv_anom]=computeclim(siv(:,ind1:ind2));

siu_clim_OS=siu_clim;
siv_clim_OS=siv_clim;

%%CM4 Hist Ensemble Member 1
model='CM4';
run='CM4_historical';

dataDir=[work_dir,'/raw_data/',model,'/',run];
var=['siu_trueE'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
siu_E1=data*100;%cm/s
x_CM=x;
y_CM=y;
ifXY_CM=ifXY;

var=['siv_trueN'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
siv_E1=data*100;%cm/s

var=['siconc'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
sic=data;

%%CM4 Hist Ensemble Member 2
model='CM4';
run='CM4_historical_noBling_Run895Rst290';

dataDir=[work_dir,'/raw_data/',model,'/',run];
var=['siu_trueE'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
siu_E2=data*100;

%var=['siv'];
var=['siv_trueN'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
siv_E2=data*100;

%%CM4 Hist Ensemble Member 3
model='CM4';
run='CM4_historical_noBling_Run895Rst332';

dataDir=[work_dir,'/raw_data/',model,'/',run];
var=['siu_trueE'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
siu_E3=data*100;

%var=['siv'];
var=['siv_trueN'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
siv_E3=data*100;

%compute clim over 2005-2014
ind1=(2005-1950)*12+1;
ind2=(2015-1950)*12;

[siu_E1_clim siu_E1_anom]=computeclim(siu_E1(:,ind1:ind2));
[siv_E1_clim siv_E1_anom]=computeclim(siv_E1(:,ind1:ind2));
[siu_E2_clim siu_E2_anom]=computeclim(siu_E2(:,ind1:ind2));
[siv_E2_clim siv_E2_anom]=computeclim(siv_E2(:,ind1:ind2));
[siu_E3_clim siu_E3_anom]=computeclim(siu_E3(:,ind1:ind2));
[siv_E3_clim siv_E3_anom]=computeclim(siv_E3(:,ind1:ind2));
[sic_clim_CM sic_anom]=computeclim(sic(:,ind1:ind2));

siu_clim_CM=1/3*(siu_E1_clim+siu_E2_clim+siu_E3_clim);
siv_clim_CM=1/3*(siv_E1_clim+siv_E2_clim+siv_E3_clim);

%Load CORE IAF run

model='OM4p25';
run='OM4p25_IAF_BLING_CFC_csf';

dataDir=[work_dir,'/raw_data/',model,'/',run];
var=['siu_trueE'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
siu=data*100;

var=['siv_trueN'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
siv=data*100;

var=['siconc'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w','area_weights')
sic=data;

%compute clim over 1998-2007
ind1=(1998-1948)*12+1;
ind2=(2008-1948)*12;

[siu_clim siu_anom]=computeclim(siu(:,ind1:ind2));
[siv_clim siv_anom]=computeclim(siv(:,ind1:ind2));
[sic_clim_OM sic_anom]=computeclim(sic(:,ind1:ind2));

siu_clim_OM=siu_clim;
siv_clim_OM=siv_clim;
    
    nTileX=3;
    nTileY = 2;
        
    figWidth = 13.5; %6.5 inches = 39 picas, for a two column figure
    deltaX   = 0.4;
    deltaX2  = 0.7;
    deltaY   = 0.5;
    deltaY2  = 0.3;
    gapX     = 0.02;
    gapY     = 0.5;
    panel_scale=1;

ax=make_axes(nTileX,nTileY,figWidth,deltaX,deltaX2,deltaY,deltaY2,gapX,gapY,panel_scale);

SITmax=2.5;
radius=35;

SIT_vel_threshold=0.75;

axes(ax(1,1))
tmp_obs=AWI_SIT_clim(:,5);%March SIT
make_polar_plot_pcolor(AWI_x,AWI_y,AWI_ifXY,tmp_obs,radius,[''],0,0,SITmax)
colormap(gca,cmapSIT)
text(-0.3,0.67,'CryoSat-2 Obs','FontWeight','bold','FontSize',18)

SIT_mask=(tmp_obs>=SIT_vel_threshold);

tmp1=mean(siu_clim_OS(:,[1 2 3 12]),2);
tmp2=mean(siv_clim_OS(:,[1 2 3 12]),2);

make_polar_plot_vec(x_OS,y_OS,ifXY_OS,tmp1,tmp2,radius,'',K,head_length,shaft_width,scale)

text(-0.7,0.67,'A','FontWeight','bold','FontSize',18)
text(-0.7,-0.2,'March SIT','FontWeight','bold','FontSize',18,'rotation',90)

axes(ax(2,1))
tmp=SIT_clim_common_nh_EnsMean(:,3);
make_polar_plot_pcolor(x_nh,y_nh,ifXY_nh,tmp,radius,'',0,0,SITmax);

SIT_mask=(tmp>=SIT_vel_threshold);

tmp1=mean(siu_clim_CM(:,[1 2 3 12]),2);
tmp2=mean(siv_clim_CM(:,[1 2 3 12]),2);
SIC=mean(sic_clim_CM(:,[1 2 3 12]),2);
SIC_mask=SIC>=0.15;
tmp1=tmp1.*SIC_mask.*SIT_mask;
tmp2=tmp2.*SIC_mask.*SIT_mask;

make_polar_plot_vec(x_CM,y_CM,ifXY_CM,tmp1,tmp2,radius,'',K,head_length,shaft_width,scale)

text(-0.7,0.67,'B','FontWeight','bold','FontSize',18)
text(-0.1,0.67,'CM4.0','FontWeight','bold','FontSize',18)
colormap(gca,cmapSIT)

axes(ax(3,1))
tmp=SIT_clim_common_nh_p25(:,3);
make_polar_plot_pcolor(x_nh,y_nh,ifXY_nh,tmp,radius,'',0,0,SITmax);
colormap(gca,cmapSIT)

SIT_mask=(tmp>=SIT_vel_threshold);

tmp1=mean(siu_clim_OM(:,[1 2 3 12]),2);
tmp2=mean(siv_clim_OM(:,[1 2 3 12]),2);
SIC=mean(sic_clim_OM(:,[1 2 3 12]),2);
SIC_mask=SIC>=0.15;
tmp1=tmp1.*SIC_mask.*SIT_mask;
tmp2=tmp2.*SIC_mask.*SIT_mask;

make_polar_plot_vec(x_CM,y_CM,ifXY_CM,tmp1,tmp2,radius,'',K,head_length,shaft_width,scale)


text(-0.7,0.67,'C','FontWeight','bold','FontSize',18)
text(-0.3,0.67,'OM4.0 CORE IAF','FontWeight','bold','FontSize',18)

                 pos=get(gca,'position');

                 c=colorbar;

                 set(gca,'position',pos)

                 x_pos=get(c, 'Position');
                 y_pos=[x_pos(1)-0.01 x_pos(2) 0.7*x_pos(3) 1*x_pos(4)];
                 set(c,'Position',y_pos)


axes(ax(1,2))
tmp2=mean(SLP_era_clim(:,[1 2 3 12]),2);
hold on
make_polar_plot_pcolor_land(x_era,y_era,ifXY_era,tmp2,radius,'',0,995,1030);
text(-0.7,0.67,'D','FontWeight','bold','FontSize',18)
text(-0.1,0.67,'ERAi','FontWeight','bold','FontSize',18)
text(-0.7,-0.1,'DJFM SLP','FontWeight','bold','FontSize',18,'rotation',90)
colormap(gca,cmapSLP)

axes(ax(2,2))
tmp1=mean(SLP_EnsMean_clim(:,[1 2 3 12]),2);
hold on
make_polar_plot_pcolor_land(x_nh_1deg,y_nh_1deg,ifXY_nh_1deg,tmp1,radius,'',0,995,1030);
text(-0.7,0.67,'E','FontWeight','bold','FontSize',18)
text(-0.1,0.67,'CM4.0','FontWeight','bold','FontSize',18)
colormap(gca,cmapSLP)

axes(ax(3,2))
tmp1=mean(SLP_amip_clim(:,[1 2 3 12]),2);
make_polar_plot_pcolor_land(x_nh_1deg,y_nh_1deg,ifXY_nh_1deg,tmp1,radius,'',0,995,1030);
text(-0.7,0.67,'F','FontWeight','bold','FontSize',18)
text(-0.2,0.67,'AM4.0 AMIP','FontWeight','bold','FontSize',18)
colormap(gca,cmapSLP)

                 pos=get(gca,'position');

                 c=colorbar;

                 set(gca,'position',pos)

                 x_pos=get(c, 'Position');
                 y_pos=[x_pos(1)-0.01 x_pos(2) 0.7*x_pos(3) 1*x_pos(4)];
                 set(c,'Position',y_pos)


set(gcf,'Color','w')
filename=[fig_dir,'/SLP_DJFM_SIT_M2_vel_',num2str(SIT_vel_threshold)];
export_fig(gcf,filename,'-a1','-png','-r864','-painters');


