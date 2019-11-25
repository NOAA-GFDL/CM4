function script_RossPolynya(fig_dir,work_dir)

close all

addpath('./utils/export_fig/')
addpath('./utils/')
addpath('./utils/m_map/')

month_str={'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'};

doload=1;

if(doload==1)
model='CM4';
run='CM4_piControl';
var='MLD_003_sh';

dataDir=[work_dir,'/raw_data/',model,'/',run];
filename=[dataDir,'/',var,'.mat']

load(filename,'data','time','x','y','ifXY','w','area_weights')
MLD=data;

var='siconc_sh';

dataDir=[work_dir,'/raw_data/',model,'/',run];
filename=[dataDir,'/',var,'.mat']

load(filename,'data','time','x','y','ifXY','w','area_weights')

SIC=data/100;
HI=zeros(size(SIC));

disp('Computing diagnostics...')
[SIC_clim SIC_std HI_clim HI_std extent extent_clim extent_anom volume volume_clim volume_anom]=computeSeaIceDiagnostics(SIC,HI,w);
disp('Done!')

%%define regional sector for SIE computation
xmin=160;
xmax=300;

ifXY_Ross=(x>=xmin & ...
     x<=xmax) & ...
     ifXY;

A=ones(size(ifXY)).*double(ifXY_Ross);
region_mask=A(ifXY);

nD=size(SIC,1);
nT=size(SIC,2);
SIC_Ross=zeros(nD,nT);
for i=1:nT
SIC_Ross(:,i)=SIC(:,i).*region_mask;
end

HI_Ross=zeros(size(SIC_Ross));

[SIC_clim_Ross SIC_std_Ross HI_clim_Ross HI_std_Ross extent_Ross extent_clim_Ross extent_anom_Ross volume_Ross volume_clim_Ross volume_anom_Ross]=computeSeaIceDiagnostics(SIC_Ross,HI_Ross,w);

xmin=160;
xmax=300;
ymin=-90;
ymax=-60;

ifXY_Ross=(x>=xmin & ...
     x<=xmax) & ...
     y<=ymax & ...
     ifXY;

A=ones(size(ifXY)).*double(ifXY_Ross);
region_mask=A(ifXY);

nT=size(SIC,2);
MLDmax_Ross=zeros(1,nT);
MLDmean_Ross=zeros(1,nT);
for i=1:nT
MLDmax_Ross(:,i)=max(MLD(:,i).*region_mask);
end

%%save 3D ocean temperature data

parse3DOcean(work_dir)
TempHovmoller(work_dir)

end

month=9;
time=1+1/24:1/12:501-1/24;
time_Sep=time(month:12:end);
extent_anom_Ross_Sep=extent_anom_Ross(month:12:end);
extent_anom_Sep=extent_anom(month:12:end);
extent_Ross_Sep=extent_Ross(month:12:end);
extent_Sep=extent(month:12:end);

MLDmax_Ross_Sep=MLDmax_Ross(month:12:end);

nTileX = 5;
nTileY = 4;

    figWidth = 7; %6.5 inches = 39 picas, for a two column figure

    deltaX   = 0.55;%side spacing?
    deltaX2  = 0.65;%side spacing?
    deltaY   = 0.5;
    deltaY2  = 0.4;
    cbar = 0;
    gapX     = 0.02;
    gapY     = 0.85;
    panelX   = ( figWidth - deltaX - deltaX2 - ( nTileX -1 ) * gapX ) / nTileX;
    panelY   = 1*panelX ;%aspect ratio

scale=5;

    posn     = [ 0, ...
                 0, ...
                 nTileX * panelX + ( nTileX - 1 ) * gapX + deltaX + deltaX2, ...
                 nTileY * panelY + ( nTileY - 1 ) * gapY + deltaY + deltaY2 ];

    fig = figure( 'units', 'inches', ...
                  'paperPosition', posn, ...
                  'position', posn, ...
                  'defaultAxesNextPlot', 'add', ...
                  'defaultAxesBox', 'on', ...
                  'defaultAxesFontSize', 12, ...
                  'defaultTextFontSize', 12, ...
                  'defaultAxesTickDir', 'out', ...
                  'defaultAxesTickLength', [ 0.01 0 ], ...
                  'defaultAxesFontName', 'times', ...
                  'defaultTextFontName', 'times', ...
                  'defaultAxesLineWidth', 2, ...
                  'defaultAxesLayer', 'top' );

  ax = zeros( nTileX, nTileY );

    for iAx = 1 : nTileX
        for jAx = 1 : 1
            ax( iAx, jAx ) = axes( 'units', 'inches', ...
                                   'position', [ deltaX + ( iAx - 1 ) * ( panelX + gapX ), ...
                                    deltaY + ( nTileY - jAx ) * ( panelY + gapY ), ...
                                    panelX, panelY ] );
        end
    end
    for iAx = 1 : nTileX
        for jAx = 2 : 2
            ax( iAx, jAx ) = axes( 'units', 'inches', ...
                                   'position', [ deltaX + ( iAx - 1 ) * ( panelX + gapX ), ...
                                    deltaY + ( nTileY - jAx ) * ( panelY + gapY )+0.7*gapY, ...
                                    panelX, panelY ] );
        end
    end
    
   for iAx = 1 : 1
        for jAx = nTileY-1 : nTileY-1
            ax( iAx, jAx ) = axes( 'units', 'inches', ...
                                   'position', [ deltaX + ( iAx - 1 ) * ( panelX + gapX ), ...
                                    deltaY + ( nTileY - jAx ) * ( panelY + gapY ) + 0.6*gapY, ...
                                    scale*panelX, 1.5*panelY ] );
        end
    end
    
    for iAx = 1 : 1
        for jAx = nTileY : nTileY
            ax( iAx, jAx ) = axes( 'units', 'inches', ...
                                   'position', [ deltaX + ( iAx - 1 ) * ( panelX + gapX ), ...
                                    deltaY + ( nTileY - jAx ) * ( panelY + gapY ), ...
                                    scale*panelX, 1.5*panelY ] );
        end
    end

years=[442:1:451];%Note: 500 year PI control is years 151-650 of this run.
years=years-150;

labels={'A','B','C','D','E','F','G','H','I','J'};

for row=1:2
for col=1:5

k=col+(row-1)*nTileX;

i=(years(k)-1)*12+9;%Sept ind

axes(ax(col,row))

radius=40;
iscbar=0;
cmin=-0.5;
cmax=0.5;

month=mod(i,12);
if(month==0)
month=12;
end

data=SIC(:,i)-SIC_clim(:,month);

title_str=['SIC Anom ',month_str{month},' Year ',num2str(1+floor((i-1)/12))];

make_polar_plot_pcolor_sh(x,y,ifXY,data,radius,'',iscbar,cmin,cmax)

text(-0.6,0.85,[labels{k},': ','Year ',num2str(1+floor((i-1)/12))],'FontWeight','bold')

if(col==1 && row==2)
text(-1.1,-0.2,['September SIC anomalies'],'Rotation',90)
end
if(col==5 && row==1)
                 pos=get(gca,'position');

                 c=colorbar('vertical');

                 set(gca,'position',pos)

                 x_pos=get(c, 'Position');
                 y_pos=[x_pos(1)+0.06 x_pos(2)-0.15 1.5*x_pos(3) 3*x_pos(4)];
                 set(c,'Position',y_pos)
end

end
end

fact=0.0005;
shift=2.0;

axes(ax(1,nTileY-1))
hold on
plot(time_Sep,1e-12*extent_Ross_Sep,'k')
plot(time_Sep,fact*MLDmax_Ross_Sep+shift,'r')
title('Ross, Amundsen and Bellingshausen September Sea Ice Extent')
xlabel('Time (years)')
ylabel('Sea Ice Extent (M km^2)')
xlim([0 501])
ylim([shift 8.5])
line([years(1)-0.5 years(1)-0.5],[shift-0.05,8.5],'Color','b','LineStyle','--')
line([years(end)+0.5 years(end)+0.5],[shift-0.05,8.5],'Color','b','LineStyle','--')
set(gca,'yTick',2:0.5:8,'yTickLabel',{'2','','3','','4','','5','','6','','7','','8'})
text(-30,9,'K','FontWeight','bold')

for i=2:0.5:5
text(507,i-0.03,[num2str((i-shift)/(fact*1000))],'Color','r');
end
text(525,8.5,'Maximum MLD (1000 m)','Rotation',270,'Color','r')


model='CM4';
run='CM4_piControl';
dataPath=[work_dir,'/raw_data/',model,'/',run,'/']; 
load([dataPath,'TempHov_Ross.mat'],'TempHov','T','Z','z')

time=1:1:500;
[T,Z]=meshgrid(time,z);
nZ=size(TempHov,1);
nT=length(time);

TempHovWinter=zeros(nZ,nT);

for i=1:nT

ind1=7+(i-1)*12;
ind2=9+(i-1)*12;

TempHovWinter(:,i)=mean(TempHov(:,ind1:ind2),2);

end 

Z=Z/1000;

axes(ax(1,nTileY))
pcolor(T,Z,TempHovWinter)
shading interp

                 pos=get(gca,'position');

                 c=colorbar('vertical');

                 set(gca,'position',pos)

                 x_pos=get(c, 'Position');
                 y_pos=[x_pos(1)+0.09 x_pos(2) 0.5*x_pos(3) 1*x_pos(4)];
                 set(c,'Position',y_pos)
caxis([-1.5 1.5])

set(gca,'yDir','rev')
ylim([0 5])
title('Ross Subsurface Temperatures (\circ C)')
xlabel('Time (years)')
ylabel('Depth (1000m)')
xlim([0 501])
text(-30,-0.1,'L','FontWeight','bold')

line([years(1)-0.5 years(1)-0.5],[0 5],'Color','b','LineStyle','--')
line([years(end)+0.5 years(end)+0.5],[0 5],'Color','b','LineStyle','--')

set(gcf,'Color','w')
filename=[fig_dir,'/polynya2'];
export_fig(gcf,filename,'-a1','-png','-r600','-nocrop','-painters');

