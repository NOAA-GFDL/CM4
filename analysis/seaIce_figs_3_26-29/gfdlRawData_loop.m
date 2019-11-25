function gfdlRawData_loop(model,run,comp,var,rootdir,year_start,year_end,delta_year,work_dir) 

if(strcmp(var,'CN'))
var_name='sic';
else
var_name=var;
end

close all

grid=fullfile(rootdir,comp,[comp,'.static.nc'])

if(strcmp(comp,'ice'))

x=ncget(grid,'GEOLON');
y=ncget(grid,'GEOLAT');
area_weights=ncget(grid,'CELL_AREA');

ifOcean=(area_weights>0);

%area is proportion of the earths surface. Convert to m^2
radius=6.371e6;
area_weights=area_weights*4*pi*radius^2; 
end

if(strcmp(comp,'ocean_monthly') )

x=ncget(grid,'geolon');
y=ncget(grid,'geolat');
area_weights=ncget(grid,'areacello');

wet=ncget(grid,'wet');
ifOcean=(wet>0);

end

%x runs from -280 to 80. Convert to 0 to 360
x=mod(x+360,360);

%create NH mask
xmin=0;
xmax=360;
ymin=40;
ymax=90;

ifXY_nh=x>=xmin & ...
     x<=xmax & ...
     y>=ymin & ...
     y<=ymax & ...
     ifOcean;

%%create SH mask

%create SH mask
xmin=0;
xmax=360;
ymin=-90;
ymax=-45;

ifXY_sh=x>=xmin & ...
     x<=xmax & ...
     y>=ymin & ...
     y<=ymax & ...
     ifOcean;

nD=nnz(ifXY_nh)

w=area_weights(ifXY_nh);

year_list=year_start:delta_year:year_end;

nT=(year_end-year_start+1)*12;
nMonth=12*delta_year;

data=zeros(nD,nT);

for k=1:length(year_list)

data_path=fullfile(rootdir,comp,'ts/monthly',[num2str(delta_year),'yr'],[comp,'.',num2str(year_list(k),'%04i'),'01-',num2str(year_list(k)+delta_year-1,'%04i'),'12.',var,'.nc'])

disp('Loading Data...')
ts=ncget(data_path,var);
disp('Done!')
time_seg=ncget(data_path,'time');

size(ts)

if(strcmp(var,'CN'))
sic=squeeze(sum(ts,3));
ts=sic;
end

shift=nMonth*(k-1);

for i=1:nMonth
dataTmp=ts(:,:,i);
data(:,shift+i)=dataTmp(ifXY_nh);
end

time(shift+1:shift+nMonth)=time_seg;

end

dataDir=[work_dir,'/raw_data/',model,'/',run];

if(~isdir(dataDir))
mkdir(dataDir)
end

filename=[dataDir,'/',var_name,'_nh.mat']

ifXY=ifXY_nh;
w=area_weights(ifXY_nh);

save(filename,'data','time','x','y','ifXY','w','area_weights','year_start','year_end','-v7.3')


nD=nnz(ifXY_sh)

data=zeros(nD,nT);

for k=1:length(year_list)

data_path=fullfile(rootdir,comp,'ts/monthly',[num2str(delta_year),'yr'],[comp,'.',num2str(year_list(k),'%04i'),'01-',num2str(year_list(k)+delta_year-1,'%04i'),'12.',var,'.nc'])

disp('Loading Data...')
ts=ncget(data_path,var);
disp('Done!')
time_seg=ncget(data_path,'time');

size(ts)

if(strcmp(var,'CN'))
sic=squeeze(sum(ts,3));
ts=sic;
end

shift=nMonth*(k-1);

for i=1:nMonth
dataTmp=ts(:,:,i);
data(:,shift+i)=dataTmp(ifXY_sh);
end

time(shift+1:shift+nMonth)=time_seg;

end

dataDir=[work_dir,'/raw_data/',model,'/',run];

if(~isdir(dataDir))
mkdir(dataDir)
end

filename=[dataDir,'/',var_name,'_sh.mat']

ifXY=ifXY_sh;
w=area_weights(ifXY_sh);

save(filename,'data','time','x','y','ifXY','w','area_weights','year_start','year_end','-v7.3')



