function TempHovmoller(work_dir)

varname='thetao_sh';

model='CM4';
run='CM4_piControl';

dataPath=[work_dir,'/raw_data/',model,'/',run,'/']; 
load([dataPath,varname,'.mat'],'x','y','z','ifXY','area_weights')
nZ=length(z);
ifXY=squeeze(ifXY(:,:,1));

w=area_weights(ifXY);

xmin=160;
xmax=230;
ymin=-90;
ymax=-60;

ifXY_Ross=(x>=xmin & ...
     x<=xmax) & ...
     y<=ymax & ...
     ifXY;

depth=num2str(ceil(z(1)),'%04i')
var=[depth,'mTemp_sh']
load([dataPath,var,'.mat'],'data')

nT=size(data,2);

TempHov=zeros(nZ,nT);

for i=1:nZ

depth=num2str(ceil(z(i)),'%04i')

if(strcmp(varname,'thetao_nh'))
var=[depth,'mTemp_nh']
end
if(strcmp(varname,'thetao_sh'))
var=[depth,'mTemp_sh']
end
if(strcmp(varname,'so_nh'))
var=[depth,'mSalt_nh']
end
if(strcmp(varname,'so_sh'))
var=[depth,'mSalt_sh']
end

load([dataPath,var,'.mat'],'data')

[regionalMean regionalMeanclim regionalMeansig]=computeRegionalMean(data,w,ifXY,ifXY_Ross);

TempHov(i,:)=regionalMean;

end

time=1+1/24:1/12:501-1/24;

[T,Z]=meshgrid(time,z);

save([dataPath,'TempHov_Ross.mat'],'TempHov','T','Z','z')


