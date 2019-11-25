function parse3DOcean(work_dir)

%this script produces 2D fields of interest using 3D ocean data

varname='thetao_sh';

model='CM4';
run='CM4_piControl';

dataPath=[work_dir,'/raw_data/',model,'/',run,'/']; 

disp('Loading data...')
load([dataPath,varname,'.mat'])
disp('Done!')

size(x)
size(y)
size(z)
size(data)

nZ=size(data,2);

for i=1:nZ
nnz(ifXY(:,:,i));
end

fulldata=data;

ifXY_ocean=ifXY(:,:,1);

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

data=squeeze(fulldata(:,i,:));

ifXY_ice=ifXY(:,:,i);
ifXY_ocean=ifXY(:,:,i);

save([dataPath,var,'.mat'],'data','time','x','y','z','ifXY_ice','ifXY_ocean','area_weights')

end




