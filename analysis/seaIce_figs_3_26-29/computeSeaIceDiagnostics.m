%this function takes raw data from a control run and computes sea ice diagnostics
%Note: If no thickness data, simply pass the function HI=zeros(size(SIC));

function [SIC_clim SIC_std HI_clim HI_std extent extent_clim extent_anom volume volume_clim volume_anom SIC_std_detrend HI_std_detrend]=computeSeaIceDiagnostics(SIC,HI,w);

nD=size(SIC,1);
nT=size(SIC,2);

SIC_clim=zeros(nD,12);
SIC_std=zeros(nD,12);
SIC_std_detrend=zeros(nD,12);
HI_clim=zeros(nD,12);
HI_std=zeros(nD,12);
HI_std_detrend=zeros(nD,12);

for month=1:12
SIC_clim(:,month)=mean(SIC(:,month:12:end),2);
SIC_std(:,month)=std(SIC(:,month:12:end),0,2);
SIC_std_detrend(:,month)=std(detrend(SIC(:,month:12:end)'),0,1);
HI_clim(:,month)=mean(HI(:,month:12:end),2);
HI_std(:,month)=std(HI(:,month:12:end),0,2);
HI_std_detrend(:,month)=std(detrend(HI(:,month:12:end)'),0,1);
end

extent=zeros(nT,1);
volume=zeros(nT,1);

for i=1:nT
tmp=SIC(:,i);
tmp=double(tmp>=0.15);
extent(i)=sum(tmp.*w);

tmpSIC=SIC(:,i);
tmpHI=HI(:,i);
volume(i)=sum(tmpSIC.*tmpHI.*w);
end

extent_clim=zeros(12,1);
volume_clim=zeros(12,1);
for month=1:12
extent_clim(month)=mean(extent(month:12:end));
volume_clim(month)=mean(extent(month:12:end));
end

extent_anom=zeros(nT,1);
volume_anom=zeros(nT,1);

extent_anom=extent-repmat(extent_clim,nT/12,1);
volume_anom=volume-repmat(volume_clim,nT/12,1);









