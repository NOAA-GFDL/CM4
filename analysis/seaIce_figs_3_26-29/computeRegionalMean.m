%this function takes a matrix of data and a region mask, and returns regional mean values

function [regionalMean regionalMeanclim regionalMeansig]=computeRegionalMean(data,w,ifXY,ifXY_reg);

nD=size(data,1);
nT=size(data,2);

nReg=size(ifXY_reg,3);

regionalMean=zeros(nT,nReg);
regionalMeanclim=zeros(12,nReg);
regionalMeansig=zeros(12,nReg);

for region=1:nReg

A=ones(size(ifXY)).*double(ifXY_reg(:,:,region));

region_mask=A(ifXY);

data_tmp = bsxfun(@times, data, region_mask);

for i=1:nT
regTmp=data_tmp(:,i);
goodinds=(~isnan(regTmp));
area=sum(region_mask.*w.*goodinds);
regionalMean(i,region)=nansum(regTmp.*w)/area;
end

for month=1:12
regionalMeanclim(month,region)=mean(regionalMean(month:12:end,region));
regionalMeansig(month,region)=std(regionalMean(month:12:end,region));
end

end








