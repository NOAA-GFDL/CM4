%run this script to produce CM4 sea ice analysis (process raw data, compute diagnostics, make figures)

addpath('./utils/')

work_dir='/work/mib';
fig_dir='/home/Mitchell.Bushuk/Figures/CM4/paper_figs';


%%%%%%%process sea ice raw data

%OM4p25 COREII
model{1}='OM4p25';
run{1}='OM4p25_IAF_BLING_CFC_csf_rerun';
rootdir{1}='/archive/ogrp/CMIP6/OMIP/warsaw_201803_mom6_2018.04.06/OM4p25_IAF_BLING_csf_rerun/gfdl.ncrc3-intel16-prod/pp';
year_start{1}=1948;
year_end{1}=2007;
delta_year{1}=20;

%CM4 Historical Member 1
model{2}='CM4';
run{2}='CM4_historical';
rootdir{2}='/archive/oar.gfdl.cmip6/CM4/warsaw_201710_om4_v1.0.1/CM4_historical/gfdl.ncrc4-intel16-prod-openmp/pp';
year_start{2}=1950;
year_end{2}=2014;
delta_year{2}=5;

%CM4 Historical Member 2
model{3}='CM4';
run{3}='CM4_historical_noBling_Run895Rst290';
rootdir{3}='/archive/Michael.Winton/CM4/warsaw_201710_om4_v1.0.1/CM4_historical_noBling_Run895Rst290/gfdl.ncrc4-intel16-prod-openmp/pp';
year_start{3}=1950;
year_end{3}=2014;
delta_year{3}=5;

%CM4 Historical Member 3
model{4}='CM4';
run{4}='CM4_historical_noBling_Run895Rst332';
rootdir{4}='/archive/Michael.Winton/CM4/warsaw_201710_om4_v1.0.1/CM4_historical_noBling_Run895Rst332/gfdl.ncrc4-intel16-prod-openmp/pp';
year_start{4}=1950;
year_end{4}=2014;
delta_year{4}=5;

%CM4 PI control
model{5}='CM4';
run{5}='CM4_piControl_C';
rootdir{5}='/archive/oar.gfdl.cmip6/CM4/warsaw_201710_om4_v1.0.1/CM4_piControl_C/gfdl.ncrc4-intel16-prod-openmp/pp';
year_start{5}=151;
year_end{5}=650;
delta_year{5}=5;

comp='ice';
var_list={'siconc','sithick','siu','siv'};

for j=1:length(model)
for k=1:length(var_list)

gfdlRawData_loop(model{j},run{j},comp,var_list{k},rootdir{j},year_start{j},year_end{j},delta_year{j},work_dir)

end 
end

%process ocean data for PI control run

comp='ocean_monthly_z_1x1deg';
var_list={'thetao'};

for j=5:5
for k=1:length(var_list)

gfdlRawData_loop(model{j},run{j},comp,var_list{k},rootdir{j},year_start{j},year_end{j},delta_year{j},work_dir)

end 
end

comp='ocean_monthly';
var_list={'MLD_003'};

for j=5:5
for k=1:length(var_list)

gfdlRawData_loop(model{j},run{j},comp,var_list{k},rootdir{j},year_start{j},year_end{j},delta_year{j},work_dir)

end 
end


%%%%%%%%%compute diagnostics and make figures

%Fig 3

script_RossPolynya(fig_dir,work_dir)

%Figs 26, 27, 28

CM4_paper_figs(fig_dir,work_dir)

%Figs 29

script_SLP_SIT(fig_dir,work_dir)




















