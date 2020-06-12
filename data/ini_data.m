function [ Xcal,Xtest,ycal,ytest,wv ] = ini_data( data_url,chemical_id,ratio )
%INI_DATA Summary of this function goes here
%   Detailed explanation goes here
if nargin<3;ratio=0.6;end
if nargin<2;chemical_id=1;end

if strcmp(data_url,'data1_')
    load('data1_milk.mat');
    [Xcal,ycal,Xtest,ytest]=sampling_sub(X,protein,ratio);
    wv=linspace(4000,10000,1557);
elseif strcmp(data_url,'data2_')
    load('data2_nir_shootout_2002.mat');
    cal_outlier=false(size(calibrate_1.data,1),1);test_outlier=false(size(test_1.data,1),1);
    cal_outlier([19,122,126,127])=true;      test_outlier([11,145,267,295,294,342,313,341,343])=true;
    Xcal=calibrate_1.data(~cal_outlier,:); ycal=calibrate_Y.data(~cal_outlier,chemical_id);
    Xtest=test_1.data(~test_outlier,:);     ytest=test_Y.data(~test_outlier,chemical_id);
    wv=600:2:1898;
elseif strcmp(data_url,'data3_')
    load('data3_corn.mat');
    [Xcal,ycal,Xtest,ytest]=sampling_sub(mp5spec.data,propvals.data(:,chemical_id),ratio);
    wv=1100:2:2498;
elseif strcmp(data_url,'data4_')
    load('data4_gasoil60.mat');
    [Xcal,ycal,Xtest,ytest]=sampling_sub(NIR,octane,ratio);
    wv=900:2:1700;
elseif strcmp(data_url,'data5_')
    load('data5_shootout2008_wheat884.mat');
    ycal=ycal(:,chemical_id);ytest=ytest(:,chemical_id);
elseif strcmp(data_url,'data6_')
    load('data6_nirbeer.mat');
    wv=xaxis;
elseif strcmp(data_url,'data7_')
    load('data7_IDRCShootOut2010Transmit.mat');
    Xcal=XcalTrans;ycal=YcalTrans(:,chemical_id);
    Xtest=XvalTrans;ytest=YvalTrans(:,chemical_id);
    wv=1100:2:2498;  %the website does not provide the wavelength
elseif strcmp(data_url,'data8_')
    load('data8_IDRCShootOut2010Reflect.mat');
    Xcal=XcalReflect;ycal=YcalReflect(:,chemical_id);
    Xtest=XvalReflect;ytest=YvalReflect(:,chemical_id);
    wv=1100:2:2498; %the website does not provide the wavelength
elseif strcmp(data_url,'data9_')
    load('data9_NIRdata_tablets');
    [Xcal,ycal,Xtest,ytest]=sampling_sub(Matrix(:,4:end),Matrix(:,chemical_id),ratio);
    wv=linspace(7398.337,10507.30,404);
elseif strcmp(data_url,'data10_')
    load('data10_Ramandata_tablets.mat')
    [Xcal,ycal,Xtest,ytest]=sampling_sub(Matrix(:,3:end),Matrix(:,chemical_id),ratio);
    wv=linspace(3600.000,200.0000,3401);
elseif strcmp(data_url,'data11_')
    load('data11_marzipanWEBdata.mat')
    if chemical_id==1;[Xcal,ycal,Xtest,ytest]=sampling_sub(NIRS1,ref_moisture,ratio);
    else; [Xcal,ycal,Xtest,ytest]=sampling_sub(NIRS1,ref_sugar,ratio);end
    wv=NIRS1_axis;
elseif strcmp(data_url,'data12_')
    load('data12_NIR_sugar.mat')
    if chemical_id==1;[Xcal,ycal,Xtest,ytest]=sampling_sub(X.data,Brix.data(:,chemical_id),ratio);
    else;[Xcal,ycal,Xtest,ytest]=sampling_sub(X.data,pol.data,ratio);end
    wv=400:2:1888;
elseif strcmp(data_url,'data13_')
    load('data13_NIRsoil.mat')
    [Xcal,ycal,Xtest,ytest]=sampling_sub(soil.data,soilref.data(:,chemical_id),ratio);
    wv=linspace( 400,2500,1050);
elseif strcmp(data_url,'data14_')
    load('data14_NITSingleSeed.mat');
    Xcal=Calibration_X;ycal=Calibration_Y;
    Xtest=Validation_X;ytest=Validation_Y;
    wv=850:2:1048;
elseif strcmp(data_url,'data15_')
    load('data15_RAMANPorkFat.mat')
    [Xcal,ycal,Xtest,ytest]=sampling_sub(X.data,Y.data(:,chemical_id),ratio);
    wv=X.axisscale{2,1};
elseif contains(data_url,'data16')
    load('data16_shao_ynzy.mat')
    [row,col]=size(data.spHHPE);
    chain=ks(data.spHHPE);
    [Xcal,ycal,Xtest,ytest]=sampling_sub(data.spHHPE,data.value(:,chemical_id),ratio);
    wv=linspace(10000,4000,3001);
elseif strcmp(data_url,'data17_')
    load('data17_carra.mat')
    Xcal=Xcal;ycal=Ycal;
    Xtest=Xval;ytest=Yval;
    wv=linspace(1100,2500,699);
elseif strcmp(data_url,'data18_')
    load('data18_gasoil.mat')
    Xcal=Xcal;ycal=Ycal;
    Xtest=Xval;ytest=Yval;
    wv=linspace(4900,9000,2128);%cm-1
elseif strcmp(data_url,'data19_')
    load('data19_water.mat')
    Xcal=Xcal;ycal=Ycal;
    Xtest=Xval;ytest=Yval;
    wv=linspace(850,1048,2151);%nm
elseif strcmp(data_url,'data20_')
    load('data20_bp50gatest.mat')
    Xcal=[bp50_s1d_hl;bp50_s1d_ll_a];ycal=[bp50_y1_hl;bp50_y1_ll_a];
    Xtest=bp50_s1d_ll_b;ytest=bp50_y1_ll_b;
    wv=750:2:1550;%nm
elseif strcmp(data_url,'data21_')
    load('data21_cngatest.mat')
    Xcal=[cn_sd_hl;cn_sd_ll_a];ycal=[cn_y_hl;cn_y_ll_a];
    Xtest=cn_sd_ll_b;ytest=cn_y_ll_b;
    wv=750:2:1550;%nm
elseif strcmp(data_url,'data22_')
    load('data22_d4052gatest.mat')
    Xcal=[d_sd_hl;d_sd_ll_a];ycal=[d_y_hl;d_y_ll_a];
    Xtest=d_sd_ll_b;ytest=d_y_ll_b;
    wv=750:2:1550;%nm
elseif strcmp(data_url,'data23_')
    load('data23_freezegatest.mat')
    Xcal=[f_sd_hl;f_sd_ll_a];ycal=[f_y_hl;f_y_ll_a];
    Xtest=f_sd_ll_b;ytest=f_y_ll_b;
    wv=750:2:1550;%nm
elseif strcmp(data_url,'data24_')
    load('data24_totalgatest.mat')
    Xcal=[t_sd_hl;t_sd_ll_a];ycal=[t_y_hl;t_y_ll_a];
    Xtest=t_sd_ll_b;ytest=t_y_ll_b;
    wv=750:2:1550;%nm
elseif strcmp(data_url,'data25_')
    load('data25_viscgatest.mat')
    Xcal=[v_sd_hl;v_sd_ll_a];ycal=[v_y_hl;v_y_ll_a];
    Xtest=v_sd_ll_b;ytest=v_y_ll_b;
    wv=750:2:1550;%nm
elseif strcmp(data_url,'data26_')
    load('data26_juice.mat')
    Xcal=OJ_x_learning;
    ycal=OJ_y_learning;
    Xtest=OJ_x_test;
    ytest=OJ_y_test;
    wv=linspace(1100,2500,700);%nm
elseif strcmp(data_url,'data27_')  
    load('data27_serum_silver_wsy.mat')
    Xcal=data(cal_idx,vsel_mcuve);
    Xtest=data(test_idx,vsel_mcuve);
    ycal=prop(cal_idx,26);
    ytest=prop(test_idx,26);
    wv=wv(vsel_mcuve);
else
    error('input is not exist')
end
end

function [Xcal,ycal,Xtest,ytest]=sampling_sub(X,y,ratio)
choice='y_order';% ±¸Ñ¡ y_order;ks
switch choice
    case 'ks'
        rank=ks(X);
        cal=rank(1:fix(size(X,1)*ratio));
        test=setdiff(rank,cal);
    case 'y_order'
        [~,y_order]=sort(y);
        n_cal=fix(size(y,1)*ratio);
        cal=linspace(1,size(y,1),n_cal);
        cal=unique(y_order(round(cal)));
        test=setdiff(y_order,cal);
end
Xcal=X(cal,:);ycal=y(cal,:);
Xtest=X(test,:);ytest=y(test,:);
end
