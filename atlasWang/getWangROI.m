clear; clc; 
whichhemis = 'rh'; % specify which hemiphere

[atl,XYZWang]=spm_read_vols(spm_vol(['maxprob_vol_',whichhemis,'.nii']));

regLab={
'V1v'
'V1d'
'V2v'
'V2d'
'V3v'
'V3d'
'hV4'
'VO1'
'VO2'
'PHC1'
'PHC2'
'MST'
'hMT'
'LO2'
'LO1'
'V3b'
'V3a'
'IPS0'
'IPS1'
'IPS2'
'IPS3'
'IPS4'
'IPS5'
'SPL1'
'hFEF'};

ROIid = 18:22;

regLab(ROIid)

count=1;
for i=ROIid
    roiCorr{count}=XYZWang(:,atl(:)==i);
    count=count+1;
end

sum_tmp=[];
for i=1:length(roiCorr)
    tmp=roiCorr{i}';
    sum_tmp=[sum_tmp;tmp,ones(size(tmp,1),1)*(i-1)];
end

MNIcorr=sum_tmp;
array2table(MNIcorr)

wangDist=pdist(sum_tmp);
[min(wangDist), max(wangDist)]
wangRes=min(wangDist)

regLab(ROIid)