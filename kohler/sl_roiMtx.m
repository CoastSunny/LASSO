function ROIs = sl_roiMtx(subjId)
    if strfind(subjId,'_')
        fsSubj=[subjId(1:(strfind(subjId,'_')-1)) '_fs4'];
    else
        fsSubj=[subjId,'_fs4'];
    end
    anatDir = getpref('mrCurrent','AnatomyFolder');
    roiFile = fullfile(anatDir,subjId,'Standard','meshes','ROIs_correlation.mat');
    load(roiFile);
end

