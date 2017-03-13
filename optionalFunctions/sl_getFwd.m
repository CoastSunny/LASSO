function fwd = sl_getFwd(projectDir,subjId)
    if strfind(subjId,'_')
        fsSubj=[subjId(1:(strfind(subjId,'_')-1)) '_fs4'];
    else
        fsSubj=[subjId,'_fs4'];
    end
    freesurfDir = getpref('freesurfer','SUBJECTS_DIR');
    srcFileName = fullfile(freesurfDir,fsSubj,'bem',[ fsSubj '-ico-5p-src.fif']);
    srcSpace = mne_read_source_spaces(srcFileName);
    
    mneFwdFileName = [subjId '-fwd.fif' ];
    mneFwdFile = fullfile(projectDir,subjId,'_MNE_',mneFwdFileName);
    mneFwd = mne_read_forward_solution(mneFwdFile);
    fwd = makeForwardMatrixFromMne(mneFwd,srcSpace);
end

