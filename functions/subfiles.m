function filelist = subfiles(inputName,incl_path)
    if nargin < 1
        templist = dir;
    else
        templist = dir(inputName);
    end
    if nargin < 2
        incl_path = false;
    else
    end
    tempIdx = strfind(inputName,'/');
    if isempty(tempIdx)
        curDir = pwd;
    else
        curDir = inputName(1:max(tempIdx));
    end
    num_folders = 0;
    for t=1:length(templist)
        if templist(t).isdir ==0 && ~strcmp(templist(t).name,'.') && ~strcmp(templist(t).name,'..')  
            num_folders = num_folders+1;
            if incl_path
                filelist(num_folders,:) = {[curDir,'/',templist(t).name]};
            else
                filelist(num_folders,:) = {templist(t).name};
            end
        else
        end
    end
    if num_folders == 0
        filelist = false;
    else
    end
end

    