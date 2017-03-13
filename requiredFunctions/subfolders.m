function folderlist = subfolders(inputName,incl_path)
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
        curDir = inputName(1:(max(tempIdx)-1));
    end
    num_folders = 0;
    if ~isempty(templist)
        for t=1:length(templist)
            if templist(t).isdir ==1 && ~strcmp(templist(t).name,'.') && ~strcmp(templist(t).name,'..')  
                num_folders = num_folders+1;
                if incl_path
                    folderlist(num_folders,:) = {[curDir,'/',templist(t).name]};
                else
                    folderlist(num_folders,:) = {templist(t).name};
                end
            else
            end
        end
    else
        folderlist = [];
    end
end

    