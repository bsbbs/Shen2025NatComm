function lists = myfnames(pattern,varargin)
fullpath = 1;
for i=1:length(varargin)
        if(ischar(varargin{i}))
            switch(varargin{i})
                case 'single'
                    fullpath = 0;
            end;
        end;
end;
folderstruct = dir(pattern);
lists = {};
for i = 1:numel(folderstruct)
    if fullpath == 1
        lists{i} = fullfile(folderstruct(i).folder,folderstruct(i).name);
    else
        lists{i} = folderstruct(i).name;
    end;
end;
lists = lists';