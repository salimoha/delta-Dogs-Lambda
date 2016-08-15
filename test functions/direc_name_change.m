d=dir('*.log');
for k=1:length(d)
    fname=d(k).name;
    [pathstr, name, ext] = fileparts(fname);
    movefile(fname, fullfile(pathstr, [name '.txt']))
    % ...
end