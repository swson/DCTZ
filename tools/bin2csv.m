function bin2csv(filename)
% usage: bin2csv('rlds.bin')
%    will generate 'rlds.csv'
    fid  = fopen(filename);
    input = fread(fid,'double');
    fclose(fid);
    token = strtok(filename,'.');
    csvwrite(strcat(token,'.csv'),input);
end