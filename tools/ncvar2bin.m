function ncvar2bin(filename, varname)
% usage: ncvar2bin('cellular.nc', 'pres')
%    will generate 'pres.bin'
  vinfo = ncinfo(filename, varname);
  vardata = ncread(filename, varname);
  foutname = strcat(varname, '.bin');
  fid = fopen(foutname, 'W');
  fwrite(fid, vardata, vinfo.Datatype);
  fclose(fid);
end
