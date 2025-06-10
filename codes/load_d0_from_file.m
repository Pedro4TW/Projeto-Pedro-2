function d0 = load_d0_from_file(filename,beta,gamma,mobility,cong,ngoods)

fid=fopen(filename);
if fid==-1
    error('%s.m: the calibration file %s does not exist.\n',mfilename,filename);
end

data=textscan(fid,'beta=%f gamma=%f mobility=%.1f cong=%d ngoods=%d d0=%f','CommentStyle','%');

I=find(abs(data{1}-beta)<1e-4 & abs(data{2}-gamma)<1e-4 & data{3}==mobility & data{4}==cong & data{5}==ngoods); % identify beta,gamma with precision 1e-4

if isempty(I)
    error('%s.m: cannot find a value for d0 corresponding to beta=%2.4f gamma=%2.4f mobility=%.1f cong=%d ngoods=%d in %s.\n',mfilename,beta,gamma,mobility,cong,ngoods,filename);
end

d0=data{6}(I(1));

fclose(fid);

