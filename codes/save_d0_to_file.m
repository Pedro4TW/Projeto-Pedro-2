function save_d0_to_file(filename,beta,gamma,mobility,cong,ngoods,d0)

% Open file and retain existing data
fid=fopen(filename,'rt');
if fid~=-1 % if file exists
    data=textscan(fid,'beta=%f gamma=%f mobility=%.1f cong=%d ngoods=%d d0=%f','CommentStyle','%');
    fclose(fid);    
else
    data=cell(1,6);
end

% Update with new value
I=find(abs(data{1}-beta)<1e-4 & abs(data{2}-gamma)<1e-4 & data{3}==mobility & data{4}==cong & data{5}==ngoods); % identify beta,gamma with precision 1e-4

if isempty(I)
    data{1} = [data{1};beta];
    data{2} = [data{2};gamma];
    data{3} = [data{3};mobility];
    data{4} = [data{4};cong];
    data{5} = [data{5};ngoods];
    data{6} = [data{6};d0];
else
    data{6}(I(1))=d0;
end

% Write to file    
fid=fopen(filename,'wt'); % create or overwrite the file
fprintf(fid,'%% This file was created by %s.m. It contains the various calibrated values of d0 for values of (beta,gamma,mobility,cong,ngoods).\n',mfilename);

for n=1:length(data{1})
    fprintf(fid,'beta=%f gamma=%f mobility=%.1f cong=%d ngoods=%d d0=%f\n',data{1}(n),data{2}(n),data{3}(n),data{4}(n),data{5}(n),data{6}(n));
end

fclose(fid);

