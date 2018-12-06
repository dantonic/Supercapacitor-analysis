function [mjerenje]=readAL(fname,I)
fid=fopen(fname);
% skip header
fgetl(fid);
fgetl(fid);

% read data
m = fscanf(fid, '%f %f', [2 inf])';  % time, voltage
len = size(m,1);
mjerenje = zeros(len,3);    % time, voltage, current
mjerenje(:,1) = m(:,1) - m(1,1);    % time begin at 0
mjerenje(:,2) = m(:,2);
mjerenje(:,3) = I;

% moving average
%a = 1;
%b(1:5) = 1/5;
%mjerenje(:,2) = filter(b,a,mjerenje(:,2));


fclose(fid);
