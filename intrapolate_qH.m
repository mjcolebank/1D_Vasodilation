qh = load('qH1_8192.dat');

plot(qh)

vq = interp1(1:1:8193,qh,1.5:1:8193.5)'; 

plot(1:1:8193,qh,'o',1.5:1:8193.5,vq,':.');
% returns interpolated values of a 1-D function at specific query points 
% using linear interpolation. Vector x contains the sample points,
% and v contains the corresponding values, v(x). Vector xq contains the coordinates of the query points.

qh_2 = zeros(2*8193,1);
qh_2(1:2:end) = qh;
qh_2(2:2:end) = vq;
qh_2(end) = qh_2(end-1);

plot(qh_2);

dlmwrite('qH1_16386.dat',qh_2,'delimiter', ',','precision', '%7.5f');