% Run the 1D ST code in the 21 vessel mouse geometry
% clear; close all;
% !chmod +x sor06
% !make clean
!make
% constant f3
f1   = 0;%7e+4;%5e+6;
f2   = -10; %10
f3   = 6e4;%2.5e4;%1e+4;%8e+4;
fs1  = f1;%5e+6;%f1;
fs2  = f2;%-20;%f2;
fs3  = f3;%*10;%1e+6;%f3;
Z0  = 0;%1e2;

% Exponential stiffness
% f1   = 7e+4;%5e+6;
% f2   = -10; %10
% f3   = 1e4;%2.5e4;%1e+4;%8e+4;
% fs1  = f1;%5e+6;%f1;
% fs2  = f2;%-20;%f2;
% fs3  = f3;%*10;%1e+6;%f3;
% Z0  = 0;%1e2;

alpha = 0.88; %Alpha
beta  = 0.68; %Beta
rm   = 0.005;

% For single LRR
lrr  = 17;% 17;%10; %20
% For exponential function, if you use it
% lrr1  = 13.39;
% lrr2  = -0.007708;
vaso = 1.0; %constrict<1, dilate>1
        
%% From the mouse paper
% Structure of the network: which vessels are connected to which, based on
% the mouse geometry
conn = [1 2 3 0
        2 4 5 0
        3 6 7 0
        4 8 9 0
        5 0 0 0
        6 14 15 0
        7 0  0  0
        8 10 11 0
        9  0  0  0
        10 12 13 0
        11 0  0  0
        12 0  0  0
        13 0  0  0
        14 16 17 0
        15 0  0  0
        16 18 19 0
        17 0  0  0
        18 20 21 0
        19 0  0  0
        20 0  0  0
        21 0  0  0];
    
% Indices for which vessels should becoupled to structured tree bondaru
% conditions
terminal = [5 7 9 11 12 13 15 17 19 20 21];
%%
% The dimensions of the geometry:
% Order is [Length Radius_in Radius_out]
% If you want to dilate/narrow, apply to both these radii values, e.g. you
% could write R_new = R_in(terminal).*(VASODILATION)

% CONTROL GEOMETRY

% L = [0.410 0.445 0.372 0.241 0.052 0.202 0.212 0.311 0.177 ...
%      0.262 0.069 0.140 0.062 0.081 0.184 0.083 0.302 0.469 ...
%      0.177 0.178 0.055];
% R_in = [0.047 0.026 0.037 0.024 0.013 0.032 0.017 0.023 0.017 ...
%         0.020 0.016 0.015 0.014 0.026 0.019 0.025 0.015 0.024 ...
%         0.015 0.022 0.018];
% R_out = R_in;

% HYPOXIC GEOMETRY
L = [0.358 0.403 0.308 0.292 0.065 0.160 0.093 0.206 0.051 ...
     0.237 0.088 0.127 0.051 0.120 0.155 0.071 0.168 0.355...
     0.186 0.224 0.107];
R_in = [0.051 0.026 0.037 0.025 0.017 0.028 0.019 0.024 0.017 ...
        0.022 0.017 0.019 0.015 0.027 0.019 0.026 0.018 0.024 ...
        0.018 0.023 0.019];
R_out = R_in;

%
dim_mat  = [L' R_in' R_out'];
%%
tot_ves = max(conn(:));
tot_term = length(terminal);
write_conn = max(conn-1,0);



dlmwrite('connectivity.txt',write_conn,'\t');
dlmwrite('terminal_vessels.txt',terminal-1,'\t');
dlmwrite('Dimensions.txt',round(dim_mat,3),'\t');

pars = [f1 f2 f3 fs1 fs2 fs3 alpha beta lrr rm Z0 vaso];
pars_str = mat2str(pars);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
out = unix(sprintf('sor06.exe  %s',pars_str(2:end-1)));
toc
%%

if out == 0
    fname = strcat('pu_ALL.2d');
    data = load(fname);
    [t,x,p,q,a,c] = gnuplot(data);
end

figure; plot(p);
figure; plot(q);


save('pqa','p','q','a','dim_mat','pars','terminal')
