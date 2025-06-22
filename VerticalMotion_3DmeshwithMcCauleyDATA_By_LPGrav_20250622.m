%2025, June, 22, for 6DoFSVC1_LPGrav
%copied from IN1_Exp
%mccauleydata.xlsxを読み込み，全体をプロットする．
%
%Ferris Wheel Illusion using IN1 model (release)
% last updated April, 10, 2024 for IEEE SMC2024 with BOS
% copied from FerrisWheelIllusionByIN1_20231104.m(deleted)
% & slightlyh modified.

%→　Ferris wheel illusion can be described even in internal model.

% For simulation but input is f(=g+a)
%{
    In-1 model in 
 　1)S. Inoue, H. Liu, and T. Wada, "Revisiting Motion Sickness Models Based on SVC Theory Considering Motion Perception,", SAE Technical paper, 2023. doi: 10.4271/2023-01-0176.
%}
clc
clear

filename = 'mccauleydata.xlsx';

% オプションを明示的に設定
opts = spreadsheetImportOptions('NumVariables', 3);
opts.DataRange = 'A2';  % 2行目から読み込み
opts.VariableNames = {'ACC', 'freqarray', 'MSIdata'};
opts.VariableTypes = {'double', 'double', 'double'};

% 読み込み
T = readtable(filename, opts);

% 各列を変数に代入
Accdata = 9.81*sqrt(2)*T{:,1}; %RMSgからm/s^2振幅に換算
Freqdata = T{:,2};
MSIdata = T{:,3};

DEG2RAD=3.141592/180.0;
grav=9.81;
Dt=0.01;
Fs=floor(1/Dt);
FTIME=7200;

MABIKI = 10;
NTIME = 30*100*Fs;
fileID = fopen('MSIfile.dat','w');

%Parameters used in 　(2024年4月11日に更新した．→2023年10月11月のは誤ってIN3のパラメータ使用
% Inoue, S., Liu, H., Wada, T. Revisiting Motion Sickness Models Based on SVC Theory Considering Motion Perception. 
% SAE Technical Paper, 2023
sixdofSVCparam.Ka=0.1;  
sixdofSVCparam.Kw=0.1;%Kw=0.8
sixdofSVCparam.Kwc=10;
sixdofSVCparam.Kgc=5;  %元々Kvc
sixdofSVCparam.Kac=0.5; %changed
sixdofSVCparam.taug=100;
sixdofSVCparam.taud=7;  
%sixdofSVCparam.taua=190; 
sixdofSVCparam.tau=2;    %changed
sixdofSVCparam.P=85; 
sixdofSVCparam.tauI = 12*60; %[sec] 
sixdofSVCparam.b=0.5;

time = [0:Dt:FTIME]';       % Time data
n_time = length(time);

fmat=zeros(3,n_time);
omegamat=zeros(3,n_time);
MSIarray = zeros(n_time,1); %for File output
tarray = zeros(n_time,1);
%freq=1.2;
%MSImesh=zeros(n_conditions,n_time);

Freqcond = [0.083, 0.167, 0.18, 0.2, 0.25, 0.333, 0.417, 0.5, 0.6,0.7];
Ampcond = [0.0278, 0.055, 0.111, 0.17, 0.222, 0.234, 0.333, 0.444, 0.555];
Ampcond = sqrt(2)*9.81*Ampcond;

n_freqcond = length(Freqcond);
n_ampcond = length(Ampcond);

FinalMSI = zeros(n_freqcond,n_ampcond); ...

for iFreqCond = 1:n_freqcond
    for jAmpCond = 1:n_ampcond
% Vertical motion
    dottheta=0;  %[rad/s]
    theta=0;     %[rad] 0=upright
for i=1:1:n_time
    omegamat(:,i) = [0,0,0]; %[dottheta,0,0];
    fmat(:,i) = [0.0, 0.0, 9.81+Ampcond(jAmpCond)*sin(2*pi*Freqcond(iFreqCond)*time(i))]; 
    %fmat(:,i) = [0.0, 9.81*sin(dottheta*t(i)), 9.81*cos(dottheta*t(i))]; 
end

uinput = [omegamat; fmat]';

  fDummy=[0;0;9.81];
  xscc=[0;0;0];
  x_hscc = [0;0;0];
  v_s=fDummy;  %  v_s = [0;0;9.81];
  hv_s = fDummy;  %  hv_s = [0;0;9.81];
  tildeg = fDummy;
  g = fDummy; %  g = [0;0;9.81];
  dotMsi = 0;
  Msi = 0;

x=[xscc; %1-3
 x_hscc; %4-6
 v_s;   %7-9
 hv_s;  %10-12
 tildeg %13-15
 g;     %16-18
 dotMsi;    %19
 Msi];      %20

%ここにあった「初期値求めるルーチン」消した．
%この時点で状態変数x(0)～x(20)の初期値が定まった．
 
%以降はデータを用いたMSI計算．

xarray=zeros(20,n_time);
conflict=zeros(n_time,1);

for count = 1:n_time
%	uinput = [omegamat(:,count); amat(:,count)]; for simulation

nextx = sixdofSVC_LPGravEXP(x, uinput(count,:), Dt, sixdofSVCparam);
    x = nextx;
    xarray(:,count) = x;
    dumv_s = [x(7),x(8),x(9)];
    dumhv_s = [x(10), x(11), x(12)];
    conflict(count) = norm(dumv_s - dumhv_s);
    tarray(count) = count*Dt;
    %MSImesh(iCond,count) = x(20);
end

FinalMSI(iFreqCond, jAmpCond) = x(20); 
    end
end

% サーフェス描画
[X, Y] = meshgrid(Freqcond, Ampcond);
Z = FinalMSI';  % 転置して [length(amp) × length(freq)]
figure;
surf(X, Y, Z);
xlabel('Frequency (Hz)');
ylabel('Amplitude [m/s2]');
zlabel('MSI');
title('MSI Surface with Highlighted Data Points');
shading interp;
colorbar;
hold on;

% 各 datapoint(i) を (freq(i), amp(i), datapoint(i)) にプロット
for i = 1:length(MSIdata)
    x = Freqdata(i);
    y = Accdata(i);
    z = MSIdata(i);
    plot3(x, y, z, 'ro', 'MarkerSize', 10, 'LineWidth', 2);  % 赤丸マーカー
end

rotate3d on;