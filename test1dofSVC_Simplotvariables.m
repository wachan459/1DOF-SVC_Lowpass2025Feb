%2025, Feb, 2
%Copied from onedofSVC_BosBles1998
%2025, Jan. 30
%Copied Jan, 30, 2025 from test 6dofSVC_In1_Simplotvariables.m (2023, June,5 version)
% For simulation but input is f(=g+a)
%{
 Bos & Bles, 1998
%}

clc
clear

DEG2RAD=3.141592/180.0;
Dt=0.01;
Fs=floor(1/Dt);
FTIME=7200;

MABIKI = 10;
NTIME = 30*60*Fs;
fileID = fopen('MSIfile.dat','w');


onedofSVCparam.Kgc=5;
onedofSVCparam.tau=5;
onedofSVCparam.Kac=1;
onedofSVCparam.tauL = 7200; %重力へのアダプテーション．適当．
onedofSVCparam.P=85; 
onedofSVCparam.tauI = 12*60; %[sec] 
onedofSVCparam.b=0.7;


%Newlly add Jun 5, 2023, for simulation.

t = [0:Dt:FTIME]';       % Time data
n_time = length(t);

amat=zeros(1,n_time);
MSIarray = zeros(n_time,1); %for File output
tarray = zeros(n_time,1);
freq=0.2;

%入力準備
for i=1:1:n_time
    grav=9.81;
    amat(i) = grav + sqrt(2)*0.5*9.81*sin(2*pi*freq*t(i)); %0.3G(RMS)
end
uinput = [amat]';

%xstateの初期値代入（自動キャリブレーションは解除）
x=zeros(5,1);
x(1) = 9.81; %g_s
x(2) = 9.81; %hg_s
x(3) = 9.81; %tldg
x(4) = 0; %dotMSI
x(5) = 0; %MSI

%以降はデータを用いたMSI計算．
% From here on, MSI calculations will be performed using the data.

xarray=zeros(5,n_time);
conflict=zeros(n_time,1);

for count = 1:n_time
nextx = onedofSVC_Lowpass_EXP2025Feb(x, uinput(count,:), Dt, onedofSVCparam);
    x = nextx;
    xarray(:,count) = x;
    conflict(count) = norm(x(1) - x(2));
    MSIarray(count) = x(5);
    tarray(count) = count*Dt;

    if rem(count, MABIKI) == 0
        fprintf(fileID,'%f, %f \n', count*Dt, x(5));
    end

end

f = uinput(:);
g_s = xarray(1,:)';
hg_s = xarray(2,:)';
tldg = xarray(3,:)';
dotMsi = xarray(4,:)';
Msi=xarray(5,:)';
prm=onedofSVCparam;

hf = (tldg+(prm.Kgc-prm.Kac)*g_s-(prm.Kgc-prm.Kac)*hg_s +prm.Kac*f)/(1+prm.Kac);

time=tarray;

Tbl = table(time, f, g_s, hg_s, tldg, hf, conflict, Msi);

sp= stackedplot(Tbl, 'XVariable', 'time');
sp.AxesProperties(1).YLimits = [0 20]; % 1つ目のプロットのy軸範囲を設
sp.AxesProperties(2).YLimits = [-10 10]; % 1つ目のプロットのy軸範囲を設
sp.AxesProperties(3).YLimits = [-10 10]; % 1つ目のプロットのy軸範囲を設
sp.AxesProperties(4).YLimits = [0 20]; % 1つ目のプロットのy軸範囲を設
sp.AxesProperties(5).YLimits = [0 20]; % 1つ目のプロットのy軸範囲を設
sp.AxesProperties(6).YLimits = [0 20]; % 1つ目のプロットのy軸範囲を設
sp.AxesProperties(7).YLimits = [0 100]; % 1つ目のプロットのy軸範囲を設


MSIatFinaltime = Msi(end);

