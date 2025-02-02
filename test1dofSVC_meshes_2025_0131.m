%2025, Feb,2
%Test for onedofSVC_lowpass
%2025, Jan. 31
%Mesh calculation imitating Bos&Bles(1998)
%
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

fileID = fopen('MSIfile.dat','w');


%Parameters used in 
% Inoue, S., Liu, H., Wada, T. Revisiting Motion Sickness Models Based on SVC Theory Considering Motion Perception. 
% SAE Technical Paper, 2023
onedofSVCparam.Kgc=1; %changed from 5
onedofSVCparam.tau=2; %changed from 5 
onedofSVCparam.P=85; 
onedofSVCparam.tauI = 12*60; %[sec] 
onedofSVCparam.b=0.7;
onedofSVCparam.Kac=1;
onedofSVCparam.tauL = 72; %�d�͂ւ̃A�_�v�e�[�V�����D�K�� changed from 7200�D

%Newlly add Jun 5, 2023, for simulation.

t = [0:Dt:FTIME]';       % Time data
n_time = length(t);

fmat=zeros(1,n_time);
MSIarray = zeros(n_time,1); %for File output
tarray = zeros(n_time,1);
%freq=0.2;

%0.1 - 1.0 Hz (10 levesl)
%0.1 - 0.4 RMSG (4 levels)

%Freq=[0.15, 0.2];
Freq=[0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
Amp=[0.05, 0.1, 0.2, 0.3, 0.4]
MSI_Final_Mat=zeros(length(Freq), length(Amp))

for ifreq=1:length(Freq)
    for jamp=1:length(Amp)
    
%���͏���
for i=1:1:n_time
    grav=9.81;
    fmat(i) = grav + sqrt(2)*Amp(jamp)*9.81*sin(2*pi*Freq(ifreq)*t(i)); %0.3G(RMS)
end
uinput = [fmat]';

%xstate�̏����l����i�����L�����u���[�V�����͉����j
x=zeros(5,1);
x(1) = 9.81; %g_s
x(2) = 9.81; %hg_s
x(3) = 9.81; %hf
x(4) = 0; %dotMSI
x(5) = 0; %MSI

%�ȍ~�̓f�[�^��p����MSI�v�Z�D
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

end

f = uinput(:);
g_s = xarray(1,:)';
hg_s = xarray(2,:)';
hf = xarray(3,:)';
dotMsi = xarray(4,:)';
Msi=xarray(5,:)';


MSI_Final_Mat(ifreq,jamp) = Msi(end);
    end
end

Z = MSI_Final_Mat;
%[xSize, ySize] = size(Z);
%[X, Y] = meshgrid(1:ySize, 1:xSize); % X, Y ���쐬


% X, Y �̃��b�V���O���b�h�쐬
[X, Y] = meshgrid(Amp, Freq); % X=�U��, Y=���g��

% 3D���b�V���v���b�g
figure;
mesh(X, Y, MSI_Final_Mat);

% �����x��
xlabel('Amplitude RMS[G]', 'FontSize', 14, 'FontWeight', 'bold'); % X�� (Amp)
ylabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold'); % Y�� (Freq)
zlabel('MSI Value', 'FontSize', 14, 'FontWeight', 'bold'); % Z�� (MSI_Final_Mat)

% �^�C�g��
title('3D Mesh Plot of MSI Values', 'FontSize', 16, 'FontWeight', 'bold');

% �J���[�o�[�ǉ�
colorbar;


