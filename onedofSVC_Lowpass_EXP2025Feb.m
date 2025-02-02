%2025, Feb, 2, 1DoF-SVC with gravity outside version
%2025, Jan. 30 BOS's 1DOF-SVC model (Bos & Bles, 1998)
%%2025, Jan. 30 copied from sixdofSVC_in1_exp.m #From GitHub
%
% Vestibular input only. 
% For experimental data
% To accept experimental data, 
% 1) f (=g+a) is input
% 2) initialization algorithm of x0 was deleted.
 
% xnext = onedofSVC_BosBles1998_EXP(xstate, uinput, Dt, prm)

%[OUTPUT]
% xnext = state at next time step (x_i+1) in 5 dimension.
%
%[INPUT]
% xstate = current state (x_i)
%{
	g_s = xstate(1);
	hg_s = xstate(2);
    hf = xstate(3);
	dotMsi = xstate(4);
	Msi=xstate(5); 
%}
% uinput = current input (u_i) 1 dim.
%           
% Dt = sampling time [s]
% prm: Struct for model parameters


function	xnext = onedofSVC_Lowpass_EXP2025Feb(xstate, uinput, Dt, prm)

% uinput=[f]; 
    g_s = xstate(1);
	hg_s = xstate(2);
    tldg = xstate(3);
	dotMsi = xstate(4);
	Msi=xstate(5); 
    f = uinput'; % for data
%a = f - g; % % for data

    dumm = (norm(g_s-hg_s)/prm.b)^2;
    x_H = dumm/(1+dumm);

fnc = [
(-g_s + f)/prm.tau;
((prm.Kgc-prm.Kac)*g_s - (1+prm.Kgc)*hg_s + tldg + prm.Kac*f )/prm.tau/(1+prm.Kac );
(-tldg + g_s)/prm.tauL
(-2*prm.tauI*dotMsi - Msi + prm.P*x_H)/prm.tauI/prm.tauI;
dotMsi];

dumxnext = xstate + fnc*Dt;
xnext = dumxnext;

end

