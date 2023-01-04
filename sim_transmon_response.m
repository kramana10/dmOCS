clear, clc, close all

%% Universal Constants

planck = 4.135E-15; % eV.s
kb = 8.617E-5; % eV/K
eVtoJ = 1.6E-19;
ec = 1.6E-19;  % coulomb

%% Aluminum specific
DOS = 5; % /Ry/Cell in weird units
cellV = 0.0662/1E9; % um^3 volume of unit cell
DOS = DOS/13.6/cellV; % /eV/um3 DOS at Fermi level in the normal state 
% Superconducting gap for OCS and also for "left" and "right" sides of QCD
delta = 2E-4; % eV

%% Experiment Fixed
T = 100E-3; %K Temperature
Rn = 20E3; % ohms normal state resistance
cpl = 6.5e-17;% F/um gate capacitance per micron
clen = 20; % um gate capacitor length in micron
Cg = clen*cpl; % F gate capacitance
etaph = 0.3;
tqp = 0.5E-3; %s
tabs = 1E-3; %s

% Transmon
Cj = 1E-15; % F junction capacitance
Cshunt = 40E-15;
Csigma = Cj+Cg+Cshunt; % F total capacitance
[Ec, Ej] = computeEcEj(Csigma, delta, Rn);
curlyN = DOS*sqrt(2*pi*delta*kb*T); % /um3
disp(['Transmon Ej/Ec ~ ' num2str(Ej/Ec)]);

%% Calculate rates etc.

n0 = 0.02; %um-3
F=0.2; %Fano factor

E_dep = 0.01:0.01:100; % eV
Nqp = E_dep.*etaph./delta;
V = [100 200 500 1000 2000 5000 10000 30000];
[X, Y] = meshgrid(Nqp,V);

K = 16.*Ej.*kb.*T./curlyN./delta./planck; % Hz um3
fakeK = 2.*kb.*T./curlyN./ec./Rn;
disp(['K: ' num2str(K) ' Hz.um3']);

% Signal
S = K.*X.*tqp./Y;

% Background
tau = 4.*tabs;
B = K.*n0.*tau;

% Sigma
sigma_B = sqrt((K.*tau+16.*K.*K.*tau.*tau./Y).*n0);
sigma_S = sqrt((K.*tqp./Y + K.*K.*tqp.*tqp.*F./Y./Y).*X);
sig_tot = max(1,sqrt(sigma_B.^2 + sigma_S.^2));

%Threshold
[val, idx] = min(abs(S+B-(B+3.*sigma_B)),[],2);
thrsh = K.*Nqp(idx).*tqp./V + B;

Bene = B*delta/tqp/etaph/K;

%% Plots

FigHandleA = figure;
set(FigHandleA, 'Position', [100, 100, 1400, 600]);
subplot(1,2,1)
loglog(E_dep,S+B,'LineWidth',2);
ax = gca;
ax.ColorOrderIndex = 1;
hold on
loglog(E_dep(idx),thrsh,'k--','LineWidth',2);
%loglog(E_dep,3.*sig_tot,'LineWidth',1);
hold off
xlabel('Deposited Energy [eV]','Interpreter','latex','FontSize',18);
ylabel('E[flips]','Interpreter','latex','FontSize',18);
set(gca,'TickLabelInterpreter','latex','FontSize',14);
set(gca,'YMinorTick','on','XMinorTick','on');
legendCell = cellstr(num2str(V', '%-d $um^3$'));
legendCell{length(V)+1}='3$\sigma_B$ threshold';
leg = legend(legendCell,'location','best','Interpreter','latex','FontSize',12);
legend box off

subplot(1,2,2)
loglog(E_dep,sig_tot.*delta.*Y./K./tau./etaph,'LineWidth',2);
hold on
plot(E_dep,Bene+E_dep./3,'k--','LineWidth',2);
xlabel('Energy [eV]','Interpreter','latex','FontSize',18);
ylabel('$\sigma_{\rm tot}$ (eV)','Interpreter','latex','FontSize',18);
set(gca,'TickLabelInterpreter','latex','FontSize',14);
set(gca,'YMinorTick','on');
set(gca,'XMinorTick','on');
legendCell = cellstr(num2str(V', '%-d $um^3$'));
leg = legend(legendCell,'location','best','Interpreter','latex','FontSize',12);
legend box off
