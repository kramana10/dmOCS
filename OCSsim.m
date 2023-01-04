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

fermi = 11.6; % eV Fermi level
Tc = 1.2; % AlMn 0.1, Al 1.2 

% Superconducting gap for OCS and also for "left" and "right" sides of QCD
delta = 2E-4; % eV
deltaL = 2.2E-4; % eV 
deltaR = 1.8E-4; % eV

%% Experiment Fixed

T = 100E-3; %K Temperature
Rn = 20.0E3; % ohms normal state resistance
cpl = 6.5e-17;% F/um gate capacitance per micron
clen = 20; % um gate capacitor length in micron
Cg = 1E-15;%clen*cpl; % F gate capacitance

%% Computed quantities

e2R = ec*ec*Rn/eVtoJ; % I^2*R*t*t = e2R % E*t --> eV*s

% Transmon
Cj = 4.43E-16; % F junction capacitance
Cshunt = 5.30E-14;
Csigma = Cj+Cg+Cshunt; % F total capacitance
[Ec, Ej] = computeEcEj(Csigma, delta, Rn);
curlyN = DOS*sqrt(2*pi*delta*kb*T); % /um3

disp(['Transmon Ej/Ec ~ ' num2str(Ej/Ec)])

%% Energy diagram
% Diagonalizes Hamiltonian to generate qubit energy levels

u=linspace(0,1,500); % Offset Charge
[EE, EO, DE] = solvesystem(Ec,Ej,u, delta, delta);

EE_f = (EE - min(EE(:,1)))./planck./1e9;
EO_f = (EO - min(EO(:,1)))./planck./1e9;

FigHandleA = figure;
set(FigHandleA, 'Position', [100, 100, 800, 600]);
%plot(u,EE,'LineWidth',2);
plot(u,EE_f,'LineWidth',2);
ax = gca;
ax.ColorOrderIndex = 1;
hold on
%plot(u,EO,'LineWidth',2,'LineStyle','--');
plot(u,EO_f,'LineWidth',2,'LineStyle','--');
hold off
xlabel('Offset Charge [CgVg/2e]','Interpreter','latex','FontSize',25);
ylabel('f$_{0j}$ [GHz]','Interpreter','latex','FontSize',25);
title(['$E_J/E_C=$ ' num2str(Ej/Ec)],'Interpreter','latex','FontSize',25);
set(gca,'TickLabelInterpreter','latex','FontSize',25);
set(gca,'YMinorTick','on');
set(gca,'XMinorTick','on');

disp(['Transmon dE/Ec ratio: ' num2str(DE(1)/Ec)]);
fr3 = min(EO_f(:,4));

%% Tunneling Rates

% We are going to ignore PAPS or coupling between 0, 1 states etc. Assume
% that "thermal current" term dominates. Might be dangerous in PAPS case.

nqp = 1; % /um3

% Transmon Case
% Tunneling in & out set by DOS integral
% Assume ~same volume absorbers and island

dmu = kb.*T.*log(1+(nqp./curlyN).*exp(delta/kb./T)); % assume chemical potential shift in Right is same.
gamma_in_transmon = (1./e2R).*quadgk(@(E)integmartinis(E,DE(1),kb.*T,delta,dmu),delta,5.*delta);
gamma_in_transmon_occupied =(1./e2R).*quadgk(@(E)integmartinis(E,-DE(1),kb.*T,delta,dmu),delta+DE(1),5.*delta);

tun_ratio = (gamma_in_transmon_occupied)*100/gamma_in_transmon;

disp('-----');
disp(['nqp (density) assumed: ' num2str(nqp) ' /um3'])
disp(['K Transmon: ' num2str(gamma_in_transmon) ' Hz'])
disp(['2nd qp Transmon ratio: ' num2str(tun_ratio) ' %'])

%% Dispersive Shift

g = 200E6; % MHz example
fr = fr3.*1E9 - 50E6; % 3rd level minus 50 MHz
nlevels = 5;
matrixelems = zeros(length(u),nlevels);
chi_ip = zeros(length(u),nlevels);
for i=1:length(u)
    [fullmat,chi_ip(i,:)]=dispermatrix(Ec,Ej,u(i)+0.5,g,fr,nlevels);
    matrixelems(i,:) = fullmat(1,:);
end

fr_range = fr:1E6:(fr+2E8);
chi_scan = zeros(1,length(fr_range));
ind = 1;
for fr_span = fr_range % scan 200 MHz around it 
    [~,chi_scan_1]=dispermatrix(Ec,Ej,u(1)+0.5,g,fr_span,nlevels);
    [~,chi_scan_2]=dispermatrix(Ec,Ej,u(250)+0.5,g,fr_span,nlevels);
    chi_scan(ind) = chi_scan_1(1) - chi_scan_2(1);
    ind = ind+1;
end

FigHandleB = figure;
set(FigHandleB, 'Position', [100, 100, 800, 600]);
semilogy(u(1:250),matrixelems(1:250,2:end),'LineWidth',2,'LineStyle','-');
ylim([1E-5 2E0]);
xlabel('Offset Charge [CgVg/2e]','Interpreter','latex','FontSize',25);
ylabel('$\langle j,o|\hat{n}|0,o\rangle$','Interpreter','latex','FontSize',25);
set(gca,'TickLabelInterpreter','latex','FontSize',25);
set(gca,'YMinorTick','on');
set(gca,'XMinorTick','on');
legendCell = cellstr(num2str([1:nlevels-1]','j=%-d'));
leg = legend(legendCell,'location','best','Interpreter','latex','FontSize',25);
legend box off

FigHandleC = figure;
set(FigHandleC, 'Position', [100, 100, 800, 600]);
plot(u(1:125),chi_ip(1:125,1:2)./1e6,'LineWidth',2,'LineStyle','-');
hold on
plot(u(250:375)-u(250),chi_ip(250:375,1:2)./1e6,'LineWidth',2,'LineStyle','-');
hold off
xlim([0 0.25]);
xlabel('Offset Charge [CgVg/2e]','Interpreter','latex','FontSize',25);
ylabel('$\chi_{i,p} [\frac{g^2}{200^2\,\rm{MHz}}]$ [MHz]','Interpreter','latex','FontSize',25);
set(gca,'TickLabelInterpreter','latex','FontSize',25);
set(gca,'YMinorTick','on');
set(gca,'XMinorTick','on');
leg = legend({'$|0,o\rangle$','$|1,o\rangle$','$|0,e\rangle$','$|1,e\rangle$'},'location','best','Interpreter','latex','FontSize',25);
legend box off

FigHandleD = figure;
set(FigHandleD, 'Position', [100, 100, 800, 600]);
plot(fr_range./1e9,chi_scan./1e6,'LineWidth',2,'LineStyle','-');
xlabel('Frequency [GHz]','Interpreter','latex','FontSize',25);
ylabel('$\Delta\chi_{0} [\frac{g^2}{200^2\,\rm{MHz}}]$ [MHz]','Interpreter','latex','FontSize',25);
set(gca,'TickLabelInterpreter','latex','FontSize',25);
set(gca,'YMinorTick','on');
set(gca,'XMinorTick','on');
