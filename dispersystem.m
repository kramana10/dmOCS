%function [matrixelems, chi_ip] = dispersystem(Ec,Ej)
%clear, clc, close all

Ec = planck*356E6;
Ej = planck*6.14E9;
%Ec = planck*190E6;
%Ej = planck*7.81E9;

    u = linspace(0,2,1000); % offset charge
    g = 200E6; % MHz example
    wr = 9E9; % GHz example
    nlevels = 5;
    matrixelems = zeros(length(u),nlevels);
    chi_ip = zeros(length(u),nlevels);
    for i=1:length(u)
        [fullmat,chi_ip(i,:)]=dispermatrix(Ec,Ej,u(i)+0.5,g,wr,nlevels);
        matrixelems(i,:) = fullmat(1,:);
    end
    %matrixelems=matrixelems(:,2:end);

    FigHandleA = figure;
    set(FigHandleA, 'Position', [100, 100, 1400, 600]);
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

    FigHandleB = figure;
    set(FigHandleB, 'Position', [100, 100, 1400, 600]);
    plot(u(1:125),chi_ip(1:125,1:2)./1e6,'LineWidth',2,'LineStyle','-');
    hold on
    plot(u(250:375)-u(250),chi_ip(250:375,1:2)./1e6,'LineWidth',2,'LineStyle','-');
    hold off
    %ylim([0 4]);
    xlim([0 0.25]);
    xlabel('Offset Charge [CgVg/2e]','Interpreter','latex','FontSize',25);
    ylabel('$\chi_{i,p} [\frac{g^2}{200^2\,\rm{MHz}}]$ [MHz]','Interpreter','latex','FontSize',25);
    set(gca,'TickLabelInterpreter','latex','FontSize',25);
    set(gca,'YMinorTick','on');
    set(gca,'XMinorTick','on');
    leg = legend({'$|0,o\rangle$','$|1,o\rangle$','$|0,e\rangle$','$|1,e\rangle$'},'location','best','Interpreter','latex','FontSize',25);
    legend box off
