% Taken from Martinis (2009), Houzet, Catelani etc.
function y=integmartinis(E,dE,kbT,delta,dmu)
    
    EL = E;
    ER = E+dE;

    rhoL = EL./(sqrt(EL.*EL - delta.^2));
    rhoR = ER./(sqrt(ER.*ER - delta.^2));
    uminusv_sq_term = (1./2).*(EL.*ER - delta.*delta)./EL./ER;

    y = 2.*rhoL.*2.*rhoR.*uminusv_sq_term.*fe(EL-dmu,kbT).*(1-fe(ER-dmu,kbT));
    