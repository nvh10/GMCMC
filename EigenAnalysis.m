function [kR, phix, phiz, kL, phiy] = EigenAnalysis(Ax, Ay, Az, Gx, Gy, Gz, Mx, My, Mz, Bxz, nelmt, ww)

kR = zeros(2*nelmt,1);
phix = zeros(nelmt,2*nelmt);
phiz = zeros(nelmt,2*nelmt);

kL = zeros(nelmt,1);
phiy = zeros(nelmt,nelmt);
    
Ab = [Ax(1:nelmt,1:nelmt) zeros(nelmt)
    transpose(Bxz(1:nelmt,1:nelmt)) Az(1:nelmt,1:nelmt)];

Cb = [Gx(1:nelmt,1:nelmt)-ww^2*Mx(1:nelmt,1:nelmt) Bxz(1:nelmt,1:nelmt)
    zeros(nelmt) Gz(1:nelmt,1:nelmt)-ww^2*Mz(1:nelmt,1:nelmt)];


Z = zeros(2*nelmt);
Y = zeros(2*nelmt);

% [Z,D] = eig(Cb,-Ab);
[Z,D] = eig(full(Cb),full(-Ab));

for itmp = 1 : 2 * nelmt
    
    kR(itmp) = sqrt(D(itmp,itmp));

    if imag(kR(itmp)) > 0
        
        kR(itmp) = -kR(itmp);
        
    end
    
    phix(:,itmp) = Z(1:nelmt,itmp);
    phiz(:,itmp) = Z(nelmt+1:2*nelmt,itmp) / kR(itmp);
    
    Y(1:nelmt,itmp) = kR(itmp) * phix(:,itmp);
    Y(nelmt+1:2*nelmt,itmp) = phiz(:,itmp);
        
    fac = transpose(Y(:,itmp)) * Ab * Z(:,itmp);
    fac = sqrt(kR(itmp) / fac);
    
    phix(:,itmp) = fac * phix(:,itmp);
    phiz(:,itmp) = fac * phiz(:,itmp);
    
end

Z = zeros(nelmt);
Y = zeros(nelmt);

Ab = Ay(1:nelmt,1:nelmt);

Cb = Gy(1:nelmt,1:nelmt) - ww^2 * My(1:nelmt,1:nelmt);

% [Z,D] = eig(Cb,-Ab);
[Z,D] = eig(full(Cb),full(-Ab));

for itmp = 1 : nelmt
    
    kL(itmp) = sqrt(D(itmp,itmp));

    if imag(kL(itmp)) > 0
        
        kL(itmp) = -kL(itmp);
        
    end
    
    phiy(:,itmp) = Z(1:nelmt,itmp);
    
    fac = transpose(phiy(:,itmp)) * Ab * phiy(:,itmp);
    
    phiy(:,itmp) = sqrt(1/fac) * phiy(:,itmp);

end

end