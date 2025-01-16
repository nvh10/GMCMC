function [Uz] = SurfaceDisplacementFreqDm(nelmt, la, mu, rho, sdamp, h, etype, ww, rr, RDisk, qDisk)

nfreq = size(ww,2);
Nmeasured = size(rr,2);

Ur = zeros(nfreq,Nmeasured);
Uz = zeros(nfreq,Nmeasured);

for ifreq = 1 : nfreq
    
    [Ax, Ay, Az, Gx, Gy, Gz, Mx, My, Mz, Bxz] = DynamicStiffness(nelmt, la, mu, rho, sdamp, h, etype, ww(ifreq));

    [kR, phix, phiz, kL, phiy] = EigenAnalysis(Ax, Ay, Az, Gx, Gy, Gz, Mx, My, Mz, Bxz, nelmt, ww(ifreq));
    
    for ir = 1 : Nmeasured
    
        for ieig = 1 : 2 * nelmt

            if (rr(ir) < RDisk)

                I1R = pi / (2*1i*kR(ieig)) * besselj(0,kR(ieig)*rr(ir)) * besselh(1,2,kR(ieig)*RDisk) - 1 / (RDisk*kR(ieig)^2);
                I2R = pi / (2*1i) * besselj(1,kR(ieig)*rr(ir)) * besselh(1,2,kR(ieig)*RDisk);
                
            end
            
            if (rr(ir) >= RDisk)
                
                I1R = pi / (2*1i*kR(ieig)) * besselj(1,kR(ieig)*RDisk) * besselh(0,2,kR(ieig)*rr(ir));
                I2R = pi / (2*1i) * besselj(1,kR(ieig)*RDisk) * besselh(1,2,kR(ieig)*rr(ir));
                
            end
            
            Ur(ifreq,ir) = Ur(ifreq,ir) + qDisk * RDisk * phix(1,ieig) * phiz(1,ieig) * I2R / kR(ieig);
            Uz(ifreq,ir) = Uz(ifreq,ir) + qDisk * RDisk * phiz(1,ieig) * phiz(1,ieig) * I1R;
                    
        end
        
    end
    
end