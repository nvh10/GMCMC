function [Uz_estimated, dUz_estimated] = ForwardModel(theta, nelmt, nu, rho, sdamp, h, etype, ww, rr, RDisk, qDisk, CR, facUz)
    
Cs = 100*[theta(1:nelmt-4).*ones(nelmt-4,1)
    theta(nelmt-3)*ones(4,1)];

mu = rho .* Cs.^2;
la = 2 * mu .* nu ./ (1 - 2 .* nu);

nfreq = size(ww,2);
Nmeasured = size(rr,2);

Const = nelmt - 3;

Uz_estimated = zeros(Nmeasured,nfreq);
dUz_estimated = zeros(Nmeasured,nfreq,nelmt-3);

parfor ifreq = 1 : nfreq
     kR_rr = zeros(2*nelmt, Nmeasured);
    [Ax, Ay, Az, Gx, Gy, Gz, Mx, My, Mz, Bxz] = DynamicStiffness(nelmt, la, mu, rho, sdamp, h, etype, ww(ifreq));

    [kR, phix, phiz, kL, phiy] = EigenAnalysis(Ax, Ay, Az, Gx, Gy, Gz, Mx, My, Mz, Bxz, nelmt, ww(ifreq));
    
 for ii = 1:Nmeasured
        kR_rr(:, ii) = kR(:) * rr(ii);
    end
    kR_RDisk = kR * RDisk;
    bj1_kR_RDisk = besselj(1, kR_RDisk);
    bh02_kR_RDisk = besselh(0, 2, kR_RDisk);
    bj2_kR_RDisk = besselj(2, kR_RDisk);
    bh02_kR_rr = besselh(0, 2, kR_rr);
    bh12_kR_rr = besselh(1, 2, kR_rr);

    for ir = 1 : Nmeasured

        for ieig = 1 : 2 * nelmt

            if (rr(ir) < RDisk)

                I1R = pi / (2*1i*kR(ieig)) * besselj(0,kR(ieig)*rr(ir)) * besselh(1,2,kR(ieig)*RDisk) - 1 / (RDisk*kR(ieig)^2);
                
            end
            
            if (rr(ir) >= RDisk)
                
                I1R = pi / (2 * 1i * kR(ieig)) * bj1_kR_RDisk(ieig) * bh02_kR_rr(ieig, ir);

                % I1R = pi / (2*1i*kR(ieig)) * besselj(1,kR(ieig)*RDisk) * besselh(0,2,kR(ieig)*rr(ir));
                
            end
            
            Uz_estimated(ir,ifreq) = Uz_estimated(ir,ifreq) + qDisk * RDisk * phiz(1,ieig) * phiz(1,ieig) * I1R;
                    
        end

    end
    
    for ieig = 1 : 2 * nelmt
        
       	phi = [phix(:,ieig)
            phiz(:,ieig)];
         K1_11 = kR(ieig)^2 * Ax(1:nelmt, 1:nelmt) + Gx(1:nelmt, 1:nelmt) - ww(ifreq)^2 * Mx(1:nelmt, 1:nelmt);
            K1_12 = kR(ieig) * Bxz(1:nelmt, 1:nelmt);
            K1_21 = transpose(K1_12);
            K1_22 = kR(ieig)^2 * Az(1:nelmt, 1:nelmt) + Gz(1:nelmt, 1:nelmt) - ww(ifreq)^2 * Mz(1:nelmt, 1:nelmt);
            K1 = [K1_11, K1_12; K1_21, K1_22];

    
        % K1 = kR(ieig)^2 * [Ax(1:nelmt,1:nelmt) zeros(nelmt)
        %                    zeros(nelmt)        Az(1:nelmt,1:nelmt)] ...
        %    + kR(ieig) * [zeros(nelmt)                     Bxz(1:nelmt,1:nelmt)
        %                  transpose(Bxz(1:nelmt,1:nelmt))  zeros(nelmt)] ...
        %    + [Gx(1:nelmt,1:nelmt) zeros(nelmt)
        %       zeros(nelmt)        Gz(1:nelmt,1:nelmt)] ...
        %    - ww(ifreq)^2 * [Mx(1:nelmt,1:nelmt) zeros(nelmt)
        %                     zeros(nelmt)        Mz(1:nelmt,1:nelmt)];
    
          K2_11 = 2 * kR(ieig) * Ax(1:nelmt, 1:nelmt);
            K2_12 = Bxz(1:nelmt, 1:nelmt);
            K2_21 = transpose(K2_12);
            K2_22 = 2 * kR(ieig) * Az(1:nelmt, 1:nelmt);
            K2 = [K2_11, K2_12; K2_21, K2_22];

        % K2 = 2*kR(ieig) * [Ax(1:nelmt,1:nelmt) zeros(nelmt)
        %                    zeros(nelmt)        Az(1:nelmt,1:nelmt)] ...
        %    + [zeros(nelmt)                     Bxz(1:nelmt,1:nelmt)
        %       transpose(Bxz(1:nelmt,1:nelmt))  zeros(nelmt)];
      
    	K3 = [Ax(1:nelmt,1:nelmt) zeros(nelmt)
              zeros(nelmt)        Az(1:nelmt,1:nelmt)];
      
        Q1 = [K1                K2*phi
              transpose(phi)*K2 transpose(phi)*K3*phi-1];
          
        invQ1 = inv(Q1);

        for itheta = 1 : Const
            
            if itheta ~= nelmt-3

                dmu = 100 * 2 * rho .* [
                    zeros(itheta-1,1)
                    100*theta(itheta)
                    zeros(nelmt-4-itheta,1)
                    zeros(4,1)];
    
                [dAx, dAy, dAz, dGx, dGy, dGz, dMx, dMy, dMz, dBxz] = dDynamicStiffness1(nelmt, la, mu, rho, sdamp, h, etype, ww(ifreq), itheta, dmu, 0);

            end
    
            if itheta == nelmt-3

                [dAx, dAy, dAz, dGx, dGy, dGz, dMx, dMy, dMz, dBxz] = dDynamicStiffness1(nelmt, la, mu, rho, sdamp, h, etype, ww(ifreq), itheta, dmu, 1);

            end
            
              K1_11 = kR(ieig)^2 * dAx(1:nelmt, 1:nelmt) + dGx(1:nelmt, 1:nelmt) - ww(ifreq)^2 * dMx(1:nelmt, 1:nelmt);
                K1_12 = kR(ieig) * dBxz(1:nelmt, 1:nelmt);
                K1_21 = transpose(K1_12);
                K1_22 = kR(ieig)^2 * dAz(1:nelmt, 1:nelmt) + dGz(1:nelmt, 1:nelmt) - ww(ifreq)^2 * dMz(1:nelmt, 1:nelmt);
                K1 = [K1_11, K1_12; K1_21, K1_22];

            % K1 = kR(ieig)^2 * [dAx(1:nelmt,1:nelmt) zeros(nelmt)
            %                    zeros(nelmt)         dAz(1:nelmt,1:nelmt)] ...
            %    + kR(ieig) * [zeros(nelmt)                     dBxz(1:nelmt,1:nelmt)
            %                  transpose(dBxz(1:nelmt,1:nelmt)) zeros(nelmt)] ...
            %    + [dGx(1:nelmt,1:nelmt) zeros(nelmt)
            %       zeros(nelmt)         dGz(1:nelmt,1:nelmt)] ...
            %    - ww(ifreq)^2 * [dMx(1:nelmt,1:nelmt) zeros(nelmt)
            %                     zeros(nelmt)         dMz(1:nelmt,1:nelmt)];
            % 

             K2_11 = 2 * kR(ieig) * dAx(1:nelmt, 1:nelmt);
                K2_12 = dBxz(1:nelmt, 1:nelmt);
                K2_21 = transpose(K2_12);
                K2_22 = 2 * kR(ieig) * dAz(1:nelmt, 1:nelmt);
                K2 = [K2_11, K2_12; K2_21, K2_22];

            % K2 = 2*kR(ieig) * [dAx(1:nelmt,1:nelmt) zeros(nelmt)
            %                    zeros(nelmt)         dAz(1:nelmt,1:nelmt)] ...
            %    + [zeros(nelmt)                     dBxz(1:nelmt,1:nelmt)
            %       transpose(dBxz(1:nelmt,1:nelmt)) zeros(nelmt)];
               
            Q2 = -[K1*phi
                0.5*transpose(phi)*K2*phi];
            
            Q2 = invQ1 * Q2;
            dphiz = Q2(nelmt+1:2*nelmt);
            dk = Q2(2*nelmt+1);
            
            for ir = 1 : Nmeasured
                
                if (rr(ir) < RDisk)

                    I1R = pi / (2*1i*kR(ieig)) * besselj(0,kR(ieig)*rr(ir)) * besselh(1,2,kR(ieig)*RDisk) - 1 / (RDisk*kR(ieig)^2);

                    dI1R = (1i*pi / (2*kR(ieig)) * (RDisk * besselj(0,kR(ieig)*rr(ir)) * besselh(2,2,kR(ieig)*RDisk) ...
                        + rr(ir) * besselj(1,kR(ieig)*rr(ir)) * besselh(1,2,kR(ieig)*RDisk)) + 2 / (RDisk*kR(ieig)^3)) * dk;
                
                end
            
                if (rr(ir) >= RDisk)

                         I1R = pi / (2 * 1i * kR(ieig)) * bj1_kR_RDisk(ieig) * bh02_kR_rr(ieig, ir);

                        dI1R = 1i * pi / (2 * kR(ieig)) * (RDisk * bj2_kR_RDisk(ieig) * bh02_kR_rr(ieig, ir) ...
                            +rr(ir) * bj1_kR_RDisk(ieig) * bh12_kR_rr(ieig, ir)) * dk;

                    % 
                    % I1R = pi / (2*1i*kR(ieig)) * besselj(1,kR(ieig)*RDisk) * besselh(0,2,kR(ieig)*rr(ir));
                    % 
                    % dI1R = 1i*pi / (2*kR(ieig)) * (RDisk * besselj(2,kR(ieig)*RDisk) * besselh(0,2,kR(ieig)*rr(ir)) ...
                    %     + rr(ir) * besselj(1,kR(ieig)*RDisk) * besselh(1,2,kR(ieig)*rr(ir))) * dk;
                    % 
                end
                                        
                dUz_estimated(ir,ifreq,itheta) = dUz_estimated(ir,ifreq,itheta) + 2*qDisk*RDisk * phiz(1,ieig) * dphiz(1) * I1R ...
                    + qDisk*RDisk * phiz(1,ieig) * phiz(1,ieig) * dI1R;
                
            end
            
        end
        
    end
        
end

Uz_estimated = reshape(Uz_estimated,[],1) / facUz;

Uz_estimated(nfreq*Nmeasured+1:nfreq*Nmeasured+nelmt-4) = diff(theta);

dUz_estimated = reshape(dUz_estimated,nfreq*Nmeasured,nelmt-3) / facUz;

for i = 1 : nelmt-4

    dUz_estimated(nfreq*Nmeasured+i,i) = -1;
    dUz_estimated(nfreq*Nmeasured+i,i+1) = 1;

end



