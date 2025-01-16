function [Ax, Ay, Az, Gx, Gy, Gz, Mx, My, Mz, Bxz] = DynamicStiffness(nelmt, la, mu, rho, beta, h, etype, ww)

Ax = sparse(nelmt+1,nelmt+1);
Ay = sparse(nelmt+1,nelmt+1);
Az = sparse(nelmt+1,nelmt+1);

Gx = sparse(nelmt+1,nelmt+1);
Gy = sparse(nelmt+1,nelmt+1);
Gz = sparse(nelmt+1,nelmt+1);

Mx = sparse(nelmt+1,nelmt+1);
My = sparse(nelmt+1,nelmt+1);
Mz = sparse(nelmt+1,nelmt+1);

Bxz = sparse(nelmt+1,nelmt+1);

% Ax = zeros(nelmt+1);
% Ay = zeros(nelmt+1);
% Az = zeros(nelmt+1);
% 
% Gx = zeros(nelmt+1);
% Gy = zeros(nelmt+1);
% Gz = zeros(nelmt+1);
% 
% Mx = zeros(nelmt+1);
% My = zeros(nelmt+1);
% Mz = zeros(nelmt+1);
% 
% Bxz = zeros(nelmt+1);

for ielmt = 1 : nelmt
    
    if etype(ielmt) == 0      % regular finite element
        
    Ax(ielmt  ,ielmt  ) = Ax(ielmt  ,ielmt  ) + h(ielmt)/3 * (la(ielmt) + 2 * mu(ielmt)) * (1 + 2*1i*beta(ielmt));
    Ax(ielmt  ,ielmt+1) = Ax(ielmt  ,ielmt+1) + h(ielmt)/6 * (la(ielmt) + 2 * mu(ielmt)) * (1 + 2*1i*beta(ielmt));
    Ax(ielmt+1,ielmt  ) = Ax(ielmt+1,ielmt  ) + h(ielmt)/6 * (la(ielmt) + 2 * mu(ielmt)) * (1 + 2*1i*beta(ielmt));
    Ax(ielmt+1,ielmt+1) = Ax(ielmt+1,ielmt+1) + h(ielmt)/3 * (la(ielmt) + 2 * mu(ielmt)) * (1 + 2*1i*beta(ielmt));
    
    Ay(ielmt  ,ielmt  ) = Ay(ielmt  ,ielmt  ) + h(ielmt)/3 * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    Ay(ielmt  ,ielmt+1) = Ay(ielmt  ,ielmt+1) + h(ielmt)/6 * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    Ay(ielmt+1,ielmt  ) = Ay(ielmt+1,ielmt  ) + h(ielmt)/6 * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    Ay(ielmt+1,ielmt+1) = Ay(ielmt+1,ielmt+1) + h(ielmt)/3 * mu(ielmt) * (1 + 2*1i*beta(ielmt));

    Az(ielmt  ,ielmt  ) = Az(ielmt  ,ielmt  ) + h(ielmt)/3 * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    Az(ielmt  ,ielmt+1) = Az(ielmt  ,ielmt+1) + h(ielmt)/6 * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    Az(ielmt+1,ielmt  ) = Az(ielmt+1,ielmt  ) + h(ielmt)/6 * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    Az(ielmt+1,ielmt+1) = Az(ielmt+1,ielmt+1) + h(ielmt)/3 * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    
    Gx(ielmt  ,ielmt  ) = Gx(ielmt  ,ielmt  ) + 1/h(ielmt) * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    Gx(ielmt  ,ielmt+1) = Gx(ielmt  ,ielmt+1) - 1/h(ielmt) * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    Gx(ielmt+1,ielmt  ) = Gx(ielmt+1,ielmt  ) - 1/h(ielmt) * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    Gx(ielmt+1,ielmt+1) = Gx(ielmt+1,ielmt+1) + 1/h(ielmt) * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    
    Gy(ielmt  ,ielmt  ) = Gy(ielmt  ,ielmt  ) + 1/h(ielmt) * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    Gy(ielmt  ,ielmt+1) = Gy(ielmt  ,ielmt+1) - 1/h(ielmt) * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    Gy(ielmt+1,ielmt  ) = Gy(ielmt+1,ielmt  ) - 1/h(ielmt) * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    Gy(ielmt+1,ielmt+1) = Gy(ielmt+1,ielmt+1) + 1/h(ielmt) * mu(ielmt) * (1 + 2*1i*beta(ielmt));

    Gz(ielmt  ,ielmt  ) = Gz(ielmt  ,ielmt  ) + 1/h(ielmt) * (la(ielmt) + 2 * mu(ielmt)) * (1 + 2*1i*beta(ielmt));
    Gz(ielmt  ,ielmt+1) = Gz(ielmt  ,ielmt+1) - 1/h(ielmt) * (la(ielmt) + 2 * mu(ielmt)) * (1 + 2*1i*beta(ielmt));
    Gz(ielmt+1,ielmt  ) = Gz(ielmt+1,ielmt  ) - 1/h(ielmt) * (la(ielmt) + 2 * mu(ielmt)) * (1 + 2*1i*beta(ielmt));
    Gz(ielmt+1,ielmt+1) = Gz(ielmt+1,ielmt+1) + 1/h(ielmt) * (la(ielmt) + 2 * mu(ielmt)) * (1 + 2*1i*beta(ielmt));
    
    Mx(ielmt  ,ielmt  ) = Mx(ielmt  ,ielmt  ) + h(ielmt)/3 * rho(ielmt);
    Mx(ielmt  ,ielmt+1) = Mx(ielmt  ,ielmt+1) + h(ielmt)/6 * rho(ielmt);
    Mx(ielmt+1,ielmt  ) = Mx(ielmt+1,ielmt  ) + h(ielmt)/6 * rho(ielmt);
    Mx(ielmt+1,ielmt+1) = Mx(ielmt+1,ielmt+1) + h(ielmt)/3 * rho(ielmt);
    
    My(ielmt  ,ielmt  ) = My(ielmt  ,ielmt  ) + h(ielmt)/3 * rho(ielmt);
    My(ielmt  ,ielmt+1) = My(ielmt  ,ielmt+1) + h(ielmt)/6 * rho(ielmt);
    My(ielmt+1,ielmt  ) = My(ielmt+1,ielmt  ) + h(ielmt)/6 * rho(ielmt);
    My(ielmt+1,ielmt+1) = My(ielmt+1,ielmt+1) + h(ielmt)/3 * rho(ielmt);

    Mz(ielmt  ,ielmt  ) = Mz(ielmt  ,ielmt  ) + h(ielmt)/3 * rho(ielmt);
    Mz(ielmt  ,ielmt+1) = Mz(ielmt  ,ielmt+1) + h(ielmt)/6 * rho(ielmt);
    Mz(ielmt+1,ielmt  ) = Mz(ielmt+1,ielmt  ) + h(ielmt)/6 * rho(ielmt);
    Mz(ielmt+1,ielmt+1) = Mz(ielmt+1,ielmt+1) + h(ielmt)/3 * rho(ielmt);
    
    Bxz(ielmt  ,ielmt  ) = Bxz(ielmt  ,ielmt  ) + 0.5 * (la(ielmt) - mu(ielmt)) * (1 + 2*1i*beta(ielmt));
    Bxz(ielmt  ,ielmt+1) = Bxz(ielmt  ,ielmt+1) - 0.5 * (la(ielmt) + mu(ielmt)) * (1 + 2*1i*beta(ielmt));
    Bxz(ielmt+1,ielmt  ) = Bxz(ielmt+1,ielmt  ) + 0.5 * (la(ielmt) + mu(ielmt)) * (1 + 2*1i*beta(ielmt));
    Bxz(ielmt+1,ielmt+1) = Bxz(ielmt+1,ielmt+1) - 0.5 * (la(ielmt) - mu(ielmt)) * (1 + 2*1i*beta(ielmt));
    
    end

    if etype(ielmt) ~= 0      % complex-length finite element or perfectly matched discrete layer
        
        if etype(ielmt) == 21
            
            nu = la(ielmt) / 2 / (la(ielmt) + mu(ielmt));
            CP = sqrt((la(ielmt) + 2 * mu(ielmt)) / rho(ielmt));
            CS = sqrt(mu(ielmt) / rho(ielmt));
            CR = (0.862 + 1.14*nu) / (1 + nu) * CS;
            aP = CP / CR;
            aS = CS / CR;
            h(ielmt) = 2 * CP / ww / sqrt(aP^2 - 1);
            
        end
        
        if etype(ielmt) == 22
            
            nu = la(ielmt) / 2 / (la(ielmt) + mu(ielmt));
            CP = sqrt((la(ielmt) + 2 * mu(ielmt)) / rho(ielmt));
            CS = sqrt(mu(ielmt) / rho(ielmt));
            CR = (0.862 + 1.14*nu) / (1 + nu) * CS;
            aP = CP / CR;
            aS = CS / CR;
            h(ielmt) = 2 * CS / ww / sqrt(aS^2 - 1);
            
        end
        
        if etype(ielmt) == 31
            
            nu = la(ielmt) / 2 / (la(ielmt) + mu(ielmt));
            CP = sqrt((la(ielmt) + 2 * mu(ielmt)) / rho(ielmt));
            CS = sqrt(mu(ielmt) / rho(ielmt));
            CR = (0.862 + 1.14*nu) / (1 + nu) * CS;
            aP = CP / CR;
            aS = CS / CR;
            h(ielmt) = -2*1i * CP / ww;
            
        end
        
        if etype(ielmt) == 32
            
            nu = la(ielmt) / 2 / (la(ielmt) + mu(ielmt));
            CP = sqrt((la(ielmt) + 2 * mu(ielmt)) / rho(ielmt));
            CS = sqrt(mu(ielmt) / rho(ielmt));
            CR = (0.862 + 1.14*nu) / (1 + nu) * CS;
            aP = CP / CR;
            aS = CS / CR;
            h(ielmt) = -2*1i * CS / ww;
            
        end
        
	Ax(ielmt  ,ielmt  ) = Ax(ielmt  ,ielmt  ) + h(ielmt)/4 * (la(ielmt) + 2 * mu(ielmt)) * (1 + 2*1i*beta(ielmt));
    Ax(ielmt  ,ielmt+1) = Ax(ielmt  ,ielmt+1) + h(ielmt)/4 * (la(ielmt) + 2 * mu(ielmt)) * (1 + 2*1i*beta(ielmt));
    Ax(ielmt+1,ielmt  ) = Ax(ielmt+1,ielmt  ) + h(ielmt)/4 * (la(ielmt) + 2 * mu(ielmt)) * (1 + 2*1i*beta(ielmt));
    Ax(ielmt+1,ielmt+1) = Ax(ielmt+1,ielmt+1) + h(ielmt)/4 * (la(ielmt) + 2 * mu(ielmt)) * (1 + 2*1i*beta(ielmt));
    
    Ay(ielmt  ,ielmt  ) = Ay(ielmt  ,ielmt  ) + h(ielmt)/4 * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    Ay(ielmt  ,ielmt+1) = Ay(ielmt  ,ielmt+1) + h(ielmt)/4 * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    Ay(ielmt+1,ielmt  ) = Ay(ielmt+1,ielmt  ) + h(ielmt)/4 * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    Ay(ielmt+1,ielmt+1) = Ay(ielmt+1,ielmt+1) + h(ielmt)/4 * mu(ielmt) * (1 + 2*1i*beta(ielmt));

    Az(ielmt  ,ielmt  ) = Az(ielmt  ,ielmt  ) + h(ielmt)/4 * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    Az(ielmt  ,ielmt+1) = Az(ielmt  ,ielmt+1) + h(ielmt)/4 * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    Az(ielmt+1,ielmt  ) = Az(ielmt+1,ielmt  ) + h(ielmt)/4 * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    Az(ielmt+1,ielmt+1) = Az(ielmt+1,ielmt+1) + h(ielmt)/4 * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    
    Gx(ielmt  ,ielmt  ) = Gx(ielmt  ,ielmt  ) + 1/h(ielmt) * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    Gx(ielmt  ,ielmt+1) = Gx(ielmt  ,ielmt+1) - 1/h(ielmt) * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    Gx(ielmt+1,ielmt  ) = Gx(ielmt+1,ielmt  ) - 1/h(ielmt) * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    Gx(ielmt+1,ielmt+1) = Gx(ielmt+1,ielmt+1) + 1/h(ielmt) * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    
    Gy(ielmt  ,ielmt  ) = Gy(ielmt  ,ielmt  ) + 1/h(ielmt) * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    Gy(ielmt  ,ielmt+1) = Gy(ielmt  ,ielmt+1) - 1/h(ielmt) * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    Gy(ielmt+1,ielmt  ) = Gy(ielmt+1,ielmt  ) - 1/h(ielmt) * mu(ielmt) * (1 + 2*1i*beta(ielmt));
    Gy(ielmt+1,ielmt+1) = Gy(ielmt+1,ielmt+1) + 1/h(ielmt) * mu(ielmt) * (1 + 2*1i*beta(ielmt));

    Gz(ielmt  ,ielmt  ) = Gz(ielmt  ,ielmt  ) + 1/h(ielmt) * (la(ielmt) + 2 * mu(ielmt)) * (1 + 2*1i*beta(ielmt));
    Gz(ielmt  ,ielmt+1) = Gz(ielmt  ,ielmt+1) - 1/h(ielmt) * (la(ielmt) + 2 * mu(ielmt)) * (1 + 2*1i*beta(ielmt));
    Gz(ielmt+1,ielmt  ) = Gz(ielmt+1,ielmt  ) - 1/h(ielmt) * (la(ielmt) + 2 * mu(ielmt)) * (1 + 2*1i*beta(ielmt));
    Gz(ielmt+1,ielmt+1) = Gz(ielmt+1,ielmt+1) + 1/h(ielmt) * (la(ielmt) + 2 * mu(ielmt)) * (1 + 2*1i*beta(ielmt));
    
    Mx(ielmt  ,ielmt  ) = Mx(ielmt  ,ielmt  ) + h(ielmt)/4 * rho(ielmt);
    Mx(ielmt  ,ielmt+1) = Mx(ielmt  ,ielmt+1) + h(ielmt)/4 * rho(ielmt);
    Mx(ielmt+1,ielmt  ) = Mx(ielmt+1,ielmt  ) + h(ielmt)/4 * rho(ielmt);
    Mx(ielmt+1,ielmt+1) = Mx(ielmt+1,ielmt+1) + h(ielmt)/4 * rho(ielmt);
    
    My(ielmt  ,ielmt  ) = My(ielmt  ,ielmt  ) + h(ielmt)/4 * rho(ielmt);
    My(ielmt  ,ielmt+1) = My(ielmt  ,ielmt+1) + h(ielmt)/4 * rho(ielmt);
    My(ielmt+1,ielmt  ) = My(ielmt+1,ielmt  ) + h(ielmt)/4 * rho(ielmt);
    My(ielmt+1,ielmt+1) = My(ielmt+1,ielmt+1) + h(ielmt)/4 * rho(ielmt);

    Mz(ielmt  ,ielmt  ) = Mz(ielmt  ,ielmt  ) + h(ielmt)/4 * rho(ielmt);
    Mz(ielmt  ,ielmt+1) = Mz(ielmt  ,ielmt+1) + h(ielmt)/4 * rho(ielmt);
    Mz(ielmt+1,ielmt  ) = Mz(ielmt+1,ielmt  ) + h(ielmt)/4 * rho(ielmt);
    Mz(ielmt+1,ielmt+1) = Mz(ielmt+1,ielmt+1) + h(ielmt)/4 * rho(ielmt);
    
    Bxz(ielmt  ,ielmt  ) = Bxz(ielmt  ,ielmt  ) + 0.5 * (la(ielmt) - mu(ielmt)) * (1 + 2*1i*beta(ielmt));
    Bxz(ielmt  ,ielmt+1) = Bxz(ielmt  ,ielmt+1) - 0.5 * (la(ielmt) + mu(ielmt)) * (1 + 2*1i*beta(ielmt));
    Bxz(ielmt+1,ielmt  ) = Bxz(ielmt+1,ielmt  ) + 0.5 * (la(ielmt) + mu(ielmt)) * (1 + 2*1i*beta(ielmt));
    Bxz(ielmt+1,ielmt+1) = Bxz(ielmt+1,ielmt+1) - 0.5 * (la(ielmt) - mu(ielmt)) * (1 + 2*1i*beta(ielmt));
    
    end
    
end

end