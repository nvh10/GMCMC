function [dAx, dAy, dAz, dGx, dGy, dGz, dMx, dMy, dMz, dBxz] = dDynamicStiffness1(nelmt, la, mu, rho, beta, h, etype, ww, itheta, dmu, ihalf)

dAx = sparse(nelmt+1,nelmt+1);
dAy = sparse(nelmt+1,nelmt+1);
dAz = sparse(nelmt+1,nelmt+1);

dGx = sparse(nelmt+1,nelmt+1);
dGy = sparse(nelmt+1,nelmt+1);
dGz = sparse(nelmt+1,nelmt+1);

dMx = sparse(nelmt+1,nelmt+1);
dMy = sparse(nelmt+1,nelmt+1);
dMz = sparse(nelmt+1,nelmt+1);

dBxz = sparse(nelmt+1,nelmt+1);


% dAx = zeros(nelmt+1);
% dAy = zeros(nelmt+1);
% dAz = zeros(nelmt+1);
% 
% dGx = zeros(nelmt+1);
% dGy = zeros(nelmt+1);
% dGz = zeros(nelmt+1);
% 
% dMx = zeros(nelmt+1);
% dMy = zeros(nelmt+1);
% dMz = zeros(nelmt+1);
% 
% dBxz = zeros(nelmt+1);

nu = la ./ (2*(la + mu));

%for ielmt = 1 : nelmt
    
    if ihalf == 0
        
    ielmt = itheta;
    
    if etype(ielmt) == 0     % regular finite element
        
    dAx(ielmt  ,ielmt  ) = dAx(ielmt  ,ielmt  ) + h(ielmt)/3 * dmu(ielmt) * 2*(1-nu(ielmt))/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dAx(ielmt  ,ielmt+1) = dAx(ielmt  ,ielmt+1) + h(ielmt)/6 * dmu(ielmt) * 2*(1-nu(ielmt))/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dAx(ielmt+1,ielmt  ) = dAx(ielmt+1,ielmt  ) + h(ielmt)/6 * dmu(ielmt) * 2*(1-nu(ielmt))/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dAx(ielmt+1,ielmt+1) = dAx(ielmt+1,ielmt+1) + h(ielmt)/3 * dmu(ielmt) * 2*(1-nu(ielmt))/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    
    dAy(ielmt  ,ielmt  ) = dAy(ielmt  ,ielmt  ) + h(ielmt)/3 * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    dAy(ielmt  ,ielmt+1) = dAy(ielmt  ,ielmt+1) + h(ielmt)/6 * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    dAy(ielmt+1,ielmt  ) = dAy(ielmt+1,ielmt  ) + h(ielmt)/6 * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    dAy(ielmt+1,ielmt+1) = dAy(ielmt+1,ielmt+1) + h(ielmt)/3 * dmu(ielmt) * (1 + 2*1i*beta(ielmt));

    dAz(ielmt  ,ielmt  ) = dAz(ielmt  ,ielmt  ) + h(ielmt)/3 * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    dAz(ielmt  ,ielmt+1) = dAz(ielmt  ,ielmt+1) + h(ielmt)/6 * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    dAz(ielmt+1,ielmt  ) = dAz(ielmt+1,ielmt  ) + h(ielmt)/6 * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    dAz(ielmt+1,ielmt+1) = dAz(ielmt+1,ielmt+1) + h(ielmt)/3 * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    
    dGx(ielmt  ,ielmt  ) = dGx(ielmt  ,ielmt  ) + 1/h(ielmt) * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    dGx(ielmt  ,ielmt+1) = dGx(ielmt  ,ielmt+1) - 1/h(ielmt) * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    dGx(ielmt+1,ielmt  ) = dGx(ielmt+1,ielmt  ) - 1/h(ielmt) * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    dGx(ielmt+1,ielmt+1) = dGx(ielmt+1,ielmt+1) + 1/h(ielmt) * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    
    dGy(ielmt  ,ielmt  ) = dGy(ielmt  ,ielmt  ) + 1/h(ielmt) * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    dGy(ielmt  ,ielmt+1) = dGy(ielmt  ,ielmt+1) - 1/h(ielmt) * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    dGy(ielmt+1,ielmt  ) = dGy(ielmt+1,ielmt  ) - 1/h(ielmt) * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    dGy(ielmt+1,ielmt+1) = dGy(ielmt+1,ielmt+1) + 1/h(ielmt) * dmu(ielmt) * (1 + 2*1i*beta(ielmt));

    dGz(ielmt  ,ielmt  ) = dGz(ielmt  ,ielmt  ) + 1/h(ielmt) * dmu(ielmt) * 2*(1-nu(ielmt))/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dGz(ielmt  ,ielmt+1) = dGz(ielmt  ,ielmt+1) - 1/h(ielmt) * dmu(ielmt) * 2*(1-nu(ielmt))/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dGz(ielmt+1,ielmt  ) = dGz(ielmt+1,ielmt  ) - 1/h(ielmt) * dmu(ielmt) * 2*(1-nu(ielmt))/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dGz(ielmt+1,ielmt+1) = dGz(ielmt+1,ielmt+1) + 1/h(ielmt) * dmu(ielmt) * 2*(1-nu(ielmt))/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    
%     dMx(ielmt  ,ielmt  ) = dMx(ielmt  ,ielmt  ) + h(ielmt)/3 * rho(ielmt);
%     dMx(ielmt  ,ielmt+1) = dMx(ielmt  ,ielmt+1) + h(ielmt)/6 * rho(ielmt);
%     dMx(ielmt+1,ielmt  ) = dMx(ielmt+1,ielmt  ) + h(ielmt)/6 * rho(ielmt);
%     dMx(ielmt+1,ielmt+1) = dMx(ielmt+1,ielmt+1) + h(ielmt)/3 * rho(ielmt);
%     
%     dMy(ielmt  ,ielmt  ) = dMy(ielmt  ,ielmt  ) + h(ielmt)/3 * rho(ielmt);
%     dMy(ielmt  ,ielmt+1) = dMy(ielmt  ,ielmt+1) + h(ielmt)/6 * rho(ielmt);
%     dMy(ielmt+1,ielmt  ) = dMy(ielmt+1,ielmt  ) + h(ielmt)/6 * rho(ielmt);
%     dMy(ielmt+1,ielmt+1) = dMy(ielmt+1,ielmt+1) + h(ielmt)/3 * rho(ielmt);
% 
%     dMz(ielmt  ,ielmt  ) = dMz(ielmt  ,ielmt  ) + h(ielmt)/3 * rho(ielmt);
%     dMz(ielmt  ,ielmt+1) = dMz(ielmt  ,ielmt+1) + h(ielmt)/6 * rho(ielmt);
%     dMz(ielmt+1,ielmt  ) = dMz(ielmt+1,ielmt  ) + h(ielmt)/6 * rho(ielmt);
%     dMz(ielmt+1,ielmt+1) = dMz(ielmt+1,ielmt+1) + h(ielmt)/3 * rho(ielmt);
%     
    dBxz(ielmt  ,ielmt  ) = dBxz(ielmt  ,ielmt  ) + 0.5 * (4*nu(ielmt)-1)/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dBxz(ielmt  ,ielmt+1) = dBxz(ielmt  ,ielmt+1) - 0.5 * 1/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dBxz(ielmt+1,ielmt  ) = dBxz(ielmt+1,ielmt  ) + 0.5 * 1/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dBxz(ielmt+1,ielmt+1) = dBxz(ielmt+1,ielmt+1) - 0.5 * (4*nu(ielmt)-1)/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    
    end

    if etype(ielmt) == 1      % complex-length finite element or perfectly matched discrete layer
        
	dAx(ielmt  ,ielmt  ) = dAx(ielmt  ,ielmt  ) + h(ielmt)/4 * dmu(ielmt) * 2*(1-nu(ielmt))/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dAx(ielmt  ,ielmt+1) = dAx(ielmt  ,ielmt+1) + h(ielmt)/4 * dmu(ielmt) * 2*(1-nu(ielmt))/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dAx(ielmt+1,ielmt  ) = dAx(ielmt+1,ielmt  ) + h(ielmt)/4 * dmu(ielmt) * 2*(1-nu(ielmt))/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dAx(ielmt+1,ielmt+1) = dAx(ielmt+1,ielmt+1) + h(ielmt)/4 * dmu(ielmt) * 2*(1-nu(ielmt))/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    
    dAy(ielmt  ,ielmt  ) = dAy(ielmt  ,ielmt  ) + h(ielmt)/4 * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    dAy(ielmt  ,ielmt+1) = dAy(ielmt  ,ielmt+1) + h(ielmt)/4 * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    dAy(ielmt+1,ielmt  ) = dAy(ielmt+1,ielmt  ) + h(ielmt)/4 * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    dAy(ielmt+1,ielmt+1) = dAy(ielmt+1,ielmt+1) + h(ielmt)/4 * dmu(ielmt) * (1 + 2*1i*beta(ielmt));

    dAz(ielmt  ,ielmt  ) = dAz(ielmt  ,ielmt  ) + h(ielmt)/4 * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    dAz(ielmt  ,ielmt+1) = dAz(ielmt  ,ielmt+1) + h(ielmt)/4 * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    dAz(ielmt+1,ielmt  ) = dAz(ielmt+1,ielmt  ) + h(ielmt)/4 * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    dAz(ielmt+1,ielmt+1) = dAz(ielmt+1,ielmt+1) + h(ielmt)/4 * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    
    dGx(ielmt  ,ielmt  ) = dGx(ielmt  ,ielmt  ) + 1/h(ielmt) * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    dGx(ielmt  ,ielmt+1) = dGx(ielmt  ,ielmt+1) - 1/h(ielmt) * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    dGx(ielmt+1,ielmt  ) = dGx(ielmt+1,ielmt  ) - 1/h(ielmt) * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    dGx(ielmt+1,ielmt+1) = dGx(ielmt+1,ielmt+1) + 1/h(ielmt) * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    
    dGy(ielmt  ,ielmt  ) = dGy(ielmt  ,ielmt  ) + 1/h(ielmt) * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    dGy(ielmt  ,ielmt+1) = dGy(ielmt  ,ielmt+1) - 1/h(ielmt) * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    dGy(ielmt+1,ielmt  ) = dGy(ielmt+1,ielmt  ) - 1/h(ielmt) * dmu(ielmt) * (1 + 2*1i*beta(ielmt));
    dGy(ielmt+1,ielmt+1) = dGy(ielmt+1,ielmt+1) + 1/h(ielmt) * dmu(ielmt) * (1 + 2*1i*beta(ielmt));

    dGz(ielmt  ,ielmt  ) = dGz(ielmt  ,ielmt  ) + 1/h(ielmt) * dmu(ielmt) * 2*(1-nu(ielmt))/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dGz(ielmt  ,ielmt+1) = dGz(ielmt  ,ielmt+1) - 1/h(ielmt) * dmu(ielmt) * 2*(1-nu(ielmt))/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dGz(ielmt+1,ielmt  ) = dGz(ielmt+1,ielmt  ) - 1/h(ielmt) * dmu(ielmt) * 2*(1-nu(ielmt))/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dGz(ielmt+1,ielmt+1) = dGz(ielmt+1,ielmt+1) + 1/h(ielmt) * dmu(ielmt) * 2*(1-nu(ielmt))/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    
%     dMx(ielmt  ,ielmt  ) = dMx(ielmt  ,ielmt  ) + h(ielmt)/4 * rho(ielmt);
%     dMx(ielmt  ,ielmt+1) = dMx(ielmt  ,ielmt+1) + h(ielmt)/4 * rho(ielmt);
%     dMx(ielmt+1,ielmt  ) = dMx(ielmt+1,ielmt  ) + h(ielmt)/4 * rho(ielmt);
%     dMx(ielmt+1,ielmt+1) = dMx(ielmt+1,ielmt+1) + h(ielmt)/4 * rho(ielmt);
%     
%     dMy(ielmt  ,ielmt  ) = dMy(ielmt  ,ielmt  ) + h(ielmt)/4 * rho(ielmt);
%     dMy(ielmt  ,ielmt+1) = dMy(ielmt  ,ielmt+1) + h(ielmt)/4 * rho(ielmt);
%     dMy(ielmt+1,ielmt  ) = dMy(ielmt+1,ielmt  ) + h(ielmt)/4 * rho(ielmt);
%     dMy(ielmt+1,ielmt+1) = dMy(ielmt+1,ielmt+1) + h(ielmt)/4 * rho(ielmt);
% 
%     dMz(ielmt  ,ielmt  ) = dMz(ielmt  ,ielmt  ) + h(ielmt)/4 * rho(ielmt);
%     dMz(ielmt  ,ielmt+1) = dMz(ielmt  ,ielmt+1) + h(ielmt)/4 * rho(ielmt);
%     dMz(ielmt+1,ielmt  ) = dMz(ielmt+1,ielmt  ) + h(ielmt)/4 * rho(ielmt);
%     dMz(ielmt+1,ielmt+1) = dMz(ielmt+1,ielmt+1) + h(ielmt)/4 * rho(ielmt);
    
    dBxz(ielmt  ,ielmt  ) = dBxz(ielmt  ,ielmt  ) + 0.5 * (4*nu(ielmt)-1)/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dBxz(ielmt  ,ielmt+1) = dBxz(ielmt  ,ielmt+1) - 0.5 * 1/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dBxz(ielmt+1,ielmt  ) = dBxz(ielmt+1,ielmt  ) + 0.5 * 1/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dBxz(ielmt+1,ielmt+1) = dBxz(ielmt+1,ielmt+1) - 0.5 * (4*nu(ielmt)-1)/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    
    end
    
    end
    
    if ihalf == 1      % complex-length finite element or perfectly matched discrete layer
        
	for ielmt = nelmt-3 : nelmt
        
        if etype(ielmt) == 21
            
            CP = sqrt((la(ielmt) + 2 * mu(ielmt)) / rho(ielmt));
            CS = sqrt(mu(ielmt) / rho(ielmt));
            CR = (0.862 + 1.14*nu(ielmt)) / (1 + nu(ielmt)) * CS;
            aP = CP / CR;
            aS = CS / CR;
            bpe = 2 / ww / sqrt(aP^2 - 1);
            bPS = sqrt(2*(1-nu(ielmt))/(1-2*nu(ielmt)));
            
        end
        
        if etype(ielmt) == 22
            
            CP = sqrt((la(ielmt) + 2 * mu(ielmt)) / rho(ielmt));
            CS = sqrt(mu(ielmt) / rho(ielmt));
            CR = (0.862 + 1.14*nu(ielmt)) / (1 + nu(ielmt)) * CS;
            aP = CP / CR;
            aS = CS / CR;
            bpe = 2 / ww / sqrt(aS^2 - 1);
            bPS = 1;
            
        end
        
        if etype(ielmt) == 31
            
            CP = sqrt((la(ielmt) + 2 * mu(ielmt)) / rho(ielmt));
            CS = sqrt(mu(ielmt) / rho(ielmt));
            CR = (0.862 + 1.14*nu(ielmt)) / (1 + nu(ielmt)) * CS;
            aP = CP / CR;
            aS = CS / CR;
            bpe = -2*1i / ww;
            bPS = sqrt(2*(1-nu(ielmt))/(1-2*nu(ielmt)));
            
        end
        
        if etype(ielmt) == 32
            
            CP = sqrt((la(ielmt) + 2 * mu(ielmt)) / rho(ielmt));
            CS = sqrt(mu(ielmt) / rho(ielmt));
            CR = (0.862 + 1.14*nu(ielmt)) / (1 + nu(ielmt)) * CS;
            aP = CP / CR;
            aS = CS / CR;
            bpe = -2*1i / ww;
            bPS = 1;
            
        end
        
	dAx(ielmt  ,ielmt  ) = dAx(ielmt  ,ielmt  ) + 100 * 3*bpe*bPS/4 * rho(ielmt) * CS^2 * 2*(1-nu(ielmt))/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dAx(ielmt  ,ielmt+1) = dAx(ielmt  ,ielmt+1) + 100 * 3*bpe*bPS/4 * rho(ielmt) * CS^2 * 2*(1-nu(ielmt))/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dAx(ielmt+1,ielmt  ) = dAx(ielmt+1,ielmt  ) + 100 * 3*bpe*bPS/4 * rho(ielmt) * CS^2 * 2*(1-nu(ielmt))/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dAx(ielmt+1,ielmt+1) = dAx(ielmt+1,ielmt+1) + 100 * 3*bpe*bPS/4 * rho(ielmt) * CS^2 * 2*(1-nu(ielmt))/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    
    dAy(ielmt  ,ielmt  ) = dAy(ielmt  ,ielmt  ) + 100 * 3*bpe*bPS/4 * rho(ielmt) * CS^2 * (1 + 2*1i*beta(ielmt));
    dAy(ielmt  ,ielmt+1) = dAy(ielmt  ,ielmt+1) + 100 * 3*bpe*bPS/4 * rho(ielmt) * CS^2 * (1 + 2*1i*beta(ielmt));
    dAy(ielmt+1,ielmt  ) = dAy(ielmt+1,ielmt  ) + 100 * 3*bpe*bPS/4 * rho(ielmt) * CS^2 * (1 + 2*1i*beta(ielmt));
    dAy(ielmt+1,ielmt+1) = dAy(ielmt+1,ielmt+1) + 100 * 3*bpe*bPS/4 * rho(ielmt) * CS^2 * (1 + 2*1i*beta(ielmt));

    dAz(ielmt  ,ielmt  ) = dAz(ielmt  ,ielmt  ) + 100 * 3*bpe*bPS/4 * rho(ielmt) * CS^2 * (1 + 2*1i*beta(ielmt));
    dAz(ielmt  ,ielmt+1) = dAz(ielmt  ,ielmt+1) + 100 * 3*bpe*bPS/4 * rho(ielmt) * CS^2 * (1 + 2*1i*beta(ielmt));
    dAz(ielmt+1,ielmt  ) = dAz(ielmt+1,ielmt  ) + 100 * 3*bpe*bPS/4 * rho(ielmt) * CS^2 * (1 + 2*1i*beta(ielmt));
    dAz(ielmt+1,ielmt+1) = dAz(ielmt+1,ielmt+1) + 100 * 3*bpe*bPS/4 * rho(ielmt) * CS^2 * (1 + 2*1i*beta(ielmt));
    
    dGx(ielmt  ,ielmt  ) = dGx(ielmt  ,ielmt  ) + 100 * rho(ielmt) / (bpe*bPS) * (1 + 2*1i*beta(ielmt));
    dGx(ielmt  ,ielmt+1) = dGx(ielmt  ,ielmt+1) - 100 * rho(ielmt) / (bpe*bPS) * (1 + 2*1i*beta(ielmt));
    dGx(ielmt+1,ielmt  ) = dGx(ielmt+1,ielmt  ) - 100 * rho(ielmt) / (bpe*bPS) * (1 + 2*1i*beta(ielmt));
    dGx(ielmt+1,ielmt+1) = dGx(ielmt+1,ielmt+1) + 100 * rho(ielmt) / (bpe*bPS) * (1 + 2*1i*beta(ielmt));
    
    dGy(ielmt  ,ielmt  ) = dGy(ielmt  ,ielmt  ) + 100 * rho(ielmt) / (bpe*bPS) * (1 + 2*1i*beta(ielmt));
    dGy(ielmt  ,ielmt+1) = dGy(ielmt  ,ielmt+1) - 100 * rho(ielmt) / (bpe*bPS) * (1 + 2*1i*beta(ielmt));
    dGy(ielmt+1,ielmt  ) = dGy(ielmt+1,ielmt  ) - 100 * rho(ielmt) / (bpe*bPS) * (1 + 2*1i*beta(ielmt));
    dGy(ielmt+1,ielmt+1) = dGy(ielmt+1,ielmt+1) + 100 * rho(ielmt) / (bpe*bPS) * (1 + 2*1i*beta(ielmt));

    dGz(ielmt  ,ielmt  ) = dGz(ielmt  ,ielmt  ) + 100 * rho(ielmt) / (bpe*bPS) * 2*(1-nu(ielmt))/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dGz(ielmt  ,ielmt+1) = dGz(ielmt  ,ielmt+1) - 100 * rho(ielmt) / (bpe*bPS) * 2*(1-nu(ielmt))/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dGz(ielmt+1,ielmt  ) = dGz(ielmt+1,ielmt  ) - 100 * rho(ielmt) / (bpe*bPS) * 2*(1-nu(ielmt))/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dGz(ielmt+1,ielmt+1) = dGz(ielmt+1,ielmt+1) + 100 * rho(ielmt) / (bpe*bPS) * 2*(1-nu(ielmt))/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    
    dMx(ielmt  ,ielmt  ) = dMx(ielmt  ,ielmt  ) + 100 * bpe*bPS/4 * rho(ielmt);
    dMx(ielmt  ,ielmt+1) = dMx(ielmt  ,ielmt+1) + 100 * bpe*bPS/4 * rho(ielmt);
    dMx(ielmt+1,ielmt  ) = dMx(ielmt+1,ielmt  ) + 100 * bpe*bPS/4 * rho(ielmt);
    dMx(ielmt+1,ielmt+1) = dMx(ielmt+1,ielmt+1) + 100 * bpe*bPS/4 * rho(ielmt);
    
    dMy(ielmt  ,ielmt  ) = dMy(ielmt  ,ielmt  ) + 100 * bpe*bPS/4 * rho(ielmt);
    dMy(ielmt  ,ielmt+1) = dMy(ielmt  ,ielmt+1) + 100 * bpe*bPS/4 * rho(ielmt);
    dMy(ielmt+1,ielmt  ) = dMy(ielmt+1,ielmt  ) + 100 * bpe*bPS/4 * rho(ielmt);
    dMy(ielmt+1,ielmt+1) = dMy(ielmt+1,ielmt+1) + 100 * bpe*bPS/4 * rho(ielmt);

    dMz(ielmt  ,ielmt  ) = dMz(ielmt  ,ielmt  ) + 100 * bpe*bPS/4 * rho(ielmt);
    dMz(ielmt  ,ielmt+1) = dMz(ielmt  ,ielmt+1) + 100 * bpe*bPS/4 * rho(ielmt);
    dMz(ielmt+1,ielmt  ) = dMz(ielmt+1,ielmt  ) + 100 * bpe*bPS/4 * rho(ielmt);
    dMz(ielmt+1,ielmt+1) = dMz(ielmt+1,ielmt+1) + 100 * bpe*bPS/4 * rho(ielmt);
    
    dBxz(ielmt  ,ielmt  ) = dBxz(ielmt  ,ielmt  ) + 100 * rho(ielmt)*CS * (4*nu(ielmt)-1)/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dBxz(ielmt  ,ielmt+1) = dBxz(ielmt  ,ielmt+1) - 100 * rho(ielmt)*CS * 1/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dBxz(ielmt+1,ielmt  ) = dBxz(ielmt+1,ielmt  ) + 100 * rho(ielmt)*CS * 1/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    dBxz(ielmt+1,ielmt+1) = dBxz(ielmt+1,ielmt+1) - 100 * rho(ielmt)*CS * (4*nu(ielmt)-1)/(1-2*nu(ielmt)) * (1 + 2*1i*beta(ielmt));
    
    end
    
    end
    
%end