clear

halfspace1();
% halfspace2();
% halfspace2();
% pavement();
% realsite();
mu = rho .* Cs.^2;
la = 2 * mu .* nu ./ (1 - 2 .* nu);

Uz_measured = zeros(nfreq, Nmeasured);

[Uz_measured] = SurfaceDisplacementFreqDm(nelmt, la, mu, rho, sdamp, h, etype, ww, rr, RDisk, qDisk);
facUz = max(max(abs(Uz_measured)));

Uz_measured = reshape(Uz_measured.', [], 1);
Uz_measured(nfreq*Nmeasured+1:nfreq*Nmeasured+nelmt-4) = 0;
Uz_measured = Uz_measured / facUz;

theta = zeros(nelmt-3, nSample);
theta(:, 1) = theta0;

lnprop01 = ones(nSample, 1);
lnprop10 = ones(nSample, 1);
lnpost0 = ones(nSample, 1);
lnpost1 = ones(nSample, 1);
gamma = zeros(nSample, 1);

Ca_D(1:nfreq*Nmeasured*2, 1:nfreq*Nmeasured*2) = eye(nfreq*Nmeasured*2) * covUz(1, 1);
Ca_D(nfreq*Nmeasured*2+1:nfreq*Nmeasured*2+nelmt-4, nfreq*Nmeasured*2+1:nfreq*Nmeasured*2+nelmt-4) = eye(nelmt-4) * covUz(end, end);
inv_Ca_D = 1 / covUz * eye(nfreq*Nmeasured*2+nelmt-4);

for ii = 1:nelmt - 4
    inv_Ca_D(nfreq*Nmeasured*2+ii, nfreq*Nmeasured*2+ii) = 1 / covdtheta;
end
Ca_D = inv(inv_Ca_D);

Ca_M = zeros(nelmt-3);
inv_Ca_M = zeros(nelmt-3);

for ie = 1:nelmt - 3
    Ca_M(ie, ie) = (0.5 * covtheta0);
end

inv_Ca_M = inv(Ca_M);
inv_covtheta0 = inv_Ca_M;
covtheta0 = Ca_M;

da_obs = [Uz_measured(1:nfreq * Nmeasured); ...
    conj(Uz_measured(1:nfreq * Nmeasured)); ...
    Uz_measured(nfreq * Nmeasured + 1:end)];

[Uz_estimated(:, 1), dUz_estimated(:, :, 1)] = ForwardModel(theta(:, 1), nelmt, nu, rho, sdamp, h, etype, ww, rr, RDisk, qDisk, CR, facUz);
beginSample = 2;

a = 0;
for iSample = 2:nSample
    tic

    Ga = [dUz_estimated(1:nfreq * Nmeasured, :, iSample - 1); ...
        conj(dUz_estimated(1:nfreq * Nmeasured, :, iSample - 1)); ...
        dUz_estimated(nfreq * Nmeasured + 1:end, :, iSample - 1)];
    d = Uz_estimated(1:nfreq*Nmeasured, iSample-1);
    cd = conj(Uz_estimated(1:nfreq * Nmeasured, iSample - 1));
    dre = Uz_estimated(nfreq*Nmeasured+1:end, iSample-1);
    da_cal = [d; cd; dre];
    % g = real(Ga'*inv_Ca_D*(da_cal - da_obs)) + inv_Ca_M * (theta(:, iSample - 1) - theta0);
    %
    % H = real(Ga'*inv_Ca_D*Ga) + inv_Ca_M;
    %
    % dtheta = -H \ g;

    K = Ca_D + Ga * Ca_M * Ga';
    Fa = Ca_M * Ga' * inv(K);
    IFaGa = (eye(nelmt - 3) - Fa * Ga);
    COV = real(2*IFaGa*Ca_M*IFaGa'+Fa*Ca_D*Fa');
    COV = beta * 0.5 * (COV + COV)';
    COV2 = COV;

    B = chol(COV, 'lower') * normrnd(0, 1, nelmt-3, 1);

    %theta(:, iSample) =  theta(:, iSample-1) + alpha * dtheta + B;
    theta(:, iSample) = alpha * theta0 + (1 - alpha) * theta(:, iSample-1) - alpha * real(Fa*(da_cal - da_obs - Ga * (theta(:, iSample - 1) - theta0))) + B;


    lnpost0(iSample, 1) = -(0.5 * real((da_cal-da_obs)' * inv_Ca_D * (da_cal - da_obs)) ...
        +0.5 * (theta(:, iSample - 1) - theta0)' * inv_covtheta0 * (theta(:, iSample - 1) - theta0));

    x1 = theta(:, iSample) - (alpha * theta0 + (1 - alpha) * theta(:, iSample - 1) - alpha * real(Fa * (da_cal - da_obs - Ga * (theta(:, iSample - 1) - theta0))));
    % x1 = theta(:, iSample) - (theta(:, iSample - 1) + alpha * dtheta);
    lnprop01(iSample) = (-0.5 * x1' * (inv(COV)) * x1);

    [Uz_estimated(:, iSample), dUz_estimated(:, :, iSample)] = ForwardModel(theta(:, iSample), nelmt, nu, rho, sdamp, h, etype, ww, rr, RDisk, qDisk, CR, facUz);
    Ga = [dUz_estimated(1:nfreq * Nmeasured, :, iSample); ...
        conj(dUz_estimated(1:nfreq * Nmeasured, :, iSample)); ...
        dUz_estimated(nfreq * Nmeasured + 1:end, :, iSample)];
    d = Uz_estimated(1:nfreq*Nmeasured, iSample);
    cd = conj(Uz_estimated(1:nfreq * Nmeasured, iSample));
    dre = Uz_estimated(nfreq*Nmeasured+1:end, iSample);
    da_cal = [d; cd; dre];

    % H = real(Ga'*inv_Ca_D*Ga) + inv_Ca_M;
    % g = real(Ga'*inv_Ca_D*(da_cal - da_obs)) + inv_Ca_M * (theta(:, iSample) - theta0);
    % Ha1 = H \ (Ga' * inv_Ca_D * Ga);
    % Ha2 = eye(nelmt-3) - Ha1;
    % COV1 = real(2*Ha2*Ca_M*Ha2'+Ha1*inv(H)');
    % COV1 = beta * 0.5 * (COV1 + COV1)';
    % dtheta = -H \ g;

    K = Ca_D + Ga * Ca_M * Ga';
    Fa = Ca_M * Ga' * inv(K);
    IFaGa = (eye(nelmt - 3) - Fa * Ga);
    COV1 = real(2*IFaGa*Ca_M*IFaGa'+Fa*Ca_D*Fa');
    COV1 = beta * 0.5 * (COV1 + COV1)';

    lnpost1(iSample, 1) = -(0.5 * real((da_cal-da_obs)' * inv_Ca_D * (da_cal - da_obs)) ...
        +0.5 * (theta(:, iSample) - theta0)' * inv_covtheta0 * (theta(:, iSample) - theta0));

    % x1 = theta(:, iSample-1) - (theta(:, iSample) + alpha * dtheta);
    x1 = theta(:, iSample-1) - (alpha * theta0 + (1 - alpha) * theta(:, iSample) - alpha * real(Fa * (da_cal - da_obs - Ga * (theta(:, iSample) - theta0))));
    lnprop10(iSample) = (-0.5 * x1' * (inv(COV1)) * x1);

    gamma(iSample) = min(1, det(COV1 \ COV2)^0.5*exp(lnpost1(iSample) - lnpost0(iSample) + lnprop10(iSample) - lnprop01(iSample)));

    u = unifrnd(0, 1);
    a = a + 1;
    if (u > gamma(iSample))

        theta(:, iSample) = theta(:, iSample-1);
        Uz_estimated(:, iSample) = Uz_estimated(:, iSample-1);
        dUz_estimated(:, :, iSample) = dUz_estimated(:, :, iSample-1);
        a = a - 1;
    end


    if rem(iSample, 100) == 1
        close all
        figure
        if (iSample > 450)
            for i = 1:nelmt - 4
                subplot(5, 8, i)
                histogram(theta(i, 400:iSample))
                drawnow
            end
            re_er(mean((theta(:, 400:iSample)')), Cs_true(1:nelmt - 3)/100)
        end

        figure
        subplot(2, 1, 1)
        if (iSample > 450)
            MM = theta(:, 400:iSample);
            Mhat = repmat(mean(theta(:, 400:iSample), 2), 1, iSample-400+1);
            COV1 = 1 / (iSample - 400) * (MM - Mhat) * (MM - Mhat)';
            imagesc(COV1)
            colorbar
            axis square

            drawnow
        end
        subplot(2, 1, 2)
        imagesc(COV/beta)
        colorbar
        axis square
        drawnow

    end

    disp([num2str(iSample - 1), '^{th} sample, accepted ratio: ', num2str((a) / (iSample - 1) * 100), '%'])
    disp([' - gamma=',num2str(det(COV1 \ COV2)^0.5),'*exp(', num2str(lnpost1(iSample) - lnpost0(iSample) + lnprop10(iSample) - lnprop01(iSample)), ')'])
    disp([' - log(posteriror distribution)=',num2str(lnpost0(iSample))])
    
    toc
end
