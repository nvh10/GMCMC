nelmt = 44;

Cs = [300 * ones(12, 1); ...
    400 * ones(8, 1); ...
    200 * ones(16, 1); ...
    500 * ones(4, 1); ...
    500 * ones(4, 1)];

nu = 0.35 * ones(nelmt, 1);
rho = 1800 * ones(nelmt, 1);
sdamp = 0.02 * ones(nelmt, 1);
Cs_true = Cs;
h = [; ...
    0.5 * ones(nelmt - 4, 1); ...
    inf * ones(4, 1); ...
    ];

etype = [; ...
    ones(nelmt - 4, 1); ...
    21; ...
    22; ...
    31; ...
    32];

ww = 2 * pi * (5:5:50);

rr = (5:1:19);

nfreq = size(ww, 2);
Nmeasured = size(rr, 2);

for ii = 1:41
    theta0(ii, 1) = 3 + 0.05 * (ii - 1);
end
covtheta0 = (0.2 * 3).^2;

CR = 10^-5 * 100^2;

covUz = 0.05^2;
covdtheta = 0.5 * (0.05 / sqrt(2 * CR))^2;


RDisk = 0.1;
qDisk = 1 / (pi * RDisk^2);


nSample = 60000;
beta = 0.5;
alpha = 0.1;
