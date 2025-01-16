nelmt = 44;

Cs = [745*ones(4,1)
    471*ones(8,1)
    408*ones(8,1)
    527*ones(24,1)];

nu = 0.25 * ones(nelmt,1);
rho = 1800 * ones(nelmt,1);
sdamp = 0.02 * ones(nelmt,1);
Cs_true=Cs;
h = [
0.05*ones(nelmt-4,1);
inf * ones(4,1)
];
RDisk = 0.1;
qDisk = 1 / (pi*RDisk^2);
etype = [
    ones(nelmt-4,1)
    21
    22
    31
    32];

ww = 2*pi * (10:10:150);

rr = (0:0.5:5);

nfreq = size(ww,2);
Nmeasured = size(rr,2);

for ii = 1:41
    theta0(ii, 1) = 7- 0.0375 * (ii - 1);
end

covtheta0 = (0.2 * 3).^2;
CR = 10^-5 * 100^2;
covUz = 0.01^2;
covdtheta = 0.5 * (0.01 / sqrt(2 * CR))^2;


nSample = 60000;
beta = 0.5;
alpha = 0.1;
