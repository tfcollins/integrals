% Coverage Probability
% A cellular network model with Ginibre configured base stations
clc;clear all;

% Notation changes (old => new):
% j => jj
% i => ii

db2lin = @(db) 10.^(db./10);

% Fixed
thetaDB = -10:1:20;
beta = 2;
c = 1;
lambda = c/pi;
mu = 0.001;
alpha = 1;

% Sum, prod and integral limits
jj = 0:40;
ii = 0:40;
maxOutter = 700;

p_int = zeros(size(thetaDB));
p_int_traps= zeros(size(thetaDB));

for i=1:length(thetaDB)
    theta = db2lin(thetaDB(i));
    disp(thetaDB(i));
    
    
    M = @(v,s) (s.^jj.*exp(-s))./(1+theta.*(v./s).^beta);
    M_int = @(v) prod( 1./factorial(jj) .* integral( @(s) M(v,s), v, inf,'RelTol',1e-6,'AbsTol',1e-9,'ArrayValued',true) );
    
    
    S = @(v,s) (s.^ii.*exp(-s))/(1+theta.*(v./s).^beta);
    S_int = @(v) sum(  v.^ii            ./integral( @(s) S(v,s),v,inf,'RelTol',1e-6,'AbsTol',1e-9,'ArrayValued',true));
    
    p = @(v) exp(-v) .* M_int(v) .* S_int(v);
    
    tic;
    p_int(i) = integral( p, 0, maxOutter,'RelTol',1e-6,'AbsTol',1e-9, 'ArrayValued', true);
    t1 = toc;
    X = 0:1/10:maxOutter;
    Y = X;
    for k=1:length(X)
        Y(k) = p(X(k));
    end
    
    tic;
    p_int_traps(i) = trapz(X,Y);
    t2 = toc;
    disp([t1,t2]);
    disp([p_int(i),p_int_traps(i)]);
    
    
end

plot(thetaDB,p_int);

