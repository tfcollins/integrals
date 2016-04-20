% Coverage Probability for 3GPP Patterns
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
maxOutter = 500;
tol = 1e-6;

P_t = 1;
g_2 = 0.1;
omegas_3db = [90,45,20]*pi/180;

%% Generate additional terms
genVars = false;
if genVars
    [theta_1s, g1s] = deal(zeros(size(omegas_3db)));
    disp('Generating Gain Variables');
    for om = 1:length(omegas_3db)
        [theta_1s(om), g1s(om)] = find_3gppAntenna_Vars(g_2,omegas_3db(om));
    end
    disp('Done Generating Gain Variables');
    save('3gpp_gains.mat','g1s','theta_1s')
else
    l = load('3gpp_gains.mat');
    theta_1s = l.theta_1s;
    g1s = l.g1s;
end

%% PPP 3GPP Case
p_s_g3pp = gen_ppp_3gpp(lambda,thetaDB,omegas_3db,g_2,g1s,theta_1s);
disp('PPP case done');
figure;
plot(thetaDB, p_s_g3pp);

%% GPP 3GPP Case
p_int = zeros(length(thetaDB),length(omegas_3db));

for o = 1:length(omegas_3db)
    
    omega_3db = omegas_3db(o);
    theta_1 = theta_1s(o);
    g_1 = g1s(o);
    
    c_b = (3/10)*1/omega_3db^2;
    
    % Check
    func = @(theta,g1,omega_3db) g_1.*10.^(-3/10 .* (theta./omega_3db).^2 );
    G = g_2*(pi-theta_1)/pi + integral( @(theta) func(theta,g_1,omega_3db), 0, theta_1); % Should equal near 1
    
    tic;
    parfor i=1:length(thetaDB)
        
        theta = db2lin(thetaDB(i));
        
        zeta = @(v,s,g1,g2) (1+theta/g_1.*(P_t.*g1.*g2).*(v./s).^beta); % prevent divide by zero error
        
        chi = @(g) 1./(2*pi*sqrt(c_b*log(10))) .* 1./(g.*sqrt(log(g_1./g)));
        
        M = @(v,s,g1,g2) s.^(jj).*exp(-s)./zeta(v,s,g1,g2);
        S = @(v,s,g1,g2) s.^(ii).*exp(-s)./zeta(v,s,g1,g2);
        
        M_int = @(v,g1,g2) prod( 1./factorial(jj) .* integral( @(s) M(v,s,g1,g2), v, inf,'RelTol',tol,'AbsTol',tol,'ArrayValued',true));
        S_int = @(v,g1,g2) sum( v.^(ii)           ./ integral( @(s) S(v,s,g1,g2), v, inf,'RelTol',tol,'AbsTol',tol,'ArrayValued',true));
        
        u = (pi-theta_1)/pi;
        
        SM_prod =  @(v,g1) chi(g1).*S_int(v,g1,1).*M_int(v,g1,1);

        e_g = @(v) integral( @(g_o)   SM_prod(v,g_o), g_2, g_1,'RelTol',tol,'AbsTol',tol,'ArrayValued',true);
        
        psi = @(v) u.*S_int(v,g_2,1).*M_int(v,g_2,1) + e_g(v);
        
        p_r = @(v) exp(-v) .* psi(v);
        
        p_int(i,o) = integral( p_r, 0, maxOutter,'RelTol',tol,'AbsTol',tol, 'ArrayValued', true);
        
    end
    toc
end
%save('data.mat','p_int');

%% Generate PPP Omni (pattern doesnt matter)
ppp_p_3gpp_omni = gen_ppp_kh(lambda,thetaDB,2*pi,g_2);

%% Plots
% Load GPP Omni case
r = load('omni_gpp.mat');
thetaDB20 = -10:1:20;

markers = {'g-','r--','b-.',   'g-o','r--o','b-.o'};

for x=1:length(omegas_3db)
    hold on;
    plot(thetaDB,p_int(:,x),markers{x} ,thetaDB, p_s_g3pp(:,x), markers{x+3});
    hold off;
end
% Add omni cases
hold on;
plot(thetaDB20,r.ref,'k--+');
plot(thetaDB20,ppp_p_3gpp_omni,'k-+');
hold off;

% Create legend labels
legStr = {};
for j=1:length(omegas_3db)
    legStr = {legStr{:},['GPP \theta_{3db}=',num2str(omegas_3db(j)*180/pi),'^\circ']};
    legStr = {legStr{:},['PPP \theta_{3db}=',num2str(omegas_3db(j)*180/pi),'^\circ']};    
end
legStr = {legStr{:},'GPP Omni','PPP Omni'};
legend(legStr{:},'Location','SouthWest');
title(['3GPP g_2=',num2str(g_2)]);
xlabel('\gamma (dB)');ylabel('p_s(\gamma)');
grid on;


