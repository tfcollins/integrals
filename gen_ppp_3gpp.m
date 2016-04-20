function p_s_g3pp = gen_ppp_3gpp(lambda,thetaDB,omegas,g_2,g1s,theta_1s)
% Ps in PPP network

alpha = 4;

db2lin = @(db) 10.^(db./10);

% Additional terms for DA
P_t = 1;

% % Load gains for 3GPP pattern
% l = load('3gpp_gains.mat');
% theta_1s = l.theta_1s;
% g1s = l.g1s;

p_s_g3pp = zeros(length(thetaDB),length(omegas));

for o = 1:length(omegas)
        
    % 3 GPP parameters
    c_b = 3/10*1/omegas(o)^2;
    g_1_3gpp = g1s(o);
    theta_1 = theta_1s(o);
    omega_3db = omegas(o);
    % Check G == 1
    func = @(theta,g1,omega_3db) g_1_3gpp.*10.^(-3/10 .* (theta./omega_3db).^2 );
    G = g_2*(pi-theta_1)/pi + integral( @(theta) func(theta,g_1_3gpp,omega_3db), 0, theta_1); % Should equal near 1
    if abs(G-1)>1e-2
        error('Bad G');
    end
    
    for T = 1:length(thetaDB)
        
        rho = @(theta,u) 1./(1+u.^(alpha/2));
        
        rho_int = @(theta,g) (g./g_1_3gpp.*theta)^(2/alpha).*integral( @(u) rho(theta,u), (g./g_1_3gpp.*theta).^(-2/alpha), inf);
        
        L_I = @(theta,v,g) exp(-pi*lambda*v*(1+rho_int(theta,g)));
        
        p_s = @(theta,g) pi*lambda *integral( @(v) L_I(theta,v,g), 0, inf);
                
        % 3GPP
        u = (pi-theta_1)/pi;
        
        chi = @(g) 1./(2*pi*sqrt(c_b*log(10))) .* 1./(g.*sqrt(log(g_1_3gpp./g)));
        r = @(g,theta) chi(g).*p_s( theta, g );
        r_int = @(theta) integral(@(g) r(g,theta), g_2, g_1_3gpp,'RelTol',10e-6,'AbsTol',10e-6,'ArrayValued',true);
        
        p_s_g3pp(T,o) = u*p_s( db2lin( thetaDB(T) ), g_2 ) + r_int( db2lin( thetaDB(T) ) );
        
    end
end



