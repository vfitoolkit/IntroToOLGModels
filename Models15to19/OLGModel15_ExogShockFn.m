function [z_grid,pi_z]=OLGModel15_ExogShockFn(p_jobsep,m,theta,gamma,agej,Jr,rho_z,sigma_z_epsilon,n_z1)
% p_jobfind: job finding probability
% p_jobsep: job separation probability

%% idiosyncratic worker productivity
[z1_grid,pi_z1]=discretizeAR1_FarmerToda(0,rho_z,sigma_z_epsilon,n_z1);
z1_grid=exp(z1_grid);

%% employment status, z2
p_jobfind=m*theta^(1-gamma);

z2_grid=[0;1]; % 0=not-employed, 1=employed

pi_z2=[1-p_jobfind,p_jobfind;...
    p_jobsep,1-p_jobsep];
% job finding rate = probability of transitioning from not-employed to employed
% job separation probability = probability of transitioning from employed to not-employed

if agej>=Jr-1
    pi_z2=[1,0; 1,0]; % Retired are 'listed' as not-employed
end

%% put together
z_grid=[z1_grid; z2_grid];
pi_z=kron(pi_z2,pi_z1); % note reverse order

end