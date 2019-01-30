function dydt = ode_cont_attractor(t, y, p, t_vec, noise_vec)
%ODE_CONT_ATTRACTOR integrates the continuous attractor differential
%equations

% interpolation of noise vector
noise_vec_temp = cat(1,interp_faster(t_vec, noise_vec(:,1:3), t)',zeros(9,1));

% integrate
dydt = (- y + p.W*y + p.DC + p.stim_in.*p.stim_vec(t) + p.PV_in.*p.PV_pert_vec(t) + noise_vec_temp)./p.taus;