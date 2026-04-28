function [x_sols, G_num_list, x_0s, n, err_vec] = mode_i_num(FRF_list, range, f_vect)

OM_vect = 2*pi.*f_vect;
f_res = length(f_vect);
f_max = f_vect(end);
f_in_range = floor(f_res/f_max);
% x = [om_i,csi_1,A1,A2,A3,RL1,RL2,RL3,RH1,RH2,RH3];

n_ch = length(FRF_list);
G_exp = zeros(length(OM_vect), n_ch);
for jj = 1:n_ch
    G_exp(:,jj) = FRF_list{jj}(OM_vect)';
    [mode(jj).ampl,mode(jj).loc] = findpeaks(abs(G_exp(:,jj))); % Amplitude and location of the resonance peaks
end % mode_ampl(ii,jj) is mode ii of pair jj
n_modes = length(mode(1).loc);

x_sols = zeros(2 + 3*n_ch, n_modes);
x_0s = x_sols;
% MODE ii
G_num_list = cell(n_modes, n_ch);
for ii = 1:n_modes % I do this, for the same mode of different pairs
    x0_ii = zeros (2 + 3*n_ch, 1);
    x0_ii(1) = 2*pi*f_vect(mode(1).loc(ii)); % This should be the same for all the pairs

    d_dom = (angle(FRF_list{1}(OM_vect(mode(1).loc(ii)+1)))- ...
            angle(FRF_list{1}(OM_vect(mode(1).loc(ii)))))/ (OM_vect(mode(1).loc(ii)+1) ...
            -OM_vect(mode(1).loc(ii))); % derivate of phase of mode ii for the jj pair
            % Estimation of x0

    x0_ii(2) = -1/(d_dom*OM_vect(mode(1).loc(ii)));

    err_vec_ii = @(x) 0;
    for jj = 1:n_ch % For each pair given (channel)
        G_num_list{ii, jj} = @(x, OM) x(2+jj)./(-OM.^2 + 1i*2*x(2)*x(1).*OM+x(1)^2) + x(2+jj+n_ch)./OM.^2 + x(2+jj+2*n_ch);
        err_vec_old = @(x)err_vec_ii(x);
        err_vec_ii = @(x) [err_vec_old(x); G_exp(mode(jj).loc(ii)-f_in_range*range : mode(jj).loc(ii)+f_in_range*range,jj) ...
            - G_num_list{ii,jj}(x,OM_vect(mode(jj).loc(ii)-f_in_range*range : mode(jj).loc(ii) + ...
            f_in_range*range))']; % global error vector of mode ii
        x0_ii(2+jj)   = FRF_list{jj}(x0_ii(1))*2*x0_ii(2)*(x0_ii(1)).^2*1i;
        x0_ii(2+jj+n_ch) = 0;
        x0_ii(2+jj+2*n_ch) = 0;
    end

    % err_ii = @(x) err_vec_ii(x)'*err_vec_ii(x);
    x_0s  (:,ii) = x0_ii;
    x_sols(:,ii) = lsqnonlin(err_vec_ii, x0_ii);
    err_vec(:,ii) = err_vec_ii(x_sols(:,ii));
end

% Now i have a solution for each mode, the global FRF is the sum of the
% G_num of a pair, summed
n.ch = n_ch;
n.modes = n_modes;
end