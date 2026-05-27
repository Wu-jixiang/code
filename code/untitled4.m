clear; close all; clc;

% ======================
% 1. Define Berea sandstone basic parameters (from Gurevich et al., 2010)
% ======================
K_g = 39e9;          % grain bulk modulus (Pa)
mu_g = 32e9;         % grain shear modulus (Pa)
rho_g = 2653;        % grain density (kg/m^3)
rho_dry = 3200;      % dry rock density (kg/m^3)

% Experimentally measured dry rock velocity-pressure data points (Han et al., 1986)
P_data_MPa = [5, 10, 20, 30, 40, 50]; % pressure points (MPa)
Vp_dry = [3280, 3380, 3520, 3610, 3660, 3690]; % dry P-wave velocity (m/s)
Vs_dry = [2160, 2230, 2320, 2370, 2400, 2415]; % dry S-wave velocity (m/s)

% Compute dry rock moduli K_dry, mu_dry from velocity data
[K_dry_data, mu_dry_data] = calcDryModuliFromVel(Vp_dry, Vs_dry, rho_dry);

% Estimate model parameters Kh, phi_c0, P_h, alpha (based on Shapiro, 2003)
[K_h, phi_c0, P_h, alpha, theta_c] = estimateShapiroParams(P_data_MPa, K_dry_data, K_g, mu_g);

% ======================
% 2. Define fluids and frequency range
% ======================
water.Kf = 2.2e9;   water.rho_f = 1031;   water.eta = 1e-3;   water.name = 'Water';
gas.Kf   = 0.0022e9; gas.rho_f = 10.8;     gas.eta   = 11e-6;  gas.name   = 'Gas';
fluid_cases = {water, gas};


chi_values = [-0.3, 0, 0.5];  
legend_chi = {'Negative coupling scenarios', 'Gurevich', 'Positive coupling scenarios'};


P_plot_MPa = 20;
freq = logspace(0, 10, 1000); 
omega = 2 * pi * freq;

% ======================
% 3. Compute and plot
% ======================
for fluid_idx = 1:length(fluid_cases)
    fluid = fluid_cases{fluid_idx};
    
    figure('Position', [100, 100, 1200, 500], ...
        'Name', sprintf('Gurevich (2010) Fig%d: %s-saturated Berea sandstone (20 MPa) heterogeneity comparison', fluid_idx+1, fluid.name));
    
    % Preallocate
    K_sat_real = zeros(length(freq), length(chi_values));
    Qinv = zeros(length(freq), length(chi_values));
    
    P_MPa = P_plot_MPa;   % 20 MPa
    [K_dry_P, mu_dry_P, phi_c_P] = getRockPropertiesAtPressure(...
        P_MPa, P_data_MPa, K_dry_data, mu_dry_data, K_h, phi_c0, P_h);
    
    for c_idx = 1:length(chi_values)
        chi = chi_values(c_idx);
        for w_idx = 1:length(omega)
            [K_sat_complex, ~] = squirtFlowModel_Gurevich(...
                omega(w_idx), fluid.Kf, fluid.eta, ...
                K_g, K_dry_P, mu_dry_P, phi_c_P, alpha, theta_c, K_h, chi);
            K_sat_real(w_idx, c_idx) = real(K_sat_complex);
            Qinv(w_idx, c_idx) = imag(K_sat_complex) / real(K_sat_complex);
        end
    end
    
    % ----- Subplot 1: Bulk modulus dispersion -----
    subplot(1,2,1); hold on; grid on; box on;
    line_styles = {'--', '-', '-.'};  
    colors = {'k', 'k', 'k'};         
    for c_idx = 1:length(chi_values)
        plot(freq, K_sat_real(:,c_idx)/1e9, ...
            'LineStyle', line_styles{c_idx}, 'Color', colors{c_idx}, ...
            'LineWidth', 2, 'DisplayName', legend_chi{c_idx});
    end
    set(gca, 'XScale', 'log');
    xlabel('Frequency (Hz)'); ylabel('Bulk modulus {\it K}_{\rm sat} (GPa)');
    title(sprintf('(%s) %s-saturated (20 MPa): bulk modulus dispersion', char(96+fluid_idx*2-1), fluid.name));
    legend('Location','best'); xlim([1e0,1e10]);
    
    % ----- Subplot 2: Attenuation -----
    subplot(1,2,2); hold on; grid on; box on;
    for c_idx = 1:length(chi_values)
        plot(freq, Qinv(:,c_idx), ...
            'LineStyle', line_styles{c_idx}, 'Color', colors{c_idx}, ...
            'LineWidth', 2, 'DisplayName', legend_chi{c_idx});
    end
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('Frequency (Hz)'); ylabel('{\it Q}^{-1}');
    title(sprintf('(%s) %s-saturated (20 MPa): attenuation', char(96+fluid_idx*2), fluid.name));
    xlim([1e0,1e10]);
    if strcmp(fluid.name,'Water'); ylim([1e-3,1e0]); else; ylim([1e-2,1e1]); end
    legend('Location','best');
    
    sgtitle(sprintf('Berea sandstone %s-saturated (20 MPa) - effect of effective pressure source term', fluid.name), 'FontWeight','bold');
end

%% ====================== Function definitions ======================
function [K_sat, mu_sat] = squirtFlowModel_Gurevich(omega, Kf, eta, K_g, K_dry, mu_dry, phi_c, alpha, theta_c, K_h, chi)
    ka = (1/alpha) * sqrt(-3i * omega * eta / Kf);  % Eq.(22)
    if abs(ka) < 1e-10
        Kf_star = 0;
    elseif abs(ka) > 100
        Kf_star = Kf * (1 + 2i./ka);               
    else
        J0ka = besselj(0, ka);
        J1ka = besselj(1, ka);
        bracket_term = 1 - (2 * J1ka) ./ (ka .* J0ka);
        Kf_star = Kf * bracket_term;
    end
    if chi ~= 0
        Kf_star = Kf_star * (1 + chi);
    end
    
    term_A = 1 / K_h;
    diff_comp = 1/K_dry - 1/K_h;
    term_B = 1 / max(diff_comp, 1e-20);
    fluid_comp = 1/Kf_star - 1/K_g;
    term_C = 1 / max(fluid_comp * phi_c, 1e-20);
    inv_K_mf = term_A + 1 / (term_B + term_C);
    K_mf = 1 / inv_K_mf;
    
    mu_mf = 1 / (1/mu_dry - (4/15) * (1/K_dry - 1/K_mf));
    
    inv_K_sat = 1/K_g + phi_c*(1/Kf_star - 1/K_g) / ...
        (1/Kf_star + phi_c*(1/Kf_star - 1/K_g) / (1/K_mf - 1/K_g));
    K_sat = 1 / inv_K_sat;
    mu_sat = mu_mf;
end

function [K_h, phi_c0, P_h, alpha, theta_c] = estimateShapiroParams(P_MPa, K_dry_data, K_g, mu_g)
    P_Pa = P_MPa * 1e6;
    inv_K_dry = 1 ./ K_dry_data;
    try
        K_h_init = K_dry_data(end);
        inv_K_h_init = 1 / K_h_init;
        y_data = inv_K_dry - inv_K_h_init;
        P_fit = P_Pa(1:min(4,end));
        log_y = log(y_data(1:min(4,end)) + eps);
        coeff = polyfit(P_fit, log_y, 1);
        P_h_init = abs(-1/coeff(1));
        phi_c0_over_P_h_init = exp(coeff(2));
        phi_c0_init = phi_c0_over_P_h_init * P_h_init;
        
        model = @(params, P) 1/params(1) + (params(2)/params(3)) * exp(-P/params(3));
        params0 = [K_h_init, phi_c0_init, P_h_init];
        lb = [K_h_init*0.5, 1e-8, 1e3];
        ub = [K_h_init*1.5, 1e-2, 1e9];
        opts = optimoptions('lsqcurvefit', 'Display', 'off');
        params_fit = lsqcurvefit(model, params0, P_Pa, inv_K_dry, lb, ub, opts);
        K_h = params_fit(1); phi_c0 = params_fit(2); P_h = params_fit(3);
        fprintf('Parameter estimation results (nonlinear fit):\n');
    catch
        warning('Optimization Toolbox not found, using linear fit.');
        K_h = K_dry_data(end);
        inv_K_h = 1 / K_h;
        y_data = inv_K_dry - inv_K_h;
        P_fit = P_Pa(1:end-1);
        log_y = log(y_data(1:end-1) + eps);
        coeff = polyfit(P_fit, log_y, 1);
        P_h = abs(-1/coeff(1));
        phi_c0_over_P_h = exp(coeff(2));
        phi_c0 = phi_c0_over_P_h * P_h;
    end
    theta_c = K_h / P_h;
    alpha = 0.01;
    fprintf('  K_h = %.2f GPa\n', K_h/1e9);
    fprintf('  phi_c0 = %.6f\n', phi_c0);
    fprintf('  P_h = %.2f MPa\n', P_h/1e6);
    fprintf('  theta_c = %.2e\n', theta_c);
    fprintf('  Aspect ratio alpha = %.4f (literature value)\n', alpha);
end

function [K_dry, mu_dry] = calcDryModuliFromVel(Vp, Vs, rho)
    Vp = Vp(:); Vs = Vs(:);
    mu_dry = rho * Vs.^2;
    K_dry = rho .* Vp.^2 - (4/3) * mu_dry;
end

function [K_dry_P, mu_dry_P, phi_c_P] = getRockPropertiesAtPressure(P_MPa, P_data_MPa, K_dry_data, mu_dry_data, K_h, phi_c0, P_h)
    P_Pa = P_MPa * 1e6;
    if P_MPa <= max(P_data_MPa) && P_MPa >= min(P_data_MPa)
        K_dry_P = interp1(P_data_MPa, K_dry_data, P_MPa, 'pchip');
        mu_dry_P = interp1(P_data_MPa, mu_dry_data, P_MPa, 'pchip');
    else
        inv_K_dry_P = 1/K_h + (phi_c0 / P_h) * exp(-P_Pa / P_h);
        K_dry_P = 1 / inv_K_dry_P;
        mu_dry_P = mu_dry_data(end) * (K_dry_P / K_h);
    end
    phi_c_P = phi_c0 * exp(-P_Pa / P_h);
end