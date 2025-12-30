% Simulates spin to orbital angular momentum conversion in k-space Stokes vectors and Mueller matrix elements
clear all; close all;

% Constants
NA = 0.86; % Numerical aperture of the microscope objective
lambda_nm = 450; % Wavelength in nm
k0 = (2 * pi / (lambda_nm * 1e-9)); % Free-space wavevector (m^-1)
k_max = k0 * NA; % Maximum wavevector in momentum space

theta = asin(NA); % Focusing angle corresponding to NA
phi = linspace(0, 2 * pi, 500); % Azimuthal angle (full circle)

d_f = 0.6; % Linear diattenuation due to focusing
delta_f = 0.6; % Linear retardance due to focusing
d_wpc = 0.25; % Linear diattenuation of WPC
delta_wpc = 0.15; % Linear retardance of WPC

% Mueller matrix for tight focusing (azimuthal linear diattenuating retarder)
function M = mueller_focusing(phi, d_f, delta_f)
    x_f = sqrt(1 - d_f^2);
    M = [1, d_f * cos(2 * phi), d_f * sin(2 * phi), 0;
         d_f * cos(2 * phi), cos(2 * phi).^2 + x_f * sin(2 * phi).^2 * cos(delta_f), (1 - x_f * cos(delta_f)) * cos(2 * phi) * sin(2 * phi), -x_f * sin(delta_f) * sin(2 * phi);
         d_f * sin(2 * phi), (1 - x_f * cos(delta_f)) * cos(2 * phi) * sin(2 * phi), sin(2 * phi).^2 + x_f * cos(2 * phi).^2 * cos(delta_f), x_f * sin(delta_f) * cos(2 * phi);
         0, x_f * sin(delta_f) * sin(2 * phi), -x_f * sin(delta_f) * cos(2 * phi), x_f * cos(delta_f)];
end

% Mueller matrix for WPC (linear diattenuating retarder with fixed axis)
function M = mueller_wpc(d_wpc, delta_wpc)
    x_wpc = sqrt(1 - d_wpc^2);
    M = [1, d_wpc, 0, 0;
         d_wpc, 1, 0, 0;
         0, 0, x_wpc * cos(delta_wpc), x_wpc * sin(delta_wpc);
         0, 0, -x_wpc * sin(delta_wpc), x_wpc * cos(delta_wpc)];
end

% Input Stokes vectors
S_LCP = [1; 0; 0; -1]; % Left circularly polarized
S_H = [1; 1; 0; 0]; % Horizontal linear polarization
S_V = [1; -1; 0; 0]; % Vertical linear polarization
S_U = [1; 0; 0; 0]; % Unpolarized

% Initialize arrays for Stokes vectors and Mueller matrix elements
S1_LCP = zeros(size(phi));
S2_LCP = zeros(size(phi));
S3_H = zeros(size(phi));
S3_V = zeros(size(phi));
S3_U = zeros(size(phi));
M_elements = struct();
M11_raw = zeros(size(phi)); % Store raw M11 before normalization
for i = 1:4
    for j = 1:4
        M_elements.(['M', num2str(i), num2str(j)]) = zeros(size(phi));
    end
end

% Compute Mueller matrix and Stokes vectors for each azimuthal angle
for i = 1:length(phi)
    M_f = mueller_focusing(phi(i), d_f, delta_f);
    M_wpc = mueller_wpc(d_wpc, delta_wpc);
    M = M_wpc * M_f; % Combined Mueller matrix
    
    % Store raw M11 before normalization
    M11_raw(i) = M(1,1);
    
    % Normalize the Mueller matrix to ensure M00 = 1
    if M(1,1) ~= 0
        M = M / M(1,1);
    end
    
    % Compute output Stokes vectors
    S_out_LCP = M * S_LCP;
    S_out_LCP = S_out_LCP / S_out_LCP(1); % Normalize by S0
    S1_LCP(i) = S_out_LCP(2);
    S2_LCP(i) = S_out_LCP(3);
    
    S_out_H = M * S_H;
    S_out_H = S_out_H / S_out_H(1);
    S3_H(i) = S_out_H(4);
    
    S_out_V = M * S_V;
    S_out_V = S_out_V / S_out_V(1);
    S3_V(i) = S_out_V(4);
    
    S_out_U = M * S_U;
    S_out_U = S_out_U / S_out_U(1);
    S3_U(i) = S_out_U(4);
    
    % Store Mueller matrix elements
    for j = 1:4
        for k = 1:4
            M_elements.(['M', num2str(j), num2str(k)])(i) = M(j, k);
        end
    end
end

% Simulate k-space arc segments
kx = k_max * cos(phi) * 1e-6; % Convert to µm^-1
ky = k_max * sin(phi) * 1e-6;

% Plotting Stokes vectors
figure('Position', [100, 100, 1200, 800]);
subplot(2, 3, 1);
scatter(kx, ky, 50, S1_LCP, 'filled');
colormap('parula');
colorbar;
title('S1 (LCP input)');
xlabel('kx (µm^{-1})');
ylabel('ky (µm^{-1})');
axis equal;
xlim([min(kx) - 1, max(kx) + 1]);
ylim([min(ky) - 1, max(ky) + 1]);

subplot(2, 3, 2);
scatter(kx, ky, 50, S2_LCP, 'filled');
colormap('parula');
colorbar;
title('S2 (LCP input)');
xlabel('kx (µm^{-1})');
ylabel('ky (µm^{-1})');
axis equal;
xlim([min(kx) - 1, max(kx) + 1]);
ylim([min(ky) - 1, max(ky) + 1]);

subplot(2, 3, 3);
scatter(kx, ky, 50, S3_H, 'filled');
colormap('parula');
colorbar;
title('S3 (H input)');
xlabel('kx (µm^{-1})');
ylabel('ky (µm^{-1})');
axis equal;
xlim([min(kx) - 1, max(kx) + 1]);
ylim([min(ky) - 1, max(ky) + 1]);

subplot(2, 3, 4);
scatter(kx, ky, 50, S3_V, 'filled');
colormap('parula');
colorbar;
title('S3 (V input)');
xlabel('kx (µm^{-1})');
ylabel('ky (µm^{-1})');
axis equal;
xlim([min(kx) - 1, max(kx) + 1]);
ylim([min(ky) - 1, max(ky) + 1]);

subplot(2, 3, 5);
scatter(kx, ky, 50, S3_U, 'filled');
colormap('parula');
colorbar;
title('S3 (U input)');
xlabel('kx (µm^{-1})');
ylabel('ky (µm^{-1})');
axis equal;
xlim([min(kx) - 1, max(kx) + 1]);
ylim([min(ky) - 1, max(ky) + 1]);

% Plot Mueller matrix elements
figure('Position', [100, 100, 1200, 1200]);
t = tiledlayout(4, 4, 'TileSpacing', 'compact');
for i = 1:4
    for j = 1:4
        nexttile;
        if i == 1 && j == 1
            % Plot raw M11 with color scale between 0 and 1
            scatter(kx, ky, 50, M11_raw, 'filled');
            colormap('parula');
            colorbar;
            caxis([0 1]); % Set color scale between 0 and 1 as requested
        else
            % Plot other normalized Mueller matrix elements
            scatter(kx, ky, 50, M_elements.(['M', num2str(i), num2str(j)]), 'filled');
            colormap('parula');
            colorbar;
        end
        title(['M', num2str(i), num2str(j)]);
        xlabel('kx (µm^{-1})');
        ylabel('ky (µm^{-1})');
        axis equal;
        xlim([min(kx) - 1, max(kx) + 1]);
        ylim([min(ky) - 1, max(ky) + 1]);
    end
end
