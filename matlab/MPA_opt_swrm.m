%% CMPA Optimization with Particle Swarm
% Brute Forcing this
clc; clearvars; close all;

if ~exist('logs', 'dir')
    mkdir('logs');  % Create directory if it does not exist
end

logfile = nextLogFile();
fid = fopen(logfile,'w');   % open for writing
fprintf(fid, "Starting run: %s\n\n", datetime("now"));
fprintf("Available Cores: %d\n", str2double(getenv('SLURM_CPUS_PER_TASK')));
 
% if isempty(gcp('nocreate'))
%     parpool('local', str2double(getenv('SLURM_CPUS_PER_TASK')));
% end
 
delete(gcp('nocreate'))
parpool('local', 6)
 
 
%% Fixed Parameters
freq              = 2305e6;               
metal_thickness   = 0.0350e-3;               % copper thickness [m]
dielectric_thick  = 0.0994e-3;               % substrate thickness [m]
eps_r             = 4;                     % dielectric constant
tan_delta         = 0;                       % loss tangent
maxBoardW         = 50e-3;                  % for area normalization [m]
 
% Material objects
DielectricLayer = dielectric("Name","Custom", ...
    "EpsilonR",eps_r,"LossTangent",tan_delta);
Copper = metal("Name","Copper", ...
    "Conductivity",59.6e6,"Thickness",metal_thickness);
 
%% Optimization using Particle Swarm
lb = [10e-3, 0.01, -0.9];
ub = [maxBoardW, 0.5, 0.9];
 
opts = optimoptions("particleswarm", ...
    "Display","iter", ...
    "SwarmSize",400, ...
    "MaxIterations",1000, ...
    "UseParallel",true, ...
    "OutputFcn",@(varargin) outfun(varargin{:},fid));
 
fprintf("Starting optimization...\n");
% fflush(stdout);
 
[x_opt, fval] = particleswarm(@(x) costFunc(x, freq, metal_thickness, dielectric_thick, DielectricLayer, Copper, maxBoardW), ...
    3, lb, ub, opts);
 
Board_W_opt   = x_opt(1);
PatchFrac_opt = x_opt(2);
FeedFrac_opt  = x_opt(3);
 
fprintf("Optimized variables:\n");
fprintf("  Board_W   = %.3f mm\n", Board_W_opt*1e3);
fprintf("  PatchFrac = %.3f\n", PatchFrac_opt);
fprintf("  FeedFrac  = %.3f\n\n", FeedFrac_opt);
 
 
%% Evaluate optimized design
pcb_opt = buildAntenna(Board_W_opt, PatchFrac_opt, FeedFrac_opt, ...
    freq, metal_thickness, dielectric_thick, DielectricLayer, Copper);
GdBi = pattern(pcb_opt, freq, 0, 0, 'Type','gain');
Zopt     = impedance(pcb_opt, freq);
s11_opt  = sparameters(pcb_opt, freq).Parameters(1,1,1);
RLopt    = -20*log10(abs(s11_opt));
ARopt    = axialRatio(pcb_opt, freq, 0, 0);
 
fprintf(fid, "\nOptimized variables:\n");
fprintf(fid, " Board_W = %.10f mm\n", Board_W_opt*1e3);
fprintf(fid, " PatchFrac = %.10f\n", PatchFrac_opt);
fprintf(fid, " FeedFrac = %.10f\n\n", FeedFrac_opt);
 
fprintf(fid, "Optimized performance:\n");
fprintf(fid, " Gain = %.2f dBi\n", GdBi);
fprintf(fid, " |Z| = %.1f Ω (Z = %.1f + j%.1f Ω)\n", abs(Zopt), real(Zopt), imag(Zopt));
fprintf(fid, " RL = %.2f dB\n", RLopt);
fprintf(fid, " AR(b) = %.2f dB\n\n", ARopt);
 
fclose(fid);
 
%% Cost Function with S11, Area penalty, and AR target
function c = costFunc(x, freq, t_met, t_diel, dielLayer, metalLayer, maxW)
 
    % Target:
    % Target values for optimization
    targetGain = 5;  % dBi
    targetRL = -10;  % dB
    targetAR = 3;    % dB
 
 
    BW = x(1); PF = x(2); FF = x(3);
    pcb = buildAntenna(BW, PF, FF, freq, t_met, t_diel, dielLayer, metalLayer);
    P    = pattern(pcb, freq, 0:10:360, 0:10:180, 'Type','gain'); GdBi = max(P(:));
    Zmag = abs(impedance(pcb, freq));
    ARdB = axialRatio(pcb, freq, 0, 0);
    RLdB = -20*log10(abs(sparameters(pcb, freq).Parameters(1,1,1)));
    % Normalize & cost terms
    c_gain = max(0, (targetGain - GdBi)/5);      % reward if gain < 5 dBi
    c_imp  = abs(Zmag - 50)/50;         % impedance mismatch
    c_ar   = max(0, (ARdB - targetAR)/3);      % AR > 3 dB penalty
    c_rl   = max(0, (targetRL - RLdB))/10;    % RL < -10 dB penalty
    c_area = (BW/maxW)^2;               % normalized area
    % Weights
    w1 = 8;   % gain
    w2 = 1;   % impedance
    w3 = 10;   % axial ratio (most important for CP)
    w4 = 1;   % return loss
    w5 = 0.5; % area (least important)
    c = w1*c_gain + w2*c_imp + w3*c_ar + w4*c_rl + w5*c_area;
 
end
 
%% Output function for logging
function stop = outfun(varargin)
    stop = false;
    if nargin < 4 || isempty(fid)
        return;
    end
    if strcmp(state,'iter') || strcmp(state,'done')
        area_mm2 = (x(1)*1e3)^2;
        line = sprintf('Iter %3d: W=%.2fmm, PF=%.3f, FF=%.3f, Area=%.1fmm^2, fval=%.2e\n', ...
            optimValues.iteration, x(1)*1e3, x(2), x(3), area_mm2, optimValues.bestfval);
        fprintf(fid, '%s', line);   % log file
        fflush(fid);
        fprintf(line);
        fflush(stdout);
        fprintf('%s', line);        % still print to screen
    end
end
%% Antenna Builder
function pcb = buildAntenna(Board_W, PFr, FFr, freq, t_met, t_diel, dielLayer, metalLayer)
    pcb = pcbStack; pcb.BoardThickness = 2*t_met + t_diel;
    B = antenna.Rectangle; B.Center = [0 0]; B.Length = 1.5*Board_W; B.Width = B.Length;
    pcb.BoardShape = B;
    P = antenna.Rectangle; P.Center = [0 0]; P.Length = Board_W; P.Width = Board_W;
    C1 = antenna.Rectangle; C1.Center = [Board_W/2, Board_W/2];
    C1.Length = Board_W*sqrt(2)*PFr; C1.Width = C1.Length;
    C1 = rotate(C1, 45, [C1.Center,-1], [C1.Center,1]); P = P - C1;
    C2 = antenna.Rectangle; C2.Center = [-Board_W/2, -Board_W/2];
    C2.Length = C1.Length; C2.Width = C1.Length;
    C2 = rotate(C2, 45, [C2.Center,-1], [C2.Center,1]); P = P - C2;
    Ant = P;
    G = antenna.Rectangle; G.Length = B.Length; G.Width = B.Width;
    pcb.Layers = {Ant, dielLayer, G}; pcb.Conductor = metalLayer;
    pcb.FeedLocations = [[ -FFr*Board_W/2, FFr*Board_W/2], 1];
    pcb.FeedDiameter = 0.00035; pcb.ViaDiameter = 0.0003;
    pcb.FeedViaModel = 'square'; pcb.FeedVoltage = 20; pcb.FeedPhase = 0;
end
 
%% Environement setup
function logfile = nextLogFile()
    % Find next logs/run_log_n.txt file
    n = 1;
    while isfile(sprintf('logs/run_log_%d.txt',n))
        n = n + 1;
    end
    logfile = sprintf('logs/run_log_%d.txt',n);
end