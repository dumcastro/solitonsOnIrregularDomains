function out = animate_wave(sol, fname, varargin)
%ANIMATE_WAVE  Write an animation of the disc wave field.
%
%   OUT = ANIMATE_WAVE(SOL,FNAME) writes an animated GIF of the solution to
%   FNAME (a '.gif' path).  The frames are rasterised directly from the
%   modal field onto a Cartesian pixel grid and written as indexed images,
%   so the routine needs no toolbox, no screen and no GETFRAME: it produces
%   identical output on a headless cluster node and on a desktop.
%
%   ANIMATE_WAVE(SOL,FNAME,'Name',Value,...) accepts:
%     'Npix'       frame size in pixels                          (520)
%     'Times'      times to render; default is SOL.t
%     'FPS'        playback rate                                  (20)
%     'Gamma'      colour stretch, as in PLOT_SOLUTION          (0.65)
%     'CLim'       colour limits; default symmetric at SOL.Wdisp
%     'Colormap'   name for WAVE_COLORMAP                     ('wave')
%     'ProgressBar' draw a time bar along the bottom            (true)
%     'Video'      also write a video beside the GIF, using
%                  VIDEOWRITER when it is available             (true)
%
%   OUT lists the files written.
%
%   Example
%     sol = wave_solver_disc('T', 8, 'nT', 241);
%     animate_wave(sol, 'figures/wave_disc.gif');
%
%   See also PLOT_SOLUTION, WAVE_COLORMAP, DEMO_WAVE_DISC.

p = struct('Npix', 520, 'Times', [], 'FPS', 20, 'Gamma', 0.65, ...
           'CLim', [], 'Colormap', 'wave', 'ProgressBar', true, 'Video', true);
f = fieldnames(p);
for i = 1:2:numel(varargin)
    ix = find(strcmpi(f, varargin{i}), 1);
    if isempty(ix)
        error('animate_wave:unknownOption', 'Unknown option ''%s''.', varargin{i});
    end
    p.(f{ix}) = varargin{i+1};
end
if isempty(p.Times), p.Times = sol.t; end
tv = p.Times(:)';
nt = numel(tv);

A = sol.Wdisp;
if ~isempty(p.CLim), A = max(abs(p.CLim)); end
if ~(A > 0), A = 1; end

% ---- 256-entry palette: field colours, then background, rim, bar --------
ncol = 252;
cmap = [wave_colormap(p.Colormap, ncol)
        1.00 1.00 1.00        % 253 background outside the disc
        0.12 0.12 0.16        % 254 rim
        0.16 0.18 0.24        % 255 progress bar, elapsed
        0.88 0.88 0.90];      % 256 progress bar, track
iBG = 253; iRIM = 254; iBAR = 255; iTRK = 256;

% ---- pixel grid --------------------------------------------------------
R  = sol.opt.R;
N  = p.Npix;
xs = linspace(-1.02*R, 1.02*R, N);
[Xp, Yp] = meshgrid(xs, -xs);            % row 1 is the top of the image
rq  = hypot(Xp, Yp);
thq = mod(atan2(Yp, Xp), 2*pi);
inside = rq <= R;
dx  = xs(2) - xs(1);
rim = abs(rq - R) <= 0.9*dx;

thc = [sol.thp, 2*pi];                   % closed angular grid of SOL.X/.Y
rqc = min(rq, R);

d = fileparts(fname);
if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end

RGBs = cell(1, nt);
fprintf('    rendering %d frames at %dx%d ...\n', nt, N, N);
for k = 1:nt
    W = frame_field(sol, tv(k));
    V = interp2(thc, sol.rp, W, thq, rqc, 'linear');
    V(~inside) = 0;

    C = sign(V).*(min(abs(V)/A, 1).^p.Gamma);       % in [-1,1]
    idx = 1 + round((C + 1)/2*(ncol - 1));
    idx = min(max(idx, 1), ncol);
    idx(~inside) = iBG;
    idx(rim)     = iRIM;

    if p.ProgressBar
        rows = (N-9):(N-5);
        idx(rows, :) = iTRK;
        if nt > 1, frac = (k - 1)/(nt - 1); else, frac = 1; end
        nb = max(1, round(frac*N));
        idx(rows, 1:nb) = iBAR;
    end

    A8 = uint8(idx - 1);
    if k == 1
        imwrite(A8, cmap, fname, 'gif', 'WriteMode', 'overwrite', ...
                'LoopCount', Inf, 'DelayTime', 1/p.FPS);
    else
        imwrite(A8, cmap, fname, 'gif', 'WriteMode', 'append', ...
                'DelayTime', 1/p.FPS);
    end
    if p.Video
        RGBs{k} = reshape(cmap(idx(:), :), [N N 3]);
    end
end
out = {fname};
fprintf('    wrote %s (%d frames)\n', fname, nt);

% ---- optional video ----------------------------------------------------
if p.Video
    [d, b] = fileparts(fname);
    made = '';
    for prof = {'MPEG-4', 'Motion JPEG AVI'}
        try
            if strcmp(prof{1}, 'MPEG-4'), ext = '.mp4'; else, ext = '.avi'; end
            vf = fullfile(d, [b ext]);
            v  = VideoWriter(vf, prof{1}); %#ok<TNMLP>
            v.FrameRate = p.FPS;
            open(v);
            for k = 1:nt, writeVideo(v, RGBs{k}); end
            close(v);
            made = vf;
            break
        catch
            % profile unavailable in this release; try the next one
        end
    end
    if ~isempty(made)
        out{end+1} = made;
        fprintf('    wrote %s\n', made);
    else
        fprintf('    (VideoWriter unavailable here - GIF only)\n');
    end
end
end

% ======================================================================
function W = frame_field(sol, t)
if ~isempty(sol.W)
    [dt, i] = min(abs(sol.t - t));
    if dt < 1e-12*max(1, abs(t)), W = sol.W(:, :, i); return, end
end
W = sol.field(t);
end
