function out = scan_corr(gh,gv,varargin)
% out = scan_corr(gh,gv), correct scanning artifacts given an image pair gh
% and gv. By default, it is expected that gh is scanned horizontally and gv
% is scanned vertically. gv will be rotated 90 degrees counterclockwise
% to align it with gh (see also the options below).
%
% out is the output structure with the following fields:
% im          : the corrected image 0.5*(gh(x+uh) + gv(x+uv))
% res         : the residual image (gh(x+uh) - gv(x+uv))
% r           : the rms of res
% uh          : the row displacements applied to gh (in x and y)
% uv          : the column displacements applied to gv (in x and y)
% ind_drift   : the mean strain due to drift [hx, hy ; vx, vy]
% ind_linestd : the std of the remaining artifacts [hx, hy ; vx, vy]
% urbt        : the rbt correction [ux, uy, C] where C is the
%               cross-correlation coefficient
% where x and y are the horizontal and vertical directions respectively,
% with zero in the top left corner.
%
% Main Concept: 
% ------------------------- 
% The correction tool works on a pair of images, the first scanned
% horizontally (H), the second scanned vertically (V).
% 
% The basic idea of the artifact correction is that both images are
% incorrect, but, there exists a third unknown image that represents the
% sample without scanning artifacts. The hypothesis is that each scan line
% is less perturbed within the line compared to perturbations between
% lines. Consequently, the pixel rows of H are moved each as one, while the
% pixel columns of V are moved each as one while minimizing the difference
% between the two. The average between the finally obtained corrected
% images will be the best representation of this unknown image.
% 
% This algorithm performs a few steps sequentially: 
% - Preparation: 
% Here V is rotated and the center square of both images are cropped.
% 
% - Rigid Body Translation: 
% A single rigid body translation is obtained using FFT cross-correlation.
% Using this translation H and V are re-cropped (from the original image
% shifting H horizontally and V vertically with integer pixel shifts.
% 
% - Steps: 
% Finding the corrected image iteratively can be challenging depending on
% the input images. To robustly converge to the solution a pyramid approach
% is applied. Here, the problem is treated stepwise where to solution
% obtained for each step will be the initial guess for the next starting
% with a zero initial guess for the first step. Only the final step will be
% important for the quality of the result. The other steps serve only to
% robustly get to this final step.
%  
% In a pyramid approach, the complexity of the problem is gradually
% increased. This is done in two ways:
%   - Gradually increase the number of degrees of freedom 
%   - Gradually increase the image details
% In the options detailed below, there are 9 options that can be given as
% an array. The longest of these arrays will be used to deduce the number
% of desired steps. The shorter arrays will be padded with their last value
% to match this length.
% Of these 9 options, 4 embody directly this pyramid strategy. The first
% two, i.e. <PolyOrder> and <LineThickness> allow controlling the number of
% degrees of freedom. The following two options <Bin> and <Blur> control
% the amount of detail in the image. In general all of these 4 options
% should decrease from one step to the next.
% 
% A note on judging the quality of the corrected image:
% ------------------------- 
% The types of images that can be given to the tool can be very diverse.
% Consequently, it is not guaranteed that the tool will converge to the
% correct solution. In that case, the corrected image may contain worse
% artifacts than either of the input images. We have tried our best to
% minimize these by making the tool as robust as possible, but there are
% limits.
% The responsibility remains with the user to verify if the correction was
% a success. Here the residual images are of utmost importance. The tool
% provides options that will store the residuals in various forms. It is
% highly advised to store these residuals and to analyze them before using
% the corrected images. However, analyzing residuals can be challenging, so
% we provide some typical issues to look for here:
% - Anything that is not white noise. The perfect residual image only
% contains white (Gaussian) noise. Any deviation from this may indicate an
% issue. However, there is one exception. One or both of the images may be
% perturbed with a low (spatial) frequency brightness or contrast change
% without affecting the quality of the corrected image. This perturbation
% will be correctly removed by the tool. A typical example of this is when
% the microscopist has zoomed on the image to adjust focus and thereby
% created some local carbon deposition central in the final image. If the
% zoomed rectangle is not too small on the final image, this large-scale
% change in gray level will not affect the result too much. - Stipes and
% Bands: if a pixel row (or group of rows) failed to be correctly
% positioned it will leave a band of elevated residuals. These bands may
% intersect with similar orthogonal errors that occurred in the other image
% rendering the errors more as blocks. The texture may be seen as similar
% to plaid or tartan. - Fringes: when the material points of one image are
% not positioned correctly it will create a repeated positive/negative
% signature of the residual image. Here the missing displacement will be
% orthogonal to the fringes. - Dual vision: this is a severe case of the
% former. If a part of the sample is placed completely at the wrong
% location, it will appear twice in the residual image, once positive and
% once negative. Typically very visual when there is a bright or dark spot
% on the sample surface.
%
% Options:
% ------------------------- 
% out = scan_corr(gh,gv,opt), a third input structure can be given to
% change the options, the following fields are recognized:
%
% Verbose       : (2) change the amount of text printed to the screen, set
%                 to 0 to disable
% Plot          : (true) use false to disable plotting which will be
%                 faster. Alternatively, close the figure during the
%                 computations if so desired.
%
% Savename      : ('') the basename for all the saved files any given 
%                 extention will be removed. If empty, no files will be
%                 saved.
% TIFdepth      : (8) set to 16 to save the tif files as 16-bit
% SaveLOG       : (true) use false to disable saving the logfile
% SaveCOR       : (true) use false to disable saving the corrected iamge
% SaveRES       : (true) use false to disable saving the residual image
%                 the residual image will be shifted such that negative
%                 values will be black, zero will be gray and positive
%                 values white. Also the amplitudes are amplified with a
%                 factor of 10.
% SavePNG       : (true) use false to disable saving the final figure
% SaveRBT       : (true) use false to disable saving the two images as they
%                 are after RBT correction
% 
% RotateH       : (0) the number of 90 degrees counter clockwise rotations
%                 for the horizontal image
% RotateH       : (&) the number of 90 degrees counter clockwise rotations
%                 for the vertical image
% DatabarRemove : (true) use false to disable automatic databar removal
% ROI           : (50) the number of pixels from each border to ignore.
%                 Alternatively use [left, right, top, bottom] to define
%                 each border individually. Note, that the algorithm will
%                 still give results outside of this border, however, these
%                 pixels are note used in the optimization algorithm.
% ROI_NAN       : (false) use true to paint the pixels outside of the ROI
%                 with NaN values.
% 
% RBT_flag      : (true) use false to disable the Rigid Body Translation 
%                 correction
% RBT_Window    : (0.1) the radius of the applied Tukey window, see also
%                 the matlab function TUKEYWIN
% 
% Steps_flag    : (true) use false to disable the Pyramid Step section
% PolyOrder     : ([2, 2, 0, 0, 0]) the order of the applied polynomial
%                 regularization. Set to 0 to disable polynomial
%                 regularizatino. It is adviced to combine polynomial
%                 regularization with a LineThickness of 1.
% LineThickness : ([1, 1, 10, 3, 1]) group lines togheter and force them to
%                 move as one.
% Bin           : ([4, 4, 2, 2, 1]) bin the image into superpixels of size
%                 NxN.
% Blur          : ([8, 4, 2, 1, 0]) blur the image with a gaussian blur of
%                 radius R
% ConvCrit      : ([1e-3, 1e-3, 1e-3, 1e-3, 1e-4]) consider the step
%                 converged if the RMS of the update in the displacement is
%                 smaller than this criterium.
% IterMax       : ([30, 30, 30, 30, 40]) the maximum number of iterations
%                 per step. Note that one iteration entails a correction of
%                 the horizontal image AND a correction of the vertical
%                 image, each may contain mutliple sub-iterations.
% SubIterMax    : (3) the maximum number of sub-iterations. Note, the
%                 corrected image is not reconstructed during these
%                 sub-iterations.
% TrustRegion   : (2) the maximum distance to allow for the displacement
%                 correction per iteration
% InterpMethod  : ([3, 3, 3, 3, 5]) the degree of the BSpline image
%                 interpolation. The interpolator used in this tool is
%                 based on the work of P. Thevenaz et al. [1]
%
% Contact 
% ----------------------------------
% This tool is the result of the MATMECA workgroup collaborations. Please
% contact either of the following for information on licenses and or to
% discuss errors bugs and features. 
% Jan Neggers: jan.neggers@centralesupelec.fr 
% Eva H?ripr?: eva.heripre@centralesupelec.fr
% 
% Warranty 
% ----------------------------------
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.
% 
% References 
% -------------------------
% [1] P. Thevenaz, T. Blu, M. Unser, "Interpolation Revisited," IEEE
% Transactions on Medical Imaging, vol. 19, no. 7, pp. 739-758, July 2000.

out = [];

if nargin == 0
    error('expecting at least 2 inputs')
end

% default options
opt.Plot = true;
opt.Verbose = 2;

opt.Savename = '';
opt.TIFdepth = 8; % or 16
opt.SaveRBT = true;
opt.SaveCOR = true;
opt.SaveRES = true;
opt.SavePNG = true;
opt.SaveLOG = true;

opt.RotateH = 0; % number of 90 degrees rotations
opt.RotateV = 1; % number of 90 degrees rotations

opt.DatabarRemove = true;
opt.ROI = 50; % pixels from the boundary
opt.ROI_NAN = false;

opt.RBT_flag = true;
opt.RBT_Window = 0.1; % tukeywin

opt.Steps_flag = true;
opt.PolyOrder = [2, 2, 0, 0, 0];
opt.LineThickness = [1, 1, 10, 3, 1];
opt.Bin = [4, 4, 2, 2, 1];
opt.Blur = [8, 4, 2, 1, 0];
opt.ConvCrit = [1e-3, 1e-3, 1e-3, 1e-3, 1e-4];
opt.IterMax = [30, 30, 30, 30, 40];
opt.SubIterMax = 3;
opt.TrustRegion = 2;
opt.InterpMethod = [3, 3, 3, 3, 5];

% hidden options
opt.cmap = [];
opt.gui = false;
opt.handles = [];
opt.fontsize = 18;
opt.dpi = 200;

% overwrite options
if (nargin == 3) && ~isempty(varargin{1}) && isstruct(varargin{1})
    fields = fieldnames(varargin{1});
    for k = 1:numel(fields)
        field = (fields{k});
        if ~isempty(field)
            opt.(field) = varargin{1}.(field);
        end
    end
end

if opt.gui
    hdl_fig = opt.handles.figure;
    hdl_a = opt.handles.axes;
    hdl_i = opt.handles.images;
    hdl_p = opt.handles.plot;
    hdl_roi = opt.handles.roi;
    hdl_abort = opt.handles.abort;
    Hs = opt.handles.status;
else
    Hs = NaN;
    hdl_abort = NaN;
end

if isempty(opt.Savename)
    opt.SaveRBT = false;
    opt.SaveCOR = false;
    opt.SaveRES = false;
    opt.SavePNG = false;
    opt.SaveLOG = false;
end

if ~isempty(opt.Savename)
    [sPath, sFilename, ~] = fileparts( opt.Savename );
    opt.Savename = fullfile(sPath,sFilename);
    
    if opt.SaveLOG
        opt.Logfile = [opt.Savename  '_log.txt'];
        % clear the file
        fid = fopen(opt.Logfile,'w+t');
        fclose(fid);
    else
        opt.Logfile = [];
    end
end

if isempty(opt.cmap)
    cmap = RdGy;
else
    cmap = opt.cmap;
end


% the step options
stepopts{1,1} = 'PolyOrder';
stepopts{2,1} = 'LineThickness';
stepopts{3,1} = 'Bin';
stepopts{4,1} = 'Blur';
stepopts{5,1} = 'ConvCrit';
stepopts{6,1} = 'IterMax';
stepopts{7,1} = 'SubIterMax';
stepopts{8,1} = 'TrustRegion';
stepopts{9,1} = 'InterpMethod';

% the the longest definition
Nsteps = 1;
for i = 1:numel(stepopts)
    N = numel(opt.(stepopts{i}));
    Nsteps = max(Nsteps,N);
end

% pad all using the last value
for i = 1:numel(stepopts)
    N = numel(opt.(stepopts{i}));
    opt.(stepopts{i})(N+1:Nsteps) = opt.(stepopts{i})(N);
end

opt.Bin = round(opt.Bin);
opt.LineThickness = round(opt.LineThickness);
opt.IterMax = round(opt.IterMax);
opt.SubIterMax = round(opt.SubIterMax);
opt.InterpMethod = round(opt.InterpMethod);
if any( opt.InterpMethod < 2 ) || any( opt.InterpMethod > 9 )
    error( 'InterpMethod must be between 2 and 9 (inclusive)' )
end

% prepare a figure
if opt.Plot && ~opt.gui
    w1 = 0.1;
    w2 = 0.57;
    InterpMethod = 0.05;
    xpos(1) = 0.13;
    xpos(2) = xpos(1)+w1+InterpMethod;
    
    h1 = 0.1;
    h2 = 0.6;
    ypos(1) = 0.1;
    ypos(2) = xpos(1)+h1+InterpMethod;
    
    col1 = [     0    0.4470    0.7410];
    col2 = [0.8500    0.3250    0.0980];
    
    figpos = [50 50 830 800];
    hdl_fig = figure('Position',figpos);
    
    % set the paper position to 1 inch per 100 pixels
    set(hdl_fig,'PaperUnits','inches','PaperPosition',figpos.*[0 0 1e-2 1e-2])
    set(hdl_fig,'PaperSize',figpos(3:4).*[1e-2 1e-2])
    
    hdl_a(1) = axes('Position',[xpos(2),ypos(2),w2,h2],'FontSize',opt.fontsize);
    hdl_i = imagesc(nan);
    set(hdl_a(1),'Box','On','DataAspectRatio',[1 1 1]);
    set(hdl_a(1),'Xtick',{},'Ytick',{});
    colorbar
    colormap(hdl_a(1),cmap);
    set(hdl_a(1),'Position',[xpos(2),ypos(2),w2,h2]);
    % plot the ROI
    hdl_roi = patch('Vertices',nan(4,2),'Faces',1:4,'FaceColor','none','EdgeColor',col1,'Parent',hdl_a(1));
    
    hdl_a(2) = axes('Position',[xpos(1),ypos(2),w1,h2],'NextPlot','Add','Box','On','FontSize',opt.fontsize,'YDir','reverse');
    hdl_p(1) = plot(nan,nan,'-','Color',col1);
    hdl_p(2) = plot(nan,nan,'-','Color',col2);
    legend(hdl_p(1:2),{'$u_{hx}$','$u_{hy}$'},'Orientation','horizontal','Interpreter','Latex','Location','NorthOutside','FontSize',opt.fontsize);
    set(hdl_a(2),'Position',[xpos(1),ypos(2),w1,h2]);
    xlabel('$u_h$ [px]','Interpreter','Latex','FontSize',opt.fontsize)
    ylabel('$y$ [px]','Interpreter','Latex','FontSize',opt.fontsize)
    
    hdl_a(3) = axes('Position',[xpos(2),ypos(1),w2,h1],'NextPlot','Add','Box','On','FontSize',opt.fontsize);
    hdl_p(3) = plot(nan,nan,'-','Color',col2);
    hdl_p(4) = plot(nan,nan,'-','Color',col1);
    legend(hdl_p(3:4),{'$u_{vx}$','$u_{vy}$'},'Orientation','vertical','Interpreter','Latex','Location','EastOutside','FontSize',opt.fontsize);
    ylabel('$u_v$ [px]','Interpreter','Latex','FontSize',opt.fontsize)
    xlabel('$x$ [px]','Interpreter','Latex','FontSize',opt.fontsize)
    set(hdl_a(3),'Position',[xpos(2),ypos(1),w2,h1]);
end


% Image prep
% ==================================================

HVstr = 'HV';

if opt.Verbose > 0
    logprintf(opt.Logfile,'Preparing images \n');
    logprintf(opt.Logfile,'--------------------------------------\n');
    appendstat(Hs,'Preparing images');
end

if opt.DatabarRemove
    gh = databar_remove(gh,0.99);
    gv = databar_remove(gv,0.99);
end

% dynamic range
if isinteger(gh)
    drangeh = single(intmax(class(gh)));
else
    drangeh = max(gh(:));
end
if isinteger(gv)
    drangev = single(intmax(class(gv)));
else
    drangev = max(gv(:));
end

% convert to float
gh = single(gh);
gv = single(gv);

% convert to grayscale
gh = grayscale(gh);
gv = grayscale(gv);

% scale the image
gh = gh./drangeh;
gv = gv./drangev;

% correct the mean
gv = gv - mean(gv(:)) + mean(gh(:));

% apply rotations
gh = rot90(gh,opt.RotateH);
gv = rot90(gv,opt.RotateV);

% get the image size
sizh = size(gh);
sizv = size(gv);

% crop to the smallest length
n = min(sizh(1), sizv(1));
m = min(sizh(2), sizv(2));

% pick the centers square of length n for each image
Ih = (1:n) + floor(0.5*(sizh(1) - n));
Jh = (1:m) + floor(0.5*(sizh(2) - m));
ghc = imageblur(gh(Ih,Jh),2);
Iv = (1:n) + floor(0.5*(sizv(1) - n));
Jv = (1:m) + floor(0.5*(sizv(2) - m));
gvc = imageblur(gv(Iv,Jv),2);

res = ghc - gvc;
r = rms(res(:));

if opt.Verbose > 0
    logprintf(opt.Logfile,'   r:%10.3e\n',r);
    appendstat(Hs,'   r:%10.3e',r);
end

if opt.Plot && ishandle(hdl_fig)
    set(hdl_a(1),'XLim',[1 n],'YLim',[1 n]);
    set(hdl_i,'CData',res);
    title(hdl_a(1),'residual after cropping','FontSize',opt.fontsize);
    clim = get_clim(res(:),[0.01,0.99]);
    set(hdl_a(1),'CLim',clim);
    
    set(hdl_roi,'Vertices',nan(4,2));
    set(hdl_a(2),'YLim',[1 n]);
    set(hdl_a(3),'XLim',[1 n]);
    set(hdl_p(1),'YData',1:n,'XData',zeros(1,n));
    set(hdl_p(2),'YData',1:n,'XData',zeros(1,n));
    set(hdl_p(3),'XData',1:n,'YData',zeros(1,n));
    set(hdl_p(4),'XData',1:n,'YData',zeros(1,n));
    drawnow;
end


% Part 1 -  RBT
% ==================================================
if opt.RBT_flag
    
    if opt.Verbose > 0
        logprintf(opt.Logfile,'Part 1 - RBT \n');
        logprintf(opt.Logfile,'--------------------------------------\n');
        appendstat(Hs,'Part 1 - RBT');
    end
    
    % FFT DIC
    % -------------------
    
    % map from index to displacement
    umap = ifftshift( (1:n) - floor(n/2) - 1);
    
    if opt.RBT_Window > 0
        % apply a window
        wn = tukeywin(n,opt.RBT_Window);
        wm = tukeywin(m,opt.RBT_Window);
        W = wn * wm.';
        
        ghc = ghc.*W;
        gvc = gvc.*W;
    else
        % 1px edge blurring
        ghc([1,end],:) = 0.5*(ghc([1,end],:) + ghc([end,1],:));
        ghc(:,[1,end]) = 0.5*(ghc(:,[1,end]) + ghc(:,[end,1]));
        gvc([1,end],:) = 0.5*(gvc([1,end],:) + gvc([end,1],:));
        gvc(:,[1,end]) = 0.5*(gvc(:,[1,end]) + gvc(:,[end,1]));
    end
    
    % zero-normalize the facet (ZNCC)
    ghc = ghc - mean(ghc(:));
    gvc = gvc - mean(gvc(:));
    ghc = ghc./std(ghc(:));
    gvc = gvc./std(gvc(:));
    
    % Cross-Correlation function
    cc = ifft2( conj(fft2(ghc)) .* fft2(gvc) ) ./ (n*n);
    
    % find the maximum in the ccf (location of best correlation)
    [C, I] = max(cc(:));
    
    % get the displacement
    [Ii, Ij] = ind2sub([n,n],I);
    urbt = umap([Ij,Ii]);
    
    if C < 0.4
        warning('fft dic: low corrlation index (C = %g), perhaps the orientations of gh and gv are incorrect?',C);
        appendstat(Hs,'WARNING: fft dic: low corrlation index (C = %g), perhaps the orientations of gh and gv are incorrect?',C);
    end
    
    % Do a smart re-crop using this knowledge
    % -------------------
    
    % shift the image gh (apply neg horizontal displacement)
    Jh = Jh - urbt(1);
    
    % shift the image gv (apply vertical displacement)
    Iv = Iv + urbt(2);
    
    % extra crop if the boundary wasn't sufficiently small
    Jv = Jv(Jh >= 1 & Jh <= sizh(2));
    Ih = Ih(Iv >= 1 & Iv <= sizv(1));
    
    % avoid going outside of the image
    Jh = Jh(Jh >= 1 & Jh <= sizh(2));
    Iv = Iv(Iv >= 1 & Iv <= sizv(1));
    
    ghc = gh(Ih,Jh);
    gvc = gv(Iv,Jv);
    
    [n, m] = size(ghc);
    
    res = ghc - gvc;
    r = rms(res(:));
    
    if opt.Verbose > 0
        logprintf(opt.Logfile,'   r:%10.3e, ux: %gpx, uy: %gpx, C:%10.3e\n',r,urbt(1),urbt(2),C);
        appendstat(Hs,'   r:%10.3e, ux: %gpx, uy: %gpx, C:%10.3e',r,urbt(1),urbt(2),C);
    end
    
    if opt.Plot && ishandle(hdl_fig)
        set(hdl_a(1),'XLim',[1 m],'YLim',[1 n]);
        set(hdl_i,'CData',res);
        title(hdl_a(1),'residual after RBT Correction','FontSize',opt.fontsize);
        clim = get_clim(res(:),[0.01,0.99]);
        set(hdl_a(1),'CLim',clim);
        drawnow;
        
        set(hdl_a(2),'YLim',[1 n]);
        set(hdl_a(3),'XLim',[1 m]);
        set(hdl_p(1),'YData',1:n,'XData',urbt(1)*ones(1,n));
        set(hdl_p(4),'XData',1:m,'YData',urbt(2)*ones(1,m));
    end
    
end

if opt.SaveRBT
    if opt.TIFdepth == 8
        gh_im = uint8(double(intmax('uint8'))*ghc);
        gv_im = uint8(double(intmax('uint8'))*gvc);
    else
        gh_im = uint16(double(intmax('uint16'))*ghc);
        gv_im = uint16(double(intmax('uint16'))*gvc);
    end
    imwrite(gh_im,sprintf('%s_rbt_h.tif',opt.Savename),'TIFF','Compression','none');
    imwrite(gv_im,sprintf('%s_rbt_v.tif',opt.Savename),'TIFF','Compression','none');
end


% overwrite the original images with the cropped ones
gh = ghc;
gv = gvc;
[n, m] = size(gh);


% the ROI
if numel(opt.ROI) == 1
    opt.ROI = [opt.ROI, m-opt.ROI+1, opt.ROI, n-opt.ROI+1];
elseif numel(opt.ROI) == 4
    opt.ROI = [opt.ROI(1), m-opt.ROI(2)+1, opt.ROI(3), n-opt.ROI(4)+1];
else
    error('ROI definition is the incorrect size')
end

% force the ROI to be inside the image
opt.ROI = clamp(opt.ROI,1,[m,m,n,n]);

% Part 2 - Steps
% ==================================================

% initialize
Ah = zeros(n,2);
Av = zeros(m,2);
itcount = zeros(Nsteps,2);

if opt.Steps_flag
    
    if opt.Verbose > 0
        logprintf(opt.Logfile,'Part 2 - Steps \n');
        logprintf(opt.Logfile,'--------------------------------------\n');
        appendstat(Hs,'Part 2 - Steps');
    end
    
    for istep = 1:Nsteps
        % Per step settings
        LineThickness = max(round(opt.LineThickness(istep)),1);
        PolyOrder = opt.PolyOrder(istep);
        Bin = opt.Bin(istep);
        Blur = opt.Blur(istep);
        InterpMethod = opt.InterpMethod(istep);
        TrustRegion = opt.TrustRegion(istep);
        IterMax = opt.IterMax(istep);
        SubIterMax = opt.SubIterMax(istep);
        ConvCrit = opt.ConvCrit(istep);
        
        if opt.Verbose > 1
            logprintf(opt.Logfile,'Step %d/%d\n',istep,Nsteps);
            logprintf(opt.Logfile,'----------------------\n');
            appendstat(Hs,'Step %d/%d\n',istep,Nsteps);
        end
        if opt.Verbose > 1
            logprintf(opt.Logfile,'      %15s: %g \n','PolyOrder',PolyOrder);
            logprintf(opt.Logfile,'      %15s: %g \n','LineThickness',LineThickness);
            logprintf(opt.Logfile,'      %15s: %g \n','Bin',Bin);
            logprintf(opt.Logfile,'      %15s: %g \n','Blur',Blur);
            logprintf(opt.Logfile,'      %15s: %g \n','ConvCrit',ConvCrit);
            logprintf(opt.Logfile,'      %15s: %g \n','IterMax',IterMax);
            logprintf(opt.Logfile,'      %15s: %g \n','SubIterMax',SubIterMax);
            logprintf(opt.Logfile,'      %15s: %g \n','TrustRegion',TrustRegion);
            logprintf(opt.Logfile,'      %15s: %g \n','InterpMethod',InterpMethod);
        end
        
        if LineThickness == 0
            % disabled step, move on to the next
            continue
        end
        
        % define the linethickness and blur relative to the original image size
        LineThickness = max(LineThickness ./ Bin,1);
        Blur = Blur ./ Bin;
        
        % bin the image
        ghb = imcoarsgrain(gh,[Bin, Bin]);
        gvb = imcoarsgrain(gv,[Bin, Bin]);
        
        % bin the displacements
        Ahb = imcoarsgrain(Ah,[Bin, 1]);
        Avb = imcoarsgrain(Av,[Bin, 1]);
        Ahb = Ahb ./ Bin;
        Avb = Avb ./ Bin;
        
        % bin the coordinates (just for the plot)
        xb = imcoarsgrain(1:m,[1, Bin]);
        yb = imcoarsgrain(1:n,[1, Bin]);
        
        % blur the image
        gvb = imageblur(gvb,Blur);
        ghb = imageblur(ghb,Blur);
        
        roi = opt.ROI./Bin;
        roih = round(roi(1:2));
        roiv = round(roi(3:4));
        
        [nb, mb] = size(gvb);
        [X, Y] = meshgrid(1:mb,1:nb);
        Iroi = X >= roi(1) & X <= roi(2) & Y >= roi(3) & Y <= roi(4);
        
        if opt.Plot && ishandle(hdl_fig)
            set(hdl_roi,'Vertices',[roi([1,2,2,1]).', roi([3,3,4,4]).']);
        end
        
        % create the h-blocks
        LTh = zeros(nb,1);
        Nrh = roiv(2) - roiv(1);
        Nh = floor(Nrh / LineThickness);
        ind = floor(linspace(1,Nh+1,Nrh+1));
        LTh( roiv(1)+(1:Nrh) ) = ind(1:Nrh);
        
        % convert Ah to ah
        phih = BlockMap(LTh,Nh);
        Ixh = 1:2:2*Nh;
        Iyh = 2:2:2*Nh;
        ah = zeros(2*Nh,1);
        ah(Ixh) = phih\Ahb(:,1);
        ah(Iyh) = phih\Ahb(:,2);
        
        % create the v-blocks
        LTv = zeros(mb,1);
        Nrv = roih(2) - roih(1);
        Nv = floor(Nrv / LineThickness);
        ind = floor(linspace(1,Nv+1,Nrv+1));
        LTv( roih(1)+(1:Nrv) ) = ind(1:Nrv);
        
        % convert Av to av
        phiv = BlockMap(LTv,Nv);
        Ixv = 1:2:2*Nv;
        Iyv = 2:2:2*Nv;
        av = zeros(2*Nv,1);
        av(Ixv) = phiv\Avb(:,1);
        av(Iyv) = phiv\Avb(:,2);
        
        [f, R] = scan_corr_im_mex(ghb,gvb,Ahb,Avb,InterpMethod);
        
        if opt.Plot && ishandle(hdl_fig)
            AlphaData =  Iroi + 0.5 * not(Iroi);
            set(hdl_i,'CData',R,'AlphaData',AlphaData);
            title(hdl_a(1),sprintf('Step %d/%d',istep,Nsteps),'FontSize',opt.fontsize);
            clim = get_clim(R(Iroi),[0.01,0.99]);
            set(hdl_a(1),'CLim',clim);
            set(hdl_a(1),'XLim',[1 mb],'YLim',[1 nb]);
            
            set(hdl_a(2),'YLim',[1 n]);
            set(hdl_a(3),'XLim',[1 m]);
            set(hdl_p(1),'YData',yb,'XData',Ahb(:,1));
            set(hdl_p(2),'YData',yb,'XData',Ahb(:,2));
            set(hdl_p(3),'XData',xb,'YData',Avb(:,1));
            set(hdl_p(4),'XData',xb,'YData',Avb(:,2));
        end
        
        Sh = SensitivityPoly(phih\(yb(:)./n),Nh,PolyOrder);
        Sv = SensitivityPoly(phiv\(xb(:)./m),Nv,PolyOrder);
        
        for it = 1:IterMax
            
            for mode = [0, 1]
                if mode == 0
                    N = Nh;
                    S = Sh;
                else
                    N = Nv;
                    S = Sv;
                end
                
                subda = zeros(2*N,1);
                for subit = 1:SubIterMax
                    
                    % count the total number of it per mode
                    itcount(istep,mode+1) = itcount(istep,mode+1) + 1;
                    
                    % perform one subiteration
                    if mode == 0
                        [M, b, r] = scan_corr_mex(f,ghb,ah,LTh,roih,mode,InterpMethod);
                    else
                        [M, b, r] = scan_corr_mex(f,gvb,av,LTv,roiv,mode,InterpMethod);
                    end
                    
                    % solve
                    da = scan_corr_solve(M,b,S);
                    
                    % trustregion
                    da = sign(da) .* min(abs(da),TrustRegion);
                    
                    if mode == 0
                        ah = ah + da;
                    else
                        av = av + da;
                    end
                    
                    % add the sub updates (for conv check)
                    subda = subda + da;
                    
                    if mode == 0
                        % convert ah to Ah
                        Ahb(:,1) = phih * ah(Ixh);
                        Ahb(:,2) = phih * ah(Iyh);
                        Ahb(:,2) = Ahb(:,2) - mean(Ahb(:,2));
                    else
                        % convert av to Av
                        Avb(:,1) = phiv * av(Ixv);
                        Avb(:,2) = phiv * av(Iyv);
                        Avb(:,1) = Avb(:,1) - mean(Avb(:,1));
                    end
                    
                    % extrapolate
                    [Ahb, Avb] = extrapolate(Ahb, Avb, LTh, LTv);
                    
                    % update the figure
                    if opt.Plot && ishandle(hdl_fig)
                        title(hdl_a(1),sprintf('Step %d/%d, it %d/%s/%d, r %10.4e',istep,Nsteps,it,HVstr(mode+1),subit,r),'FontSize',opt.fontsize);
                        set(hdl_p(1),'XData',Ahb(:,1) * Bin);
                        set(hdl_p(2),'XData',Ahb(:,2) * Bin);
                        set(hdl_p(3),'YData',Avb(:,1) * Bin);
                        set(hdl_p(4),'YData',Avb(:,2) * Bin);
                        drawnow
                    end
                    
                    if opt.Verbose > 1
                        logprintf(opt.Logfile,'   it%3d/%s/%02d, r: %10.4e, da: %10.4e\n',it,HVstr(mode+1),subit,r,rms(da));
                        appendstat(Hs,'   it%3d/%s/%02d, r: %10.4e, da: %10.4e\n',it,HVstr(mode+1),subit,r,rms(da));
                    end
                    
                    if rms(da) < ConvCrit
                        break
                    end
                    if ishandle(hdl_abort) && (hdl_abort.Value == 1)
                        appendstat(Hs,'ABORT');
                        return
                    end
                    
                end
                
                % update the reference image
                [f, R] = scan_corr_im_mex(ghb,gvb,Ahb,Avb,InterpMethod);
                
                if opt.Plot && ishandle(hdl_fig)
                    title(hdl_a(1),sprintf('Step %d/%d, it %d/%s/%d, r %10.4e',istep,Nsteps,it,HVstr(mode+1),subit,r),'FontSize',opt.fontsize);
                    AlphaData =  Iroi + 0.5 * not(Iroi);
                    set(hdl_i,'CData',R,'AlphaData',AlphaData);
                    clim = get_clim(R(Iroi),[0.01,0.99]);
                    set(hdl_a(1),'CLim',clim);
                    drawnow
                end
            end
            
            if rms(subda) < ConvCrit
                break
            end
            if ishandle(hdl_abort) && (hdl_abort.Value == 1)
                appendstat(Hs,'ABORT');
                return
            end
        end
        
        % unbin the displacements
        Ahb = Ahb * Bin;
        Avb = Avb * Bin;
        nn = Bin * nb;
        mm = Bin * mb;
        Ah(1:nn,1) = reshape(repmat(Ahb(:,1).',Bin,1),nn,1);
        Ah(1:nn,2) = reshape(repmat(Ahb(:,2).',Bin,1),nn,1);
        Av(1:mm,1) = reshape(repmat(Avb(:,1).',Bin,1),mm,1);
        Av(1:mm,2) = reshape(repmat(Avb(:,2).',Bin,1),mm,1);
        
        if opt.Verbose == 1
            logprintf(opt.Logfile,'   Step %2d/%2d, it: %3d (H%3d,V%3d) r: %10.4e\n',istep,Nsteps,sum(itcount(istep,:)),itcount(istep,:),r);
            appendstat(Hs,'   Step %2d/%2d, it: %3d (H%3d,V%3d) r: %10.4e',istep,Nsteps,sum(itcount(istep,:)),itcount(istep,:),r);
        end
        
        
    end
    
else    
    % when steps is disabled
    InterpMethod = 5;
    LineThickness = 1;
        
    roi = opt.ROI;
    roih = round(roi(1:2));
    roiv = round(roi(3:4));
    
    % create the h-blocks
    LTh = zeros(n,1);
    Nrh = roiv(2) - roiv(1);
    Nh = floor(Nrh / LineThickness);
    ind = floor(linspace(1,Nh+1,Nrh+1));
    LTh( roiv(1)+(1:Nrh) ) = ind(1:Nrh);

    % create the v-blocks
    LTv = zeros(m,1);
    Nrv = roih(2) - roih(1);
    Nv = floor(Nrv / LineThickness);
    ind = floor(linspace(1,Nv+1,Nrv+1));
    LTv( roiv(1)+(1:Nrv) ) = ind(1:Nrv);
    
end

% Output
% ==================================================

% recompute the image without binning nor blurring
[f, R] = scan_corr_im_mex(gh,gv,Ah,Av,InterpMethod);

if opt.SavePNG
    savepng(sprintf('%s_scancorr.png',opt.Savename),hdl_fig,opt.dpi)
end

out.im = f;
out.R = single(R);
out.r = r;
out.uh = Ah;
out.uv = Av;

% compute the indicators
Mh = [ones(n,1), (1:n).'];
Mv = [ones(m,1), (1:m).'];
ah = Mh(LTh>0,:) \ out.uh(LTh>0,:);
av = Mv(LTv>0,:) \ out.uv(LTv>0,:);
std_uh = std(out.uh(LTh>0,:) - Mh(LTh>0,:) * ah);
std_uv = std(out.uv(LTv>0,:) - Mv(LTv>0,:) * av);

out.ind_drift = [ah(2,:); av(2,:)];
out.ind_linestd = [std_uh; std_uv];


if opt.Verbose > 0
    logprintf(opt.Logfile,'--------------------------------------\n');
    logprintf(opt.Logfile,'Correction Summary:\n');
    logprintf(opt.Logfile,'   Iterations   : total: %d, H %d, V %d\n',sum(itcount(:)),sum(itcount,1));
    logprintf(opt.Logfile,'   Residual     : %10.4e\n',r);
    logprintf(opt.Logfile,'   Displacement : hx:%10.3e, hy:%10.3e, vx:%10.3e, vy:%10.3e \n',rms(Ah),rms(Av));
    logprintf(opt.Logfile,'   Drift        : hx:%10.3e, hy:%10.3e, vx:%10.3e, vy:%10.3e\n',out.ind_drift([1,3,2,4]));
    logprintf(opt.Logfile,'   Line STD     : hx:%10.3e, hy:%10.3e, vx:%10.3e, vy:%10.3e\n',out.ind_linestd([1,3,2,4]));
    logprintf(opt.Logfile,'--------------------------------------\n');
    appendstat(Hs,'Correction Summary:');
    appendstat(Hs,'   Iterations   : total: %d, H %d, V %d',sum(itcount(:)),sum(itcount,1));
    appendstat(Hs,'   Residual     : %10.4e\n',r);
    appendstat(Hs,'   Displacement : hx:%10.3e, hy:%10.3e, vx:%10.3e, vy:%10.3e',rms(Ah),rms(Av));
    appendstat(Hs,'   Drift        : hx:%10.3e, hy:%10.3e, vx:%10.3e, vy:%10.3e',out.ind_drift([1,3,2,4]));
    appendstat(Hs,'   Line STD     : hx:%10.3e, hy:%10.3e, vx:%10.3e, vy:%10.3e',out.ind_linestd([1,3,2,4]));
end

if opt.SaveCOR
    im = out.im;
    im(isnan(im)) = 0;
    if opt.TIFdepth == 16
        im = uint16(65535 * im);
    else
        im = uint8(255 * im);
    end
    imwrite(im,sprintf('%s_corrected.tif',opt.Savename),'TIFF','Compression','none');
end

if opt.SaveRES
    res = out.R;
    res(isnan(res)) = 0;
    if opt.TIFdepth == 16
        res = uint16(10 * 65535 * res + 32768);
    else
        res = uint8(10 * 255 * res + 128);
    end
    imwrite(res,sprintf('%s_residual.tif',opt.Savename),'TIFF','Compression','none');
end

if opt.ROI_NAN
    [X, Y] = meshgrid(1:m,1:n);
    roi = opt.ROI;
    Iroi = X >= roi(1) & X <= roi(2) & Y >= roi(3) & Y <= roi(4);
    roix = 1:m >= roi(1) & 1:m <= roi(2);
    roiy = 1:n >= roi(3) & 1:n <= roi(4);
    
    out.im(~Iroi) = NaN;
    out.R(~Iroi) = NaN;
    out.uh(~roiy,:) = NaN;
    out.uv(~roix,:) = NaN;
end

if opt.RBT_flag
    out.urbt = [urbt, C];
end

end

function phi = BlockMap(LT,N)
% mapping from pixels to blocks
n = numel(LT);
J = LT;
I = 1:n;
phi = sparse(I(LT~=0),J(LT~=0),1,n,N);
end

function S = SensitivityPoly(x,N,d)
% mapping to a polygon basis
if d <= 0
    S = [];
    return
end
Ix = 1:2:2*N;
Iy = 2:2:2*N;
S = zeros(2*N,2*(d+1));
for k = 0:d
    S(Ix,2*k+1) = x.^k;
    S(Iy,2*k+2) = x.^k;
end
end

function da = scan_corr_solve(M,b,S)
N = numel(b)/2;
if isempty(S)
    % no reprojection, direct solution
    I = repmat([1, 2; 3, 4], N, 1);
    M11 = M(I == 1);
    M12 = M(I == 2);
    M21 = M(I == 3);
    M22 = M(I == 4);
    
    b1 = b(1:2:2*N,1);
    b2 = b(2:2:2*N,1);
    
    D = M11.*M22 - M12.*M21;
    
    da = zeros(2*N,1);
    da(1:2:2*N,1) = M22.*b1./D - M12.*b2./D;
    da(2:2:2*N,1) = M11.*b2./D - M21.*b1./D;
else
    % formatting M in a large sparse matrix
    I = [(1:2*N), (1:2*N)].';
    J = [(1:2:2*N), (2:2:2*N) ; (1:2:2*N), (2:2:2*N)];
    J = reshape(J,2*N,2);
    M = sparse(I,J,M,2*N,2*N);
    
    % projecting to the new bases
    M = transpose(S) * M * S;
    b = transpose(S) * b;
    
    % sovling
    dp = M \ b;
    
    % reprojecting back to the old basis
    da = S * dp;
end

da(isnan(da)) = 0;
end

function [Ah, Av] = extrapolate(Ah, Av, binh, binv)
n = size(Ah,1);
m = size(Av,1);

% extrapolate the displacements
Ni = round(0.1*min(n,m));

% h
x = linspace(-1,1,n).';
M = [ones(n,1), x];

I1 = find(binh > 0,Ni,'first');
I2 = 1:I1(1)-1;
I3 = find(binh > 0,Ni,'last');
I4 = I3(end)+1 : n;

a = M(I1,:) \ Ah(I1,:);
Ah(I2,:) = M(I2,:) * a;

a = M(I3,:) \ Ah(I3,:);
Ah(I4,:) = M(I4,:) * a;

% v
x = linspace(-1,1,m).';
M = [ones(m,1), x];

I1 = find(binv > 0,Ni,'first');
I2 = 1:I1(1)-1;
I3 = find(binv > 0,Ni,'last');
I4 = I3(end)+1 : m;

a = M(I1,:) \ Av(I1,:);
Av(I2,:) = M(I2,:) * a;

a = M(I3,:) \ Av(I3,:);
Av(I4,:) = M(I4,:) * a;
end

function A = databar_remove(A,T)
% the threshold is the relative number of pixels that have zero gradient,
% if a pixel row has more than T, it and all rows below will be conidered
% databar.

m = size(A,2);

% count the number of pixels with zero derivative
d = sum(diff(A,1,2)==0,2);

% find the first row that surpasses the threshold
n = find(d > T*m,1) - 1;

if isempty(n)
    return
elseif n < 0.75*size(A,1)
    warning('databar removal failed');
else
    % crop out the databar
    A = A(1:n,:);
end

end

% ==================================================
function x = clamp(x,a,b)
x = max(x,a);
x = min(x,b);
end

% ==================================================
function A = imcoarsgrain(A,sps)
if any(sps < 1) || all(sps <= 1)
    return
end
if numel(sps) == 1
    sps = [sps,sps];
end
if isinteger(A)
    A = double(A);
end

% coarsegrain in each direction
siz = size(A);
for i = 1:2
    
    I = 1 : sps(i) : floor(siz(i)/sps(i))*sps(i);
    siz(i) = length(I);
    B = zeros(siz);
    for j = 1:sps(i)
        switch i
            case 1; B = B + A(I,:,:);
            case 2; B = B + A(:,I,:);
        end
        
        % go the the next subpixel in this row
        I = I + 1;
    end
    % overwrite A
    A = B / sps(i);
end
end

% ==================================================
function clim = get_clim(A,qtl)
clim = quantile(A(:),qtl);
if (clim(1) == clim(2))
    if clim(1) ~= 0
        clim = clim(1)*[-0.01,0.01];
    elseif clim(2) ~= 0
        clim = clim(2)*[-0.01,0.01];
    else
        clim = [-1,1];
    end
elseif  any(isnan(clim)) || any(isinf(clim))
    clim = [-1,1];
else
    clim = max(abs(clim))*[-1,1];
end
end

% ==================================================
function savepng(filename,varargin)
if nargin == 1
    H = gcf;
    r = 200;
elseif nargin == 2
    H = varargin{1};
    r = 200;
elseif nargin == 3
    H = varargin{1};
    r = varargin{2};
else
    error('savepng: incorrect number of inputs');
end

C = get(H,'Color');
set(H,'Color',[1 1 1]);


rstr = sprintf('-r%d',r);

% fix the extention, to be always .png (small case)
filename = [regexprep(filename,'.png$','','ignorecase') '.png'];

% save pdf
print(H,filename,'-dpng',rstr)

set(H,'Color',C);
end

% ==================================================
function logprintf(logfile,varargin)
% logprint(logfile,...), logprint is a wrapper around fprintf, instead of a
% file identifier (fid) specify the name of a file in logfile. logprint
% will open the file, print the message both to the file and to the screen
% and then close the file. If the file already exists, the message will be
% appended.
fprintf(varargin{:});
if ~isempty(logfile)
    fid = fopen(logfile,'a+t');
    fprintf(fid,varargin{:});
    fclose(fid);
end
end

% -------------------------------
function appendstat(Hs,varargin)
if ishandle(Hs)
    stat = [sprintf(varargin{:}) ; Hs.String];
    Hs.String = stat;
    drawnow;
end
end

function cmap = RdGy()
cmap = [        0.1020    0.1020    0.1020
    0.1420    0.1420    0.1420
    0.1820    0.1820    0.1820
    0.2220    0.2220    0.2220
    0.2620    0.2620    0.2620
    0.3020    0.3020    0.3020
    0.3475    0.3475    0.3475
    0.3929    0.3929    0.3929
    0.4384    0.4384    0.4384
    0.4839    0.4839    0.4839
    0.5294    0.5294    0.5294
    0.5694    0.5694    0.5694
    0.6094    0.6094    0.6094
    0.6494    0.6494    0.6494
    0.6894    0.6894    0.6894
    0.7294    0.7294    0.7294
    0.7592    0.7592    0.7592
    0.7890    0.7890    0.7890
    0.8188    0.8188    0.8188
    0.8486    0.8486    0.8486
    0.8784    0.8784    0.8784
    0.9027    0.9027    0.9027
    0.9271    0.9271    0.9271
    0.9514    0.9514    0.9514
    0.9757    0.9757    0.9757
    1.0000    1.0000    1.0000
    0.9984    0.9718    0.9561
    0.9969    0.9435    0.9122
    0.9953    0.9153    0.8682
    0.9937    0.8871    0.8243
    0.9922    0.8588    0.7804
    0.9851    0.8165    0.7263
    0.9780    0.7741    0.6722
    0.9710    0.7318    0.6180
    0.9639    0.6894    0.5639
    0.9569    0.6471    0.5098
    0.9333    0.5929    0.4682
    0.9098    0.5388    0.4267
    0.8863    0.4847    0.3851
    0.8627    0.4306    0.3435
    0.8392    0.3765    0.3020
    0.8110    0.3200    0.2753
    0.7827    0.2635    0.2486
    0.7545    0.2071    0.2220
    0.7263    0.1506    0.1953
    0.6980    0.0941    0.1686
    0.6392    0.0753    0.1592
    0.5804    0.0565    0.1498
    0.5216    0.0376    0.1404
    0.4627    0.0188    0.1310
    0.4039         0    0.1216];
end


% ---------------------------------------------
function A = grayscale(A)
if ismatrix(A)
    return
end
if (ndims(A) == 3) && (size(A,3) == 3)
    if isinteger(A)
        A = single(A);
    end
    A = 0.2989*A(:,:,1) + 0.5870*A(:,:,2) + 0.1140*A(:,:,3);
end
end

% ---------------------------------------------
function K = DiscreteGaussKernel(s,n)
% K = DiscreteGaussKernel(s), compute the discrete gauss kernel at
% integer locations using the gamma function
%
% K = DiscreteGaussKernel(s,n), specify the width of the kernel if left
% empty n = ceil(6*s)


if (nargin == 1) || isempty(n)
    n = ceil(6*s);
else
    n = round(n);
end

% create the interpolation space
x = -n:n;

% define the left and right sides of the pixels
x1 = x - 0.5;
x2 = x + 0.5;

if false
    % integrate from -inf to inf (add the values to the ends of the kernel)
    N = 2 * n + 1;
    x1(1) = -inf;
    x2(N) = +inf;
end

% integrate the left and right side of the pixel
K1 = 0.5 * erf(sqrt(2) * x1 / ( 2 * s ));
K2 = 0.5 * erf(sqrt(2) * x2 / ( 2 * s ));

% compute the integral over the pixel
K = K2 - K1;

% make sure the kernel summs to one (it may be less because the tails are
% not included in the integration)
K = K ./ sum(K);
end


% ---------------------------------------------
function A = imageblur(A,R)
% B = image_blur(A, R) blurs the image A with a gaussian kernel of radius R

if all(R < 0.1)
    return
end
if numel(R) == 1
    R = [R, R];
end


for i = 1:2
    
    % create the kernel
    % ------------------------
    K = DiscreteGaussKernel(R(i),4);
    
    % mirroring the image
    % ------------------------
    [am,an,ap] = size(A);
    km         = length(K);
    
    % mirror distance
    Nm = floor(km/2);
    
    % skip a direction if it has too few pixels
    if (Nm >= am) || (R(i) == 0)
        if ismatrix(A)
            % transpose the image and blur again
            A = A';
        elseif ndims(A) == 3
            % transpose the image and blur again
            A = permute(A,[2 3 1]);
        end
        continue
    end
    
    Amir = zeros(am+2*Nm,an,ap,'single');
    
    % image rows and columns in mirrored image
    Im = Nm+1:Nm+am;
    
    % image rows and columns in the original image
    Jm = 1:am;
    
    % mirror parts in mirrored image
    Iml = Nm:-1:1; % left columns
    Imr = 2*Nm+am:-1:Nm+1+am; % right columns
    
    % mirror parts in original image
    Jml = 2:Nm+1;
    Jmr = (am:am+Nm-1)-Nm;
    
    % Create an image with mirrored edges
    Amir([Iml,Im,Imr],:,:) = A([Jml,Jm,Jmr],:,:);
    
    % Apply the blur using convolution
    % ------------------------
    
    A = conv2(Amir,K,'valid');
    
    % transpose the image and blur again
    A = A';
end

end
