function varargout = scan_corr_GUI()

% lists of stuff
tiffdepths = {'8';'16'};
InterpMethods = {'bspline2','bspline3','bspline4','bspline5','bspline6','bspline7','bspline8','bspline9'};
PolyOrders = {'None','linear','quadratic','cubic','quartic','quintic'};

if isdeployed
    inifile = fullfile(userpath,'scan_corr_GUI.ini');
    currentpath = userpath;
else
    inifile = 'scan_corr_GUI.ini';
    currentpath = pwd;
end

% default options
defini.TIFdepth = 1; % 1:8 or 2:16
defini.SaveLOG = true;
defini.SaveCOR = true;
defini.SaveRES = true;
defini.SavePNG = true;
defini.SaveMAT = true;
defini.SaveRBT = true;

defini.RotateH = 0; % number of 90 degrees rotations
defini.RotateV = 1; % number of 90 degrees rotations
defini.DatabarRemove = true;
defini.ROI = 50; % pixels from the boundary
defini.ROI_NAN = false;

defini.RBT_flag = true;
defini.RBT_Window = 0.1; % tukeywin

defini.Steps_flag = true;
defini.LineThickness = [1, 1, 10, 3, 1];
defini.PolyOrder = [2, 2, 0];
defini.Bin = [4, 4, 2, 2, 1];
defini.Blur = [8, 4, 2, 1, 0];
defini.ConvCrit = 10.^(-1*[2, 2, 2, 2, 4]);
defini.IterMax = 30;
defini.SubIterMax = 3;
defini.TrustRegion = 2;
defini.InterpMethod = [3, 3, 3, 3, 5];

% create the ini-file if it doens't exit
if ~exist(inifile,'file')
    iniwrite(defini);
end

% read ini-file, with the stored defaults
ini = iniread();

% the step options
step_opts{1,1} = 'PolyOrder';
step_opts{2,1} = 'LineThickness';
step_opts{3,1} = 'Bin';
step_opts{4,1} = 'Blur';
step_opts{5,1} = 'ConvCrit';
step_opts{6,1} = 'IterMax';
step_opts{7,1} = 'SubIterMax';
step_opts{8,1} = 'TrustRegion';
step_opts{9,1} = 'InterpMethod';

% the the longest definition
Nsteps = 1;
for ii = 1:numel(step_opts)
    N = numel(ini.(step_opts{ii}));
    Nsteps = max(Nsteps,N);
end

% pad all using the last value
for ii = 1:numel(step_opts)
    N = numel(ini.(step_opts{ii}));
    ini.(step_opts{ii})(N+1:Nsteps) = ini.(step_opts{ii})(N);
end

SO = [];
StepList = {};
steplist_update();

ini.Bin = round(ini.Bin);
ini.LineThickness = round(ini.LineThickness);
ini.IterMax = round(ini.IterMax);
ini.SubIterMax = round(ini.SubIterMax);
ini.InterpMethod = round(ini.InterpMethod);
if any( ini.InterpMethod < 2 ) || any( ini.InterpMethod > 9 )
    error( 'InterpMethod must be between 2 and 9 (inclusive)' )
end


% colors
color.bg = [1, 1, 1];
color.button = 0.95*[1, 1, 1];
color.section = 0.95*[1, 1, 1];
color.go = [0.533 0.80 0.533];
color.abort = [0.831 0.416 0.416];
color.progress1 = [0.533 0.80 0.533];
color.progress2 = [0.333 0.667 0.333];
color.fontA = '#804515';
color.fontB = '#0D4D4D';
color.col1 = [0.133 0.40 0.40];
color.col2 = [0.667 0.424 0.224];
color.roi = [0.0 1.0 1.0];
color.cent = [1.0 1.0 0.0];


cmap = RdGy();

% Globals
% ===============================================
savepath = '';
filelist = {};
rule = repmat('=',1,40);
v = [];

% possible image types
imagetypes = {'*.png';'*.jpg';'*.jpeg';'*.gif';'*.tif';'*.tiff';'*.bmp'};
% add all caps versions of  these files
imagetypes = [imagetypes ; upper(imagetypes)];
% convert to an array as a string
imtypstr = sprintf(['%s',repmat(';%s',1,numel(imagetypes)-1)],imagetypes{:});

% Create the Figure
% ===============================================

% get the monitor size and the active monitor
Hr = groot;
mon = Hr.MonitorPositions;
Nmon = size(mon,1);
if Nmon > 1
    p0 = Hr.PointerLocation;
    Imon = find(p0(1) >= mon(:,1) & p0(1) < mon(:,1)+mon(:,3) & p0(2) >= mon(:,2) & p0(2) < mon(:,2)+mon(:,4));
    if ~isempty(Imon)
        mon = mon(Imon,:);
    else
        mon = mon(1,:);
    end
end

% determine the fontsize size
if mon(4) <= 800
    fontsize = [7,8,16];
    screenfillfactor = 1;
elseif mon(4) <= 1024
    fontsize = [8,9,18];
    screenfillfactor = 0.95;
elseif mon(4) <= 1280
    fontsize = [10,12,20];
    screenfillfactor = 0.85;
else
    fontsize = [12,14,22];
    screenfillfactor = 0.75;
end

if ispc
    fixwidthfont = 'FixedWidth';
    % fixwidthfont = 'Courier New';
else
    fixwidthfont = 'Monospaced';
end


% figure positions
height = screenfillfactor*mon(4);
width = screenfillfactor*mon(3);
left = 0.5*(mon(3)-width);
bottom = 0.5*(mon(4)-height);

% image window position
pos.fig = [left, bottom, width, height];


% The toolbar figure
Hf = figure('OuterPosition',pos.fig);
Hf.Name = 'Scan_Corr 2.01';
Hf.NumberTitle = 'off';
Hf.Color = color.bg;
Hf.MenuBar = 'none';
Hf.ToolBar = 'figure';
Hf.Tag = 'ScanCorr';

Hf.PaperUnits = 'inches';
Hf.PaperPosition = pos.fig.*[0 0 1e-2 1e-2];
Hf.PaperSize = pos.fig(3:4).*[1e-2 1e-2];

% Create the TABS
tabgp = uitabgroup(Hf,'Position',[0 0 1 1]);
tab(1) = uitab(tabgp,'Title','Preparation','BackgroundColor',color.bg);
tab(2) = uitab(tabgp,'Title','Correction','BackgroundColor',color.bg);

% populate the first tab
% ===============================================

% Create the GUI panels
% ----------------------------------
Hp1 = uipanel('Title','Images','Parent',tab(1),...
    'Position',[0/3,0,1/3,1],'BackgroundColor',color.bg);
Hp2 = uipanel('Title','Settings','Parent',tab(1),...
    'Position',[1/3,0,1/3,1],'BackgroundColor',color.bg);
Hp3 = uipanel('Title','Help','Parent',tab(1),...
    'Position',[2/3,0,1/3,1],'BackgroundColor',color.bg);

% Image Controls
x0 = linspace(0.01,0.99,4);
y0 = linspace(0.01,0.99,30);
dy = 0.98*(y0(2)-y0(1));
Hb.imlist = uicontrol('String','',...
    'Style','listbox',...
    'Units','normalized',...
    'Value',[],...
    'Max',2,...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(1),x0(4)-x0(1),y0(15)-y0(1)],...
    'Parent',Hp1,...
    'call',{@image_fun,'select'});

uicontrol('String','Add',...
    'ToolTipString','add images to the list',...
    'Style','pushbutton',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(16),x0(2)-x0(1),dy],...
    'Parent',Hp1,...
    'BackgroundColor',color.button,...
    'call',{@image_fun,'add'});

uicontrol('String','Remove',...
    'ToolTipString','remove the selected images from the list',...
    'Style','pushbutton',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(15),x0(2)-x0(1),dy],...
    'Parent',Hp1,...
    'BackgroundColor',color.button,...
    'call',{@image_fun,'del'});

uicontrol('String','Up',...
    'ToolTipString','move the selected images up',...
    'Style','pushbutton',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(2),y0(16),x0(2)-x0(1),dy],...
    'Parent',Hp1,...
    'BackgroundColor',color.button,...
    'call',{@image_fun,'up'});

uicontrol('String','Down',...
    'ToolTipString','move the selected images down',...
    'Style','pushbutton',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(2),y0(15),x0(2)-x0(1),dy],...
    'Parent',Hp1,...
    'BackgroundColor',color.button,...
    'call',{@image_fun,'down'});

uicontrol('String','Interlace',...
    'ToolTipString','Iterlace the selection [a,a,a,b,b,b] -> [a,b,a,b,a,b]',...
    'Style','pushbutton',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(3),y0(16),x0(2)-x0(1),dy],...
    'Parent',Hp1,...
    'BackgroundColor',color.button,...
    'call',{@image_fun,'interlace'});

uicontrol('String','Flip Odd/Even',...
    'ToolTipString','Interchange the odd and even files in the selection [a,b,a,b,a,b] -> [b,a,b,a,b,a]',...
    'Style','pushbutton',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(3),y0(15),x0(2)-x0(1),dy],...
    'Parent',Hp1,...
    'BackgroundColor',color.button,...
    'call',{@image_fun,'flip'});


Ha.preview = axes('Position',[x0(1),y0(17),x0(4)-x0(1),y0(end)-y0(17)],'Parent',Hp1);
Hi.preview = imagesc(NaN,'Parent',Ha.preview);
set(Ha.preview,...
    'Box','On',...
    'XTick',[],...
    'YTick',[],...
    'XColor','none',...
    'YColor','none',...
    'DataAspectRatio',[1 1 1]);
colormap(Ha.preview,gray);
Hp.roi = patch('Vertices',nan(4,2),'Faces',1:4,'EdgeColor',color.roi,'FaceColor','None','Parent',Ha.preview);
Hp.cent = patch('Vertices',nan(4,2),'Faces',1:4,'EdgeColor',color.cent,'FaceColor','None','Parent',Ha.preview);


% options
% ------------------------------
x0 = linspace(0.01,0.99,3);
dx = 0.98*(x0(2)-x0(1));
y0 = linspace(0.99,0.01,35);
dy = 0.98*(y0(1)-y0(2));

% create the controls
row = 2;
uicontrol('String','Save path:',...
    'Style','text',...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(row),dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
uicontrol('String','Browse',...
    'ToolTipString','Select a save path',...
    'Style','pushbutton',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(2),y0(row),dx,dy],...
    'Parent',Hp2,...
    'BackgroundColor',color.button,...
    'call',@set_save_path);
row = row + 1;
Hb.savepath = uicontrol('String',savepath,...
    'ToolTipString','Choose a path where to save the results',...
    'Style','edit',...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(row),x0(3)-x0(1),dy],...
    'Parent',Hp2);
row = row + 1;

Hb.SaveLOG = uicontrol('String','Save logfile (.txt)',...
    'ToolTipString','Save the status output to a .txt file',...
    'Style','checkbox',...
    'Value',ini.SaveLOG,...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(row),dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
Hb.SavePNG = uicontrol('String','Save final screen (.png)',...
    'ToolTipString','Save the final screen as a .png (may cause flicker)',...
    'Style','checkbox',...
    'Value',ini.SavePNG,...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(2),y0(row),dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
row = row + 1;
Hb.SaveCOR = uicontrol('String','Save corrected image (.tif)',...
    'ToolTipString','Save the corrected image as a .tif',...
    'Style','checkbox',...
    'Value',ini.SaveCOR,...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(row),dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
Hb.SaveMAT = uicontrol('String','Save data (.mat)',...
    'ToolTipString','Save the output as a matlab .mat data file',...
    'Style','checkbox',...
    'Value',ini.SaveMAT,...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(2),y0(row),dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
row = row + 1;
Hb.SaveRES = uicontrol('String','Save residual image (.tif)',...
    'ToolTipString','Save the residual as a .tif (amplified with a factor 10)',...
    'Style','checkbox',...
    'Value',ini.SaveRES,...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(row),dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
Hb.SaveRBT = uicontrol('String','Save RBT images (.tif)',...
    'ToolTipString','Save the cropped image after RBT corrections',...
    'Style','checkbox',...
    'Value',ini.SaveRBT,...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(2),y0(row),dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
row = row + 1;
uicontrol('String','TIF bit depth',...
    'Style','text',...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(row),dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
Hb.TIFdepth = uicontrol('String',{'8';'16'},...
    'ToolTipString','Select a bit depth at which to store the corrected image',...
    'Style','popupmenu',...
    'Value',ini.TIFdepth,...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(2),y0(row),dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);

% -------------------
row = row + 2;
uicontrol('String','Preparation',...
    'Style','text',...
    'Value',ini.RBT_flag,...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontWeight','bold',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(row),x0(3)-x0(1),dy],...
    'BackgroundColor',color.section,...
    'Parent',Hp2);
row = row + 1;
uicontrol('String','Rotate H',...
    'Style','text',...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(row),dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
Hb.RotateH = uicontrol('String',num2str(ini.RotateH),...
    'ToolTipString','How many times to rotate the H image 90 degrees counter clockwise',...
    'Style','Edit',...
    'HorizontalAlignment','center',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(2),y0(row),dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
row = row + 1;
uicontrol('String','Rotate V',...
    'Style','text',...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(row),dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
Hb.RotateV = uicontrol('String',num2str(ini.RotateV),...
    'ToolTipString','How many times to rotate the V image 90 degrees counter clockwise',...
    'Style','Edit',...
    'HorizontalAlignment','center',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(2),y0(row),dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
row = row + 1;
Hb.DatabarRemove = uicontrol('String','Auto remove databar',...
    'ToolTipString','Automatically detect and remove the databar',...
    'Style','checkbox',...
    'Value',ini.DatabarRemove,...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(row),x0(3)-x0(1),dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
row = row + 1;
uicontrol('String','ROI',...
    'Style','text',...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(row),dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
Hb.ROI = uicontrol('String',mat2str(ini.ROI),...
    'ToolTipString','The size of the border (in px) to exclude from the ROI, alternativly use [left, right, top, bottom]',...
    'Style','Edit',...
    'HorizontalAlignment','center',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(2),y0(row),0.5*dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
uicontrol('String','Select ROI',...
    'ToolTipString','Select the ROI on the image on the left',...
    'Style','PushButton',...
    'HorizontalAlignment','center',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(2)+0.5*dx,y0(row),0.5*dx,dy],...
    'call',@roi_select,...
    'BackgroundColor',color.button,...
    'Parent',Hp2);

row = row + 2;
Hb.RBT_flag = uicontrol('String','Rigid Body Translation',...
    'Style','checkbox',...
    'Value',ini.RBT_flag,...
    'ToolTipString','uncheck to disable this part',...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontWeight','bold',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(row),x0(3)-x0(1),dy],...
    'BackgroundColor',color.section,...
    'Parent',Hp2);
row = row + 1;
uicontrol('String','RBT Window',...
    'Style','text',...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(row),dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
Hb.RBT_Window = uicontrol('String',num2str(ini.RBT_Window),...
    'ToolTipString','The relative border size of the window (tukeywin) to apply for the FFT for the RBT step (set to 0 to disable)',...
    'Style','Edit',...
    'HorizontalAlignment','center',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(2),y0(row),dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);

% -----------------------------------------------------------------
row = row + 2;
Hb.Steps_flag = uicontrol('String','Pyramid Steps',...
    'Style','checkbox',...
    'Value',ini.Steps_flag,...
    'ToolTipString','uncheck to disable this part',...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontWeight','bold',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(row),x0(3)-x0(1),dy],...
    'BackgroundColor',color.section,...
    'Parent',Hp2);

row = row + 5;
Hb.StepList = uicontrol('String',StepList,...
    'Style','listbox',...
    'Units','normalized',...
    'Value',1,...
    'Max',2,...
    'FontName',fixwidthfont,...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(row),x0(3)-x0(1),y0(1)-y0(6)],...
    'Parent',Hp2,...
    'call',{@step_fun,'select'});
row = row + 1;
uicontrol('String','Add',...
    'ToolTipString','add a step to the list',...
    'Style','pushbutton',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1)+0*(x0(3)-x0(1))/5,y0(row),(x0(3)-x0(1))/5,dy],...
    'Parent',Hp2,...
    'BackgroundColor',color.button,...
    'call',{@step_fun,'add'});
uicontrol('String','Copy',...
    'ToolTipString','duplicate the currently selected steps',...
    'Style','pushbutton',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1)+1*(x0(3)-x0(1))/5,y0(row),(x0(3)-x0(1))/5,dy],...
    'Parent',Hp2,...
    'BackgroundColor',color.button,...
    'call',{@step_fun,'copy'});
uicontrol('String','Remove',...
    'ToolTipString','remove the selected steps from the list',...
    'Style','pushbutton',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1)+2*(x0(3)-x0(1))/5,y0(row),(x0(3)-x0(1))/5,dy],...
    'Parent',Hp2,...
    'BackgroundColor',color.button,...
    'call',{@step_fun,'del'});
uicontrol('String','Up',...
    'ToolTipString','remove the selected steps from the list',...
    'Style','pushbutton',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1)+3*(x0(3)-x0(1))/5,y0(row),(x0(3)-x0(1))/5,dy],...
    'Parent',Hp2,...
    'BackgroundColor',color.button,...
    'call',{@step_fun,'up'});
uicontrol('String','Down',...
    'ToolTipString','remove the selected steps from the list',...
    'Style','pushbutton',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1)+4*(x0(3)-x0(1))/5,y0(row),(x0(3)-x0(1))/5,dy],...
    'Parent',Hp2,...
    'BackgroundColor',color.button,...
    'call',{@step_fun,'down'});

row = row + 1;
uicontrol('String','PolyOrder',...
    'Style','text',...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(row),dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
Hb.PolyOrder = uicontrol('String',PolyOrders,...
    'ToolTipString','The order of the regularizing polynomial basis, use None to disable, LineThickness will be forced to 1 if enabled',...
    'Style','popupmenu',...
    'Value',ini.PolyOrder(1)+1,...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(2),y0(row),dx,dy],...
    'call',{@step_fun,'opt'},...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
row = row + 1;
uicontrol('String','LineThickness',...
    'Style','text',...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(row),dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
Hb.LineThickness = uicontrol('String',mat2str(ini.LineThickness(1)),...
    'ToolTipString','Group N lines together and force them to move as one',...
    'Style','Edit',...
    'HorizontalAlignment','center',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(2),y0(row),dx,dy],...
    'call',{@step_fun,'opt'},...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
row = row + 1;
uicontrol('String','Bin',...
    'Style','text',...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(row),dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
Hb.Bin = uicontrol('String',mat2str(ini.Bin(1)),...
    'ToolTipString','Bin the image using superpixels of size NxN (this increases speed and robustness at the cost of accuracy)',...
    'Style','Edit',...
    'HorizontalAlignment','center',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(2),y0(row),dx,dy],...
    'call',{@step_fun,'opt'},...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
row = row + 1;
uicontrol('String','Blur',...
    'Style','text',...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(row),dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
Hb.Blur = uicontrol('String',mat2str(ini.Blur(1)),...
    'ToolTipString','Blur the image with a radius of R (this increases robustness at the cost of accuracy)',...
    'Style','Edit',...
    'HorizontalAlignment','center',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(2),y0(row),dx,dy],...
    'call',{@step_fun,'opt'},...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
row = row + 1;
uicontrol('String','ConvCrit',...
    'Style','text',...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(row),dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
Hb.ConvCrit = uicontrol('String',mat2str(ini.ConvCrit(1)),...
    'ToolTipString','Consider the step converged if the iterative update (da) is smaller then X',...
    'Style','Edit',...
    'HorizontalAlignment','center',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(2),y0(row),dx,dy],...
    'call',{@step_fun,'opt'},...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
row = row + 1;
uicontrol('String','IterMax',...
    'Style','text',...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(row),dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
Hb.IterMax = uicontrol('String',mat2str(ini.IterMax(1)),...
    'ToolTipString','Maximum number of iterations per step (for both H and V)',...
    'Style','Edit',...
    'HorizontalAlignment','center',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(2),y0(row),dx,dy],...
    'call',{@step_fun,'opt'},...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
row = row + 1;
uicontrol('String','SubIterMax',...
    'Style','text',...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(row),dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
Hb.SubIterMax = uicontrol('String',mat2str(ini.SubIterMax(1)),...
    'ToolTipString','Maximum number of sub-iterations (within one H or V loop)',...
    'Style','Edit',...
    'HorizontalAlignment','center',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(2),y0(row),dx,dy],...
    'call',{@step_fun,'opt'},...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
row = row + 1;
uicontrol('String','TrustRegion',...
    'Style','text',...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(row),dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
Hb.TrustRegion = uicontrol('String',mat2str(ini.TrustRegion(1)),...
    'ToolTipString','Limit the iterative update to X pixels',...
    'Style','Edit',...
    'HorizontalAlignment','center',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(2),y0(row),dx,dy],...
    'call',{@step_fun,'opt'},...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
row = row + 1;
uicontrol('String','InterpMethod',...
    'Style','text',...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(row),dx,dy],...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);
Hb.InterpMethod = uicontrol('String',InterpMethods,...
    'ToolTipString','The BSpline interpolation order',...
    'Value',ini.InterpMethod(1)-1,...
    'Style','popupmenu',...
    'HorizontalAlignment','center',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(2),y0(row),dx,dy],...
    'call',{@step_fun,'opt'},...
    'BackgroundColor',color.bg,...
    'Parent',Hp2);

row = row + 2;
uicontrol('String','Save Defaults',...
    'ToolTipString','Save the current settings as defaults',...
    'Style','pushbutton',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(row),dx,dy],...
    'Parent',Hp2,...
    'BackgroundColor',color.button,...
    'call',{@defaults_fun,'save'});
uicontrol('String','Reset Defaults',...
    'ToolTipString','Reset the defaults',...
    'Style','pushbutton',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(2),y0(row),dx,dy],...
    'Parent',Hp2,...
    'BackgroundColor',color.button,...
    'call',{@defaults_fun,'reset'});


Hhelp = uicontrol('Style','edit',...
    'Max',2,...
    'String','Help',...
    'Enable','inactive',...
    'units','normalized',...
    'BackgroundColor','w',...
    'HorizontalAlignment','left',...
    'Position',[0.01,0.01,0.98,0.98],...
    'FontSize',fontsize(2),...
    'Tag','help',...
    'Parent',Hp3);
help_fun();

% populate the second tab
% ===============================================

fspect = pos.fig(4)./pos.fig(3);

Hp1 = uipanel('Title','Live Result','Parent',tab(2),...
    'Position',[0,0,fspect,1],'BackgroundColor',color.bg);
Hp2 = uipanel('Title','Status','Parent',tab(2),...
    'Position',[fspect,0,1-fspect,1],'BackgroundColor',color.bg);

x0 = linspace(0.01,0.99,5);
y0 = linspace(0.01,0.99,30);
dy = 0.98*(y0(2)-y0(1));
Hs = uicontrol('String',{'';'';'Status text, new lines appear at the top, read this from the bottom upwards'},...
    'Max',2,...
    'Style','edit',...
    'HorizontalAlignment','left',...
    'Enable','inactive',...
    'Units','normalized',...
    'BackgroundColor','w',...
    'Position',[x0(1),y0(2),x0(5)-x0(1),y0(20)-y0(2)],...
    'FontName',fixwidthfont,...
    'FontSize',fontsize(2),...
    'Parent',Hp2);

Hb.go = uicontrol('String','Go',...
    'ToolTipString','add images to the list',...
    'Style','pushbutton',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(1),y0(20),x0(3)-x0(1),dy],...
    'Parent',Hp2,...
    'BackgroundColor',color.go,...
    'call',@go);

Hb.abort = uicontrol('String','Abort',...
    'ToolTipString','remove the selected images from the list',...
    'Style','togglebutton',...
    'Units','normalized',...
    'FontSize',fontsize(2),...
    'Position',[x0(3),y0(20),x0(3)-x0(1),dy],...
    'Parent',Hp2,...
    'BackgroundColor',color.abort,...
    'call',@abort);

Ha.liveview = axes('Position',[x0(1),y0(21),x0(5)-x0(1),y0(end)-y0(21)],'Parent',Hp2,'NextPlot','Add','YDir','Reverse');
Hi.liveview = imagesc(NaN,'Parent',Ha.liveview);
set(Ha.liveview,...
    'Box','On',...
    'XTick',[],...
    'YTick',[],...
    'XColor','none',...
    'YColor','none',...
    'DataAspectRatio',[1 1 1]);
colormap(Ha.liveview,gray);

Ha.progress = axes('Position',[x0(1),y0(1),x0(5)-x0(1),dy],'Parent',Hp2,'NextPlot','Add','YDir','Reverse');
Hi.progress = imagesc(ones(1,1,3),'Parent',Ha.progress(1));
set(Ha.progress,...
    'Box','On',...
    'XTick',[],...
    'YTick',[],...
    'XColor','none',...
    'YColor','none',...
    'Color','none',...
    'CLim',[0,1]);
Ht.progress = text(0.01,0.5,'Progress bar...','Units','Normalize','HorizontalAlignment','Left','VerticalAlignment','Middle','Fontsize',fontsize(2),'Color','k');

w1 = 0.1;
w2 = 0.54;
s = 0.05;
xpos(1) = 0.13;
xpos(2) = xpos(1)+w1+s;

h1 = 0.1;
h2 = 0.6;
ypos(1) = 0.1;
ypos(2) = xpos(1)+h1+s;

hdl_a(1) = axes('Position',[xpos(2),ypos(2),w2,h2],'Parent',Hp1,'FontSize',fontsize(3));
hdl_i = imagesc(nan);
set(hdl_a(1),'Box','On','DataAspectRatio',[1 1 1]);
set(hdl_a(1),'Xtick',{},'Ytick',{});
colorbar
set(hdl_a(1),'Position',[xpos(2),ypos(2),w2,h2]);
colormap(hdl_a(1),cmap);
% plot the ROI
hdl_roi = patch('Vertices',nan(4,2),'Faces',1:4,'FaceColor','none','EdgeColor',color.col1,'Parent',hdl_a(1));

hdl_a(2) = axes('Position',[xpos(1),ypos(2),w1,h2],'NextPlot','Add','Box','On','Parent',Hp1,'FontSize',fontsize(3),'ydir','reverse');
hdl_p(1) = plot(nan,nan,'-','Color',color.col1);
hdl_p(2) = plot(nan,nan,'-','Color',color.col2);
legend(hdl_p(1:2),{'$u_{hx}$','$u_{hy}$'},'Orientation','horizontal','Interpreter','Latex','Location','NorthOutside','FontSize',fontsize(3));
set(hdl_a(2),'Position',[xpos(1),ypos(2),w1,h2]);
xlabel('$u_h$ [px]','Interpreter','Latex','FontSize',fontsize(3))
ylabel('$y$ [px]','Interpreter','Latex','FontSize',fontsize(3))

hdl_a(3) = axes('Position',[xpos(2),ypos(1),w2,h1],'NextPlot','Add','Box','On','Parent',Hp1,'FontSize',fontsize(3));
hdl_p(3) = plot(nan,nan,'-','Color',color.col2);
hdl_p(4) = plot(nan,nan,'-','Color',color.col1);
legend(hdl_p(3:4),{'$u_{vx}$','$u_{vy}$'},'Orientation','vertical','Interpreter','Latex','Location','EastOutside','FontSize',fontsize(3));
ylabel('$u_v$ [px]','Interpreter','Latex','FontSize',fontsize(3))
xlabel('$x$ [px]','Interpreter','Latex','FontSize',fontsize(3))
set(hdl_a(3),'Position',[xpos(2),ypos(1),w2,h1]);

% opt.handles.figure = Hf;
% opt.handles.axes = hdl_a;
% opt.handles.images = hdl_i;
% opt.handles.plot = hdl_p;
% opt.handles.roi = hdl_roi;
% opt.handles.status = Hs;
% opt.handles.abort = Hb.abort;


% Button Callbacks
% ===========================================================

% ---------------------------------------------
    function set_save_path(varargin)
        
        if isempty(savepath)
            pth = currentpath;
        else
            pth = savepath;
        end
        
        appendstat('No save folder defined, please choose one in the <Preparation> tab.');
        pth = uigetdir(pth,'Select a folder to save the results');
        if pth ~= 0
            savepath = pth;
            Hb.savepath.String = savepath;
            appendstat('New save folder: %s',savepath);
        end
    end

% ---------------------------------------------
    function go(varargin)
        appendstat('Go...');
        Nf = numel(filelist);
        
        if Nf < 2
            appendstat('Need at least 2 files, aborted');
            Hb.go.Value = 0;
            return
        end
        
        if mod(Nf,2) == 1
            appendstat('Odd number of files detected, aborted');
            Hb.go.Value = 0;
            return
        end
        
        if isempty(savepath)
            set_save_path();
        end
        
        % parse the options
        savepath = Hb.savepath.String;
        if ~isempty(savepath) && ~exist(savepath,'dir')
            ButtonName = questdlg('The directory does not exist, should it be created?', ...
                'Create directory?', ...
                'Yes', 'No', 'No');
            if strcmpi(ButtonName,'Yes')
                mkdir(savepath);
            end
            Hb.savepath.String = savepath;
        end
        
        if isempty(savepath)
            appendstat('No save folder defined, job aborted');
            return
        end
        
        if ~exist(savepath,'dir')
            appendstat('Save folder does not exist, job aborted');
            return
        end
        
        % process options
        % --------------
        
        % read the step data
        opt = ini;
        
        % fixed options
        opt.Verbose = 1;
        opt.Plot = 1;
        opt.gui = true;

        % some handles to transfer
        opt.handles.figure = Hf;
        opt.handles.axes = hdl_a;
        opt.handles.images = hdl_i;
        opt.handles.plot = hdl_p;
        opt.handles.roi = hdl_roi;
        opt.handles.status = Hs;
        opt.handles.abort = Hb.abort;        
        
        % get the options from the uicontrols
        opt.TIFdepth = tiffdepths{Hb.TIFdepth.Value};
        opt.SaveCOR = Hb.SaveCOR.Value;
        opt.SaveRES = Hb.SaveRES.Value;
        opt.SavePNG = Hb.SavePNG.Value;
        opt.SaveRBT = Hb.SaveRBT.Value;
        opt.SaveLOG = Hb.SaveLOG.Value;
        
        opt.RotateH = round(eval(Hb.RotateH.String));
        opt.RotateV = round(eval(Hb.RotateV.String));
        opt.DatabarRemove = Hb.DatabarRemove.Value;
        opt.ROI = eval(Hb.ROI.String);
        
        opt.RBT_flag = Hb.RBT_flag.Value;
        opt.RBT_Window = eval(Hb.RBT_Window.String); % tukeywin
        opt.Steps_flag = Hb.Steps_flag.Value;
        
        % number of correction pairs
        Nc = Nf/2;
        Hi.progress.CData = ones(1,Nc,3);
        Ha.progress.XLim = [0.5,Nc+0.5];
        
        appendstat(sprintf('Starting the correction of %d pairs of images',Nc));
        appendstat('','','');
        for k = 1:Nc
            
            Hi.progress.CData(1,k,:) = color.progress2;
            Ht.progress.String = sprintf('Processing image pair %d/%d',k,Nc);
            
            kh = 1 + 2*(k-1);
            kv = 2 + 2*(k-1);
            
            % generate a basename
            [~, basenameH, ~] = fileparts(filelist{kh});
            [~, basenameV, ~] = fileparts(filelist{kv});

            % use the first image as basename for the saved files
            opt.Savename = fullfile(savepath,basenameH);
            
            appendstat(rule);
            appendstat(sprintf('Correction %d/%d',k,Nc));
            appendstat(rule);
            appendstat(sprintf('H: %s',basenameH));
            appendstat(sprintf('V: %s',basenameV));
            appendstat('');
            
            % read the two images
            gh = imread(filelist{kh});
            gv = imread(filelist{kv});
            
            siz = size(gh);
            Hi.liveview.CData = gh;
            Ha.liveview.XLim = [1 siz(2)];
            Ha.liveview.YLim = [1 siz(1)];
            
            % do the hard work
            out = scan_corr(gh,gv,opt);
            
            if (Hb.abort.Value == 1)
                appendstat(sprintf('Aborted at image pair %d/%d',k,Nc));
                Ht.progress.String = sprintf('Aborted at image pair %d/%d',k,Nc);
                Hb.abort.Value = 0;
                return
            end
            
            if Hb.SaveMAT.Value
                save(fullfile(savepath,[basenameH, '.mat']),'-v7.3','-struct','out');
            end
            
            Hi.progress.CData(1,k,:) = color.progress1;
            
        end
        
        appendstat(rule);
        appendstat('Job done...');
        appendstat(rule);
        Ht.progress.String = 'Job done...';
        Hb.go.Value = 0;
    end

% ---------------------------------------------
    function abort(varargin)
    end

% ---------------------------------------------
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

% ---------------------------------------------
    function image_fun(varargin)
        % Buttons in the image panel
        
        type = varargin{3};
        if strcmpi(type,'select')
            % called upon events in the image listbox
            if isempty(filelist)
                return
            end
            selected = Hb.imlist.Value;
            if numel(selected) > 1
                return
            end
            if selected > numel(filelist)
                return
            end
            im = imread(filelist{selected});
            im = databar_remove(im,0.99);
            siz = size(im);
            Hi.preview.CData = im;
            Ha.preview.XLim = [1 siz(2)];
            Ha.preview.YLim = [1 siz(1)];
            
            xc = 0.5*(1+siz(2));
            yc = 0.5*(1+siz(1));
            L = 0.5*min(siz(1:2));
            vx = L*[-1; 1; 1; -1] + xc;
            vy = L*[-1; -1; 1; 1] + yc;
            set(Hp.cent,'Vertices',[vx,vy]);
            
        elseif strcmpi(type,'add')
            % called upon pressing the Add button
            
            % popup, select files
            [filename, filepath] = uigetfile({imtypstr,'All Image Files';'*.*','All Files' },'Select an Image',currentpath,'MultiSelect','on');
            if ~iscell(filename)
                filename = {filename};
            end
            % if cancel
            if filename{1} == 0
                appendstat('Add: no file selected');
                return
            end
            % remember the current path
            currentpath = filepath;
            for k = 1:numel(filename)
                Nf = numel(filelist);
                filelist{Nf+1,1} = fullfile(filepath,filename{k});
            end
            update_filelist();
            Hb.imlist.Value = numel(filelist);
            image_fun([],[],'select')
            if k == 1
                appendstat('Add: 1 file added');
            else
                appendstat(sprintf('Add: %d files added',k));
            end
        elseif strcmpi(type,'del')
            % called upon pressing the Del button
            
            % remove files from the list (and Delete measurements)
            selected = Hb.imlist.Value;
            if ~isempty(selected)
                % clean any existing plot objects
                Nf = numel(filelist);
                notselected = setdiff(1:Nf,selected);
                filelist = filelist(notselected);
                Hb.imlist.Value = 1;
                update_filelist();
                if numel(selected) == 1
                    appendstat('Remove: 1 file removed from the list');
                else
                    appendstat(sprintf('Remove: %d files removed from the list',numel(selected)));
                end
            end
        elseif strcmpi(type,'up')
            % called upon pressing the Up button
            selected = Hb.imlist.Value;
            if isempty(selected)
                return
            end
            Ns = numel(selected);
            Nf = numel(filelist);
            
            order = 1:Nf;
            newselect = selected;
            for k = 1:Ns
                kk = selected(k);
                if kk == 1
                    continue
                end
                if (kk-1 >= 1) && ~any(selected == order(kk-1))
                    newselect(k) = selected(k)-1;
                    order([kk-1, kk]) = order([kk, kk-1]);
                end
            end
            filelist = filelist(order);
            Hb.imlist.Value = newselect;
            update_filelist();
        elseif strcmpi(type,'down')
            % called upon pressing the Up button
            selected = Hb.imlist.Value;
            if isempty(selected)
                return
            end
            Ns = numel(selected);
            Nf = numel(filelist);
            
            order = 1:Nf;
            newselect = selected;
            for k = Ns:-1:1
                kk = selected(k);
                if kk == Nf
                    continue
                end
                if (kk+1 <= Nf) && ~any(selected == order(kk+1))
                    newselect(k) = selected(k)+1;
                    order([kk+1, kk]) = order([kk, kk+1]);
                end
            end
            filelist = filelist(order);
            Hb.imlist.Value = newselect;
            update_filelist();
        elseif strcmpi(type,'interlace')
            % called upon pressing the Interlace button
            selected = Hb.imlist.Value;
            if isempty(selected)
                return
            end
            
            Ns = numel(selected);
            if mod(Ns,2) ~= 0
                % selection must be even
                return
            end
            
            newselect = reshape(reshape(selected,Ns/2,2).',1,Ns);
            filelist(selected) = filelist(newselect);
            update_filelist();
        elseif strcmpi(type,'flip')
            % called upon pressing the Interlace button
            selected = Hb.imlist.Value;
            if isempty(selected)
                return
            end
            
            Ns = numel(selected);
            if mod(Ns,2) ~= 0
                % selection must be even
                return
            end
            
            newselect = reshape(flipud(reshape(selected,2,Ns/2)),1,Ns);
            filelist(selected) = filelist(newselect);
            update_filelist();
        end
    end

% ---------------------------------------------
    function update_filelist()
        Nf = numel(filelist);
        if numel(Hb.imlist.String) > Nf
            Hb.imlist.String = Hb.imlist.String(1:Nf);
        end
        for k = 1:Nf
            [~,filename,~] = fileparts(filelist{k});
            if mod(k,2) == 1
                hvstr = 'H';
            else
                hvstr = 'V';
            end
            if any(mod(k,4) == [1,2])
                fontcol = color.fontA;
            else
                fontcol = color.fontB;
            end
            Hb.imlist.String{k} = sprintf('<html><font color="%s">%03d %s: %s</font></html>',fontcol,ceil(k/2),hvstr,filename);
        end
    end


% ---------------------------------------------
    function steplist_update(varargin)
        SO = zeros(numel(step_opts),Nsteps);
        for k = 1:numel(step_opts)
            SO(k,:) = ini.(step_opts{k});
        end
        
        for k = 1:Nsteps
            frmt = 'Step %i [%2d,%4d,%3d,%5.1f, %7.1e,%3d,%3d,%3d, %d]';
            StepList{k,1} = sprintf(frmt,k,SO(:,k));
        end
    end


% ---------------------------------------------
    function step_fun(varargin)        
        type = varargin{3};
        if strcmpi(type,'select')
            selected = Hb.StepList.Value;
            if isempty(selected)
                return
            end
            i = selected(1);
            
            % update the shown settings
            Hb.LineThickness.String = num2str(ini.LineThickness(i));
            Hb.PolyOrder.Value = ini.PolyOrder(i)+1;
            Hb.Bin.String = num2str(ini.Bin(i));
            Hb.Blur.String = num2str(ini.Blur(i));
            Hb.IterMax.String = num2str(ini.IterMax(i));
            Hb.SubIterMax.String = num2str(ini.SubIterMax(i));
            Hb.TrustRegion.String = num2str(ini.TrustRegion(i));
            Hb.ConvCrit.String = num2str(ini.ConvCrit(i));
            Hb.InterpMethod.Value = ini.InterpMethod(i)-1;
            
            
        elseif strcmpi(type,'add')
            N = numel(StepList);
            Nsteps = N + 1;
            
            ini.LineThickness(N+1) = 1;
            ini.PolyOrder(N+1) = 0;
            ini.Bin(N+1) = 1;
            ini.Blur(N+1) = 0;
            ini.IterMax(N+1) = 30;
            ini.SubIterMax(N+1) = 3;
            ini.TrustRegion(N+1) = 2;
            ini.ConvCrit(N+1) = 1e-4;
            ini.InterpMethod(N+1) = 5;

            steplist_update();
            Hb.StepList.String = StepList;
            Hb.StepList.Value = Nsteps;
            step_fun([],[],'select');
            
        elseif strcmpi(type,'copy')
            sel = Hb.StepList.Value;
            if isempty(sel)
                return
            end
            
            N = numel(StepList);
            
            newsel = [];
            val = [];
            for k = 1:N
                newsel = [newsel, k];
                if any(k == sel)
                    newsel = [newsel, k];
                    val = [val, numel(newsel)];
                end
            end
            
            ini.LineThickness = ini.LineThickness(newsel);
            ini.PolyOrder = ini.PolyOrder(newsel);
            ini.Bin = ini.Bin(newsel);
            ini.Blur = ini.Blur(newsel);
            ini.IterMax = ini.IterMax(newsel);
            ini.SubIterMax = ini.SubIterMax(newsel);
            ini.TrustRegion = ini.TrustRegion(newsel);
            ini.ConvCrit = ini.ConvCrit(newsel);
            ini.InterpMethod = ini.InterpMethod(newsel);

            Nsteps = numel(newsel);
            
            steplist_update();
            Hb.StepList.String = StepList;
            Hb.StepList.Value = val;
            step_fun([],[],'select');
            
            
            
        elseif strcmpi(type,'del')
            % remove files from the list (and Delete measurements)
            selected = Hb.StepList.Value;
            if ~isempty(selected)
                % clean any existing plot objects
                N = numel(StepList);
                notselected = setdiff(1:N,selected);
                
                ini.LineThickness = ini.LineThickness(notselected);
                ini.PolyOrder = ini.PolyOrder(notselected);
                ini.Bin = ini.Bin(notselected);
                ini.Blur = ini.Blur(notselected);
                ini.IterMax = ini.IterMax(notselected);
                ini.SubIterMax = ini.SubIterMax(notselected);
                ini.TrustRegion = ini.TrustRegion(notselected);
                ini.ConvCrit = ini.ConvCrit(notselected);
                ini.InterpMethod = ini.InterpMethod(notselected);
                
                Nsteps = numel(notselected);
                StepList = StepList(notselected);
                steplist_update();
                Hb.StepList.String = StepList;
                Hb.StepList.Value = [];
                step_fun([],[],'select');
            end
        elseif strcmpi(type,'up')
            selected = Hb.StepList.Value;
            if isempty(selected)
                return
            end
            Ns = numel(selected);
            Nf = numel(StepList);
            
            order = 1:Nf;
            newselect = selected;
            for k = 1:Ns
                kk = selected(k);
                if kk == 1
                    continue
                end
                if (kk-1 >= 1) && ~any(selected == order(kk-1))
                    newselect(k) = selected(k)-1;
                    order([kk-1, kk]) = order([kk, kk-1]);
                end
            end
            
            ini.LineThickness = ini.LineThickness(order);
            ini.PolyOrder = ini.PolyOrder(order);
            ini.Bin = ini.Bin(order);
            ini.Blur = ini.Blur(order);
            ini.IterMax = ini.IterMax(order);
            ini.SubIterMax = ini.SubIterMax(order);
            ini.TrustRegion = ini.TrustRegion(order);
            ini.ConvCrit = ini.ConvCrit(order);
            ini.InterpMethod = ini.InterpMethod(order);
            
            StepList = StepList(order);
            
            steplist_update();
            Hb.StepList.String = StepList;
            Hb.StepList.Value = newselect;
            step_fun([],[],'select');
            
        elseif strcmpi(type,'down')
            selected = Hb.StepList.Value;
            if isempty(selected)
                return
            end
            Ns = numel(selected);
            Nf = numel(StepList);
            
            order = 1:Nf;
            newselect = selected;
            for k = Ns:-1:1
                kk = selected(k);
                if kk == Nf
                    continue
                end
                if (kk+1 <= Nf) && ~any(selected == order(kk+1))
                    newselect(k) = selected(k)+1;
                    order([kk+1, kk]) = order([kk, kk+1]);
                end
            end
            
            ini.LineThickness = ini.LineThickness(order);
            ini.PolyOrder = ini.PolyOrder(order);
            ini.Bin = ini.Bin(order);
            ini.Blur = ini.Blur(order);
            ini.IterMax = ini.IterMax(order);
            ini.SubIterMax = ini.SubIterMax(order);
            ini.TrustRegion = ini.TrustRegion(order);
            ini.ConvCrit = ini.ConvCrit(order);
            ini.InterpMethod = ini.InterpMethod(order);
            
            StepList = StepList(order);
            
            steplist_update();
            Hb.StepList.String = StepList;
            Hb.StepList.Value = newselect;
            step_fun([],[],'select');
        elseif strcmpi(type,'opt')
            sel = Hb.StepList.Value;
            if isempty(sel)
                return
            end
            
            ini.LineThickness(sel) = eval(Hb.LineThickness.String);                       
            ini.PolyOrder(sel) = Hb.PolyOrder.Value - 1;            
            ini.Bin(sel) = eval(Hb.Bin.String);
            ini.Blur(sel) = eval(Hb.Blur.String);
            ini.IterMax(sel) = eval(Hb.IterMax.String);
            ini.SubIterMax(sel) = eval(Hb.SubIterMax.String);
            ini.TrustRegion(sel) = eval(Hb.TrustRegion.String);
            ini.ConvCrit(sel) = eval(Hb.ConvCrit.String);
            ini.InterpMethod(sel) = Hb.InterpMethod.Value + 1;

            steplist_update();
            Hb.StepList.String = StepList;
            

        else
            error('unknown %s',type)

        end
        
    end

% ---------------------------------------------
    function help_fun(varargin)
        % this function creates a new window with the following text inside
        % if the help is not open, otherwise it closes the help window
        helpfile = 'scan_corr_GUI_Help.txt';
        if exist(helpfile,'file')
            fid = fopen(helpfile,'rt');
            txt = textscan(fid,'%s','Whitespace','','Delimiter','\n');
            fclose(fid);
            txt = txt{1};
        else
            txt = {'Help file is missing, please locate scan_corr_GUI_Help.txt or contact jan.neggers@centralesupelec.fr'};
        end
        set(Hhelp,'String',txt)
    end

% ----------------------------------
    function ini = iniread(varargin)
        if ~exist(inifile,'file')
            error('inifile not found')
        end
        
        % read the entire option file into memory
        fid = fopen(inifile,'rt');
        if fid == -1
            error('cannot read from file: %s',inifile)
        end
        F = textscan(fid,'%s','Delimiter','\n','Whitespace','');
        fclose(fid);
        F = F{1};
        
        % Remove empty lines
        I = regexpi(F,'.*');
        I = ~cellfun('isempty',I);
        F = F(I);
        
        % Remove comment lines
        I = regexpi(F,'^[%#\[].*');
        I = cellfun('isempty',I);
        F = F(I);
        
        % Remove inline comments
        I = regexpi(F,'\s*[%#].*');
        J = ~cellfun('isempty',I);
        n = 1:length(J);
        for k = n(J)
            % keep only the characters before the match
            F{k} = F{k}(1:I{k}-1);
        end
        
        % Remove extra white spaces
        F = strtrim(F);
        
        % process to a structure
        for k = 1:length(F)
            % one line
            A = F{k};
            try
                % evaluate the line
                S = regexpi(A,'^\s*(.*?)\s*=\s*(.*?)\s*$','tokens');
                name = S{1}{1};
                value = S{1}{2};
            catch
                fprintf(2,'!!! iniread: error processing the line [%s]\n',A);
                return
            end
            
            ini.(name) = eval(value);
        end
    end

% ----------------------------------
    function iniwrite(ini)
        str = {;
            '% These are the default settings for the scan_corr_GUI'
            '% This file will be automatically generated if it does not exist, i.e. delete it to reset the settings.'
            ''};
        
        str{end+1} = '[Saving Results]';
        str{end+1} = sprintf('%s = %d %% %s','SaveLOG',ini.SaveLOG,'0 for false, 1 for true, save the logfile');
        str{end+1} = sprintf('%s = %d %% %s','SaveCOR',ini.SaveCOR,'0 for false, 1 for true, save the corrected image');
        str{end+1} = sprintf('%s = %d %% %s','SaveRES',ini.SaveRES,'0 for false, 1 for true, save the residual image');
        str{end+1} = sprintf('%s = %d %% %s','SavePNG',ini.SavePNG,'0 for false, 1 for true, save the final screen');
        str{end+1} = sprintf('%s = %d %% %s','SaveMAT',ini.SaveMAT,'0 for false, 1 for true, save the data as .mat');
        str{end+1} = sprintf('%s = %d %% %s','SaveRBT',ini.SaveRBT,'0 for false, 1 for true, save the RBT images');
        str{end+1} = sprintf('%s = %d %% %s','TIFdepth',ini.TIFdepth,'1 for 8bit, 2 for 16bit');
        
        str{end+1} = '[Image Preparation]';
        str{end+1} = sprintf('%s = %d %% %s','RotateH',ini.RotateH,'number of 90 degree CC rotations for H');
        str{end+1} = sprintf('%s = %d %% %s','RotateV',ini.RotateV,'number of 90 degree CC rotations for V');
        str{end+1} = sprintf('%s = %d %% %s','DatabarRemove',ini.DatabarRemove,'auto remove the databar [0=false,1=true]');
        str{end+1} = sprintf('%s = %s %% %s','ROI',mat2str(ini.ROI),'the relative border around the ROI');
        str{end+1} = sprintf('%s = %d %% %s','ROI_NAN',ini.ROI_NAN,'paint pixels outside the ROI as nan in the stored iamges');
        
        str{end+1} = '[RBT correction]';
        str{end+1} = sprintf('%s = %d %% %s','RBT_flag',ini.RBT_flag,'enable (1) or disable (0) this part');
        str{end+1} = sprintf('%s = %g %% %s','RBT_Window',ini.RBT_Window,'the radius of the tukeywin');
        
        str{end+1} = '[Steps]';
        str{end+1} = sprintf('%s = %d %% %s','Steps_flag',ini.Steps_flag,'enable (1) or disable (0) this part');
        str{end+1} = sprintf('%s = %s %% %s','LineThickness',mat2str(ini.LineThickness),'group N pixel rows together to make a thick line that moves as one');
        str{end+1} = sprintf('%s = %s %% %s','PolyOrder',mat2str(ini.PolyOrder),'The order of the polynomial basis 0 <= p <= 5 (use 0 to disable)');
        str{end+1} = sprintf('%s = %s %% %s','Bin',mat2str(ini.Bin),'bin the images to superpixels of size N');
        str{end+1} = sprintf('%s = %s %% %s','Blur',mat2str(ini.Blur),'blur the images with radius R');
        str{end+1} = sprintf('%s = %s %% %s','TrustRegion',mat2str(ini.TrustRegion),'limit the iterative update to N pixels');
        str{end+1} = sprintf('%s = %s %% %s','ConvCrit',mat2str(ini.ConvCrit),'convergence criteria');
        str{end+1} = sprintf('%s = %s %% %s','IterMax',mat2str(ini.IterMax),'max number of iterations');
        str{end+1} = sprintf('%s = %s %% %s','SubIterMax',mat2str(ini.SubIterMax),'max number of sub-iterations');
        str{end+1} = sprintf('%s = %s %% %s','InterpMethod',mat2str(ini.InterpMethod),'the bspline interpolation order 2 <= s <= 9');
        
        fid = fopen(inifile,'w+t');
        if fid == -1
            error('cannot write to file: %s',inifile)
        end
        fprintf(fid,'%s\n',str{:});
        fclose(fid);
    end

% -------------------------------
    function defaults_fun(varargin)
        type = varargin{3};
        if strcmpi(type,'save')
            % default options
            ini.TIFdepth = Hb.TIFdepth.Value;
            ini.SaveCOR = Hb.SaveCOR.Value;
            ini.SaveRES = Hb.SaveRES.Value;
            ini.SavePNG = Hb.SavePNG.Value;
            ini.SaveRBT = Hb.SaveRBT.Value;
            ini.SaveLOG = Hb.SaveLOG.Value;
            ini.SaveMAT = Hb.SaveMAT.Value;
            
            ini.RotateH = round(eval(Hb.RotateH.String));
            ini.RotateV = round(eval(Hb.RotateV.String));
            ini.DatabarRemove = Hb.DatabarRemove.Value;
            ini.ROI = eval(Hb.ROI.String);
            ini.ROI_NAN = false;
            
            ini.RBT_flag = Hb.RBT_flag.Value;
            ini.RBT_Window = eval(Hb.RBT_Window.String); % tukeywin
            ini.Steps_flag = Hb.Steps_flag.Value;
                        
            iniwrite(ini)
        elseif strcmpi(type,'reset')
            iniwrite(defini)
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


% -------------------------------
    function appendstat(varargin)
        stat = [varargin(:) ; Hs.String];
        Hs.String = stat;
    end

% -------------------------------
    function roi_select(varargin)
        v = [];
        Hf.Pointer = 'crosshair';
        Hf.WindowButtonMotionFcn = @roi_move;
        Hf.WindowButtonDownFcn = @roi_click;
    end
% =====================================================
    function roi_move(hObject, varargin)
        % function used when moving the mouse with NO button pressed
        
        p = get(Ha.preview,'CurrentPoint');
        p = p(1,1:2);
        
        if size(v,1) == 1
            V = [v ; p];
            C = [min(V(:,1)), max(V(:,1)), min(V(:,2)), max(V(:,2))];
            set(Hp.roi,'Vertices',C([1,3;2,3;2,4;1,4]));
        end
    end

% =====================================================
    function roi_click(hObject, varargin)
        % function used when moving the mouse with NO button pressed
        
        p = get(Ha.preview,'CurrentPoint');
        p = p(1,1:2);
        
        v = [v ; p];
        if size(v,1) >= 2
            
            ROI = round([min(v(:,1)), max(v(:,1)), min(v(:,2)), max(v(:,2))]);
            set(Hp.roi,'Vertices',ROI([1,3;2,3;2,4;1,4]));
            c = Hp.cent.Vertices;
            
            % [left right top bottom]
            roi(1) = (ROI(1)-c(1,1));
            roi(2) = (c(2,1)-ROI(2));
            roi(3) = (ROI(3)-c(2,2));
            roi(4) = (c(3,2)-ROI(4));
            Hb.ROI.String = sprintf('[%3.1f, %3.1f, %3.1f, %3.1f]',roi);
            
            Hf.WindowButtonMotionFcn = [];
            Hf.WindowButtonDownFcn = [];
            v = [];
            set(Hf,'Pointer','arrow');
        end
    end

% =====================================================
% End of GUI
% =====================================================

% uiwait(Hf);

% prepare the outputs
if nargout == 1
    varargout{1} = [];
end

% End of main function
end

