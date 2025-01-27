clear; close all

if true
    mex('scan_corr_mex.cpp','-largeArrayDims','-DNOMINMAX');
    mex('scan_corr_im_mex.cpp','-largeArrayDims','-DNOMINMAX');
end

% filegh = '../data/Etape6-testRecharge_ech1_stepA15_0deg.tif';
% filegv = '../data/Etape6-testRecharge_ech1_step15_90deg.tif';

% filegh = '../data/MATMECA_Paris_Saclay_A01_01H.tif';
% filegv = '../data/MATMECA_Paris_Saclay_A01_02V.tif';

filegh = '../data/MATMECA_Paris_Saclay_A02_01H.tif';
filegv = '../data/MATMECA_Paris_Saclay_A02_02V.tif';

% filegh = '../data/MATMECA_Paris_Saclay_E03_02H.tif';
% filegv = '../data/MATMECA_Paris_Saclay_E03_01V.tif';

% filegh = '../virt_exp/tif/orthogonality_a010_mode3_H.tif';
% filegv = '../virt_exp/tif/orthogonality_a010_mode3_V.tif';
% filegh = '../virt_exp/tif/orthogonality_a000_mode0_H.tif';
% filegv = '../virt_exp/tif/orthogonality_a000_mode0_V.tif';


gh = double(imread(filegh))./255;
gv = double(imread(filegv))./255;

% gh = gh(400:end-400,:);
% gv = gv(:,400:end-400);

ha(1) = subplot(1,2,1);
imagesc(gh)
ha(2) = subplot(1,2,2);
imagesc(gv)

set(ha,'DataAspectRatio',[1 1 1]);


[n, m] = size(gh);
roi = 50*[1, 1, 1, 1];

% bin = 1;
% f = image_coarsegrain(f,bin);
% g = image_coarsegrain(g,bin);

opt.Verbose = 2;
opt.ROI = roi;

% opt.closefig = true;
opt.Logfile = 'test/test.log';
opt.Savename = 'test/test.tif';

opt.SaveRBT = true;
opt.SaveRES = true;
opt.SaveCOR = true;
opt.SavePNG = true;

opt.Steps_flag = true;
opt.LineThickness = [1, 1, 1, 5, 2, 1];
opt.PolyOrder = [1, 2, 3, 0];
opt.Blur = [16, 8, 4, 1, 1, 0];
opt.Bin = 2;
opt.TrustRegion = 0.01;
opt.IterMax = 30;
opt.SubIterMax = 3;
opt.ConvCrit = 1e-3;
opt.InterpMethod = 3;

tic
out = scan_corr(gh,gv,opt);
toc


