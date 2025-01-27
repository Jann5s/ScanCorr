clear; close all;

if true
    % debug = '-g';
    debug = '';
    
    cppfile = 'scan_corr_mex.cpp';
    mex(cppfile,debug,'-largeArrayDims','-DNOMINMAX')
    
end

if true
    % debug = '-g';
    debug = '';
    
    cppfile = 'scan_corr_im_mex.cpp';
    mex(cppfile,debug,'-largeArrayDims','-DNOMINMAX')
end

return

filegv = 'test/test_corrected_gv.tiff';
filegh = 'test/test_corrected_gh.tiff';
gv = single(imread(filegv));
gh = single(imread(filegh));


[n, m] = size(gv);

% gv = gv(1:round(0.9*n),:);
% gh = gh(1:round(0.9*n),:);
% [n, m] = size(gv);

gv = gv(:,1:round(0.9*m));
gh = gh(:,1:round(0.9*m));
[n,m] = size(gv);

% debug g
% g = f([2:n,n],[2:m,1]);

gv = gv ./ 255;
gh = gh ./ 255;

% gv = (gv - mean(gv(:))) ./ 255;
% gh = (gh - mean(gh(:))) ./ 255;

% roi [left, right, top, bottom]
roi = [10, m-20+1, 30, n-40+1];
roih = roi(1:2);
roiv = roi(3:4);

Iroi = find(1:n <= roi(3) | 1:n >= roi(4));
Jroi = find(1:m <= roi(1) | 1:m >= roi(2));

% gv(Iroi,:) = NaN;
% gh(Iroi,:) = NaN;
% gv(:,Jroi) = NaN;
% gh(:,Jroi) = NaN;

figure;
subplot(2,2,1)
imagesc(gv)
colorbar
subplot(2,2,2)
imagesc(gh)
colorbar
subplot(2,2,3)
imagesc(gv-gh)
colorbar
drawnow

figure;
hi = imagesc(gv-gh);
ha = gca;
colorbar
caxis(0.2*[-1,1])


trustregion = 2;

Nit = 20;
Nsubit = 5;
res = [];

blur =   [ 2,  2,  1, 1, 0.5, 0];
bsteps = [50, 20, 10, 5, 2  , 1];
Nbin = numel(bsteps);

s = 5;
threads = 0;

convcrit = [1e-4, 5e-5];

% initialize
% f = gv;
% ght = gh;
% gvt = gv;
Ah = zeros(n,2);
Av = zeros(m,2);
tic
[f, R] = scan_corr_im_mex(gh,gv,Ah,Av,s,threads);
toc

f0 = f;
R0 = R;
figure;
imagesc(f0);

figure;
hp = plot(1:n, Ah, 1:m, Av);

t = zeros(2,1);
c = zeros(2,1);

tic
for i = 1:Nbin
    
    gvb = image_blur(gv,blur(i));
    ghb = image_blur(gh,blur(i));
   
    % create the h-blocks
    binh = zeros(n,1);
    Nrh = roiv(2) - roiv(1);
    Nh = floor(Nrh / bsteps(i));
    ind = floor(linspace(1,Nh+1,Nrh+1));
    binh( roiv(1)+(1:Nrh) ) = ind(1:Nrh);
        
    % convert Ah to ah
    ah = zeros(2*Nh,1);
    Ixh = 1:2:2*Nh;
    Iyh = 2:2:2*Nh;
    for k = 1:Nh
        ah(Ixh(k)) = mean(Ah(binh == k,1));
        ah(Iyh(k)) = mean(Ah(binh == k,2));        
    end

    % create the v-blocks
    binv = zeros(m,1);
    Nrv = roih(2) - roih(1);
    Nv = floor(Nrv / bsteps(i));
    ind = floor(linspace(1,Nv+1,Nrv+1));
    binv( roih(1)+(1:Nrv) ) = ind(1:Nrv);
    
    % convert Av to av
    av = zeros(2*Nv,1);
    Ixv = 1:2:2*Nv;
    Iyv = 2:2:2*Nv;
    for k = 1:Nv
        av(Ixv(k)) = mean(Av(binv == k,1));
        av(Iyv(k)) = mean(Av(binv == k,2));        
    end
        
    for it = 1:Nit
        
        for mode = [0, 1]
            
            for subit = 1:Nsubit
                if mode == 0
                    N = Nh;
                    tic
                    [M, b, r] = scan_corr_mex(f,ghb,ah,binh,roih,mode,s,threads);
                    t(1) = t(1) + toc;
                    c(1) = c(1) + 1;
                else
                    N = Nv;
                    tic
                    [M, b, r] = scan_corr_mex(f,gvb,av,binv,roiv,mode,s,threads);
                    t(2) = t(2) + toc;
                    c(2) = c(2) + 1;
                end
                
                
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

                
                da(isnan(da)) = 0;
                
                if mode == 0
                    % force da(y) to be mean zero
                    da(Iyh) = da(Iyh) - mean(da(Iyh));
                    ah = ah + da;
                else
                    % force da(x) to be mean zero
                    da(Ixv) = da(Ixv) - mean(da(Ixv));
                    av = av + da;
                end
                
                da = sign(da) .* min(abs(da),trustregion);

                % convert ah to Ah
                for k = 1:Nh
                    Ah(binh == k,1) = ah(Ixh(k));
                    Ah(binh == k,2) = ah(Iyh(k));
                end
                
                % convert av to Av
                for k = 1:Nv
                    Av(binv == k,1) = av(Ixv(k));
                    Av(binv == k,2) = av(Iyv(k));
                end
                
                set(hp(1),'YData',Ah(:,1));
                set(hp(2),'YData',Ah(:,2));
                set(hp(3),'YData',Av(:,1));
                set(hp(4),'YData',Av(:,2));
                drawnow
                
                fprintf('%3d: %3d-%d-%02d, r: %10.4e, da: %10.4e\n',bsteps(i),it,mode,subit,r,rms(da));
                
                res(end+1) = r;
                
                if rms(da) < convcrit(1)
                    break
                end
                
            end
            
            % update the reference image
            [f, R] = scan_corr_im_mex(ghb,gvb,Ah,Av,s,threads);
            set(hi,'CData',R);
            set(ha,'CLim',5 * r *[-1, 1]);
            drawnow
            % f = 0.5*(ght + gvt);
        end
        if rms(da) < convcrit(2)
            break
        end
        
    end
    
    % convert ah to Ah
    for k = 1:Nh
        Ah(binh == k,1) = ah(Ixh(k));
        Ah(binh == k,2) = ah(Iyh(k));
    end
    
    % convert av to Av
    for k = 1:Nv
        Av(binv == k,1) = av(Ixv(k));
        Av(binv == k,2) = av(Iyv(k));
    end
    
end
toc

%% extrapolate the displacements
Ni = round(0.1*min(n,m));

% h
x = linspace(-1,1,n).';
M = [ones(n,1), x, x.^2];

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
M = [ones(m,1), x, x.^2];

I1 = find(binv > 0,Ni,'first');
I2 = 1:I1(1)-1;
I3 = find(binv > 0,Ni,'last');
I4 = I3(end)+1 : m;

a = M(I1,:) \ Av(I1,:);
Av(I2,:) = M(I2,:) * a;

a = M(I3,:) \ Av(I3,:);
Av(I4,:) = M(I4,:) * a;

% show the new data
set(hp(1),'YData',Ah(:,1));
set(hp(2),'YData',Ah(:,2));
set(hp(3),'YData',Av(:,1));
set(hp(4),'YData',Av(:,2));
drawnow

tic
[f1, R1] = scan_corr_im_mex(gh,gv,Ah,Av,s,threads);
toc

figure;
imagesc(f1);

disp(t./c)
