function batch_analyze_4ch_tiffs2()
% REMOVE % in line46 for batch analysis
% Batch analyzer for 4-channel TIFFs:
% - Select multiple TIFFs
% - For each file: build Fig1-3 (main), Fig4 (pair overlaps), Fig5 (triple overlaps)
% - Save all three figures into a per-file subfolder
% - Collect overlap coefficients into a table & save CSV at the end

    set(0,'DefaultFigureWindowStyle','docked');
    %set(0,'DefaultFigureWindowStyle','normal');  


    % ---------- File selection ----------
    [files, folder] = uigetfile({'*.tif;*.tiff','TIFF files (*.tif, *.tiff)'}, ...
                                'Select one or more 4-channel TIFFs','MultiSelect','on');
    if isequal(files,0)
        disp('No files selected.'); return;
    end
    if ischar(files), files = {files}; end

    % ---------- Results accumulator ----------
    results = [];   % struct array
    coefNamesPairs   = {'coef12','coef13','coef23','coef14','coef24','coef34'};
    coefNamesTriples = {'coef123','coef124','coef134','coef234'};

    % ---------- Loop over files ----------
    for k = 1:numel(files)
        fname = files{k};
        srcPath = fullfile(folder, fname);
        fprintf('Processing %d/%d: %s\n', k, numel(files), fname);

        try
            [coefStruct, fMain, f4, f5] = analyze_one_tiff(srcPath); %#ok<NASGU>
        catch ME
            warning('Failed on %s: %s', fname, ME.message);
            continue;
        end

% ----- Save figures (only .fig, no PNG, no subfolders) -----
[~, base, ~] = fileparts(fname);
savefigsafe(fMain, fullfile(folder, [base '_Fig1to3.fig']));
savefigsafe(f4,   fullfile(folder, [base '_Fig4.fig']));
savefigsafe(f5,   fullfile(folder, [base '_Fig5.fig']));

% Close to keep memory in check
% close([fMain, f4, f5]);


        % ----- Stash results -----
        R = struct('file', fname);
        % copy over all pair & triple coefficients if present
        for nm = coefNamesPairs,   key = nm{1}; if isfield(coefStruct,key), R.(key) = coefStruct.(key); end, end
        for nm = coefNamesTriples, key = nm{1}; if isfield(coefStruct,key), R.(key) = coefStruct.(key); end, end
        results = [results; R]; %#ok<AGROW>
    end

    % ---------- Convert to table & save CSV ----------
    if isempty(results)
        warning('No results to save.');
        return;
    end

    % Ensure all fields exist in the table
    results = fillMissingFields(results, [coefNamesPairs, coefNamesTriples]);

    T = struct2table(results);
    assignin('base','batch_overlap_results', T);  % put in workspace

    csvPath = fullfile(folder, sprintf('batch_overlap_results_%s.csv', datestr(now,'yyyymmdd_HHMMSS')));
    writetable(T, csvPath);
    fprintf('Saved CSV: %s\n', csvPath);
end


% ======================= Single-file analyzer ============================
function [coefStruct, fMain, f4, f5] = analyze_one_tiff(srcPath)
% Builds Fig1-3 (one big main figure), Fig4 (pair overlaps), Fig5 (triple overlaps)
% Returns coefficients in coefStruct and figure handles

    [~, fname, ext] = fileparts(srcPath); %#ok<NASGU>
    info = imfinfo(srcPath);
    assert(numel(info) >= 4, 'File has only %d page(s); need >=4.', numel(info));

    ch1 = imread(srcPath,1);
    ch2 = imread(srcPath,2);
    ch3 = imread(srcPath,3);
    ch4 = imread(srcPath,4);

    % ---------- Helpers ----------
    makeRGB   = @(r,g,b) cat(3,r,g,b);
    blendGray = 0.7;                  % for Gray blending into RGB

    % ---------- Fig1 (normalized [0,1]) ----------
    ch1n = mat2gray(ch1); ch2n = mat2gray(ch2);
    ch3n = mat2gray(ch3); ch4n = mat2gray(ch4);

    RGB_all_1 = min(1, makeRGB(ch1n,ch2n,ch3n) + blendGray*repmat(ch4n,[1 1 3]));
    RGB_14_1  = makeRGB(min(1,ch1n+blendGray*ch4n), blendGray*ch4n, blendGray*ch4n);
    RGB_12_1  = makeRGB(ch1n, ch2n, zeros(size(ch1n)));
    RGB_13_1  = makeRGB(ch1n, zeros(size(ch1n)), ch3n);
    RGB_23_1  = makeRGB(zeros(size(ch1n)), ch2n, ch3n);

    % ---------- Fig2 (LUT-applied; preserve 16-bit, then display as [0,1]) ----------
    contrast16 = @(img) uint16( ...
        max(0, min(1, (double(img) - prctile(double(img(:)),90)) ./ ...
        max(eps, prctile(double(img(:)),99.9) - prctile(double(img(:)),90)))) * 65535 );

    ch1c = contrast16(ch1); ch2c = contrast16(ch2);
    ch3c = contrast16(ch3); ch4c = contrast16(ch4);

    ch1c_f = double(ch1c)/65535; ch2c_f = double(ch2c)/65535;
    ch3c_f = double(ch3c)/65535; ch4c_f = double(ch4c)/65535;

    RGB_all_2 = min(1, makeRGB(ch1c_f,ch2c_f,ch3c_f) + blendGray*repmat(ch4c_f,[1 1 3]));
    RGB_14_2  = makeRGB(min(1,ch1c_f+blendGray*ch4c_f), blendGray*ch4c_f, blendGray*ch4c_f);
    RGB_12_2  = makeRGB(ch1c_f, ch2c_f, zeros(size(ch1c_f)));
    RGB_13_2  = makeRGB(ch1c_f, zeros(size(ch1c_f)), ch3c_f);
    RGB_23_2  = makeRGB(zeros(size(ch1c_f)), ch2c_f, ch3c_f);
    RGB_24_2  = makeRGB(blendGray*ch4c_f, min(1, ch2c_f + blendGray*ch4c_f), blendGray*ch4c_f);
    RGB_34_2  = makeRGB(blendGray*ch4c_f, blendGray*ch4c_f, min(1, ch3c_f + blendGray*ch4c_f));

    % ---------- Fig3 (Otsu on LUT-applied; pseudocolor merges) ----------
    I1 = im2double(ch1c); I2 = im2double(ch2c); I3 = im2double(ch3c); I4 = im2double(ch4c);
    T1 = graythresh(I1); BW1 = imbinarize(I1, T1);
    T2 = graythresh(I2); BW2 = imbinarize(I2, T2);
    T3 = graythresh(I3); BW3 = imbinarize(I3, T3);
    T4 = graythresh(I4); BW4 = imbinarize(I4, T4);

    BW1d = double(BW1); BW2d = double(BW2); BW3d = double(BW3); BW4d = double(BW4);
    Z = zeros(size(BW1d));
    addWhite = @(RGB, W) min(1, RGB + 0.5 * repmat(W, [1 1 3]));
    RGB_14    = addWhite(cat(3, BW1d, Z,    Z   ), BW4d);
    RGB_12    =            cat(3, BW1d, BW2d, Z   );
    RGB_13    =            cat(3, BW1d, Z,    BW3d);
    RGB_23    =            cat(3, Z,    BW2d, BW3d);
    RGB_123   =            cat(3, BW1d, BW2d, BW3d);
    RGB_1234  = addWhite(RGB_123, BW4d);
    RGB_24    = addWhite(cat(3, Z,    BW2d, Z   ), BW4d);
    RGB_34    = addWhite(cat(3, Z,    Z,    BW3d), BW4d);

    % ---------- Overlaps (pairs + triples) from LUT-applied channels ----------
    epsSafe = 1e-12;
    A = double(ch1c); B = double(ch2c); C = double(ch3c); D = double(ch4c);
    norm1 = A / max(sum(A(:)), epsSafe);
    norm2 = B / max(sum(B(:)), epsSafe);
    norm3 = C / max(sum(C(:)), epsSafe);
    norm4 = D / max(sum(D(:)), epsSafe);

    % Pairs
    [overlap12,coef12] = ovMap(norm1,norm2);
    [overlap13,coef13] = ovMap(norm1,norm3);
    [overlap23,coef23] = ovMap(norm2,norm3);
    [overlap14,coef14] = ovMap(norm1,norm4);
    [overlap24,coef24] = ovMap(norm2,norm4);
    [overlap34,coef34] = ovMap(norm3,norm4);

    % Triples (element-wise min across 3)
    [overlap123,coef123] = ovMap3(norm1,norm2,norm3);
    [overlap124,coef124] = ovMap3(norm1,norm2,norm4);
    [overlap134,coef134] = ovMap3(norm1,norm3,norm4);
    [overlap234,coef234] = ovMap3(norm2,norm3,norm4);

    % ---------- Build figures ----------
    % Main figure: Fig1-3 (stacked)
    fMain = figure('Color','w','Name',sprintf('All-in-one: %s', [fname ext]), ...
                   'Units','normalized','Position',[0.03 0.05 0.94 0.90]);
    mainTL = tiledlayout(fMain, 3, 1, 'Padding','compact', 'TileSpacing','compact');

    % Section 1: Fig1 (3x3)
    tl1 = tiledlayout(mainTL, 2, 6, 'Padding','compact','TileSpacing','compact'); tl1.Layout.Tile = 1;
    nexttile(tl1); imshow(ch1n,[]);     title('ch1 (Red) – norm');
    nexttile(tl1); imshow(ch2n,[]);     title('ch2 (Green) – norm');
    nexttile(tl1); imshow(ch3n,[]);     title('ch3 (Blue) – norm');
    nexttile(tl1); imshow(ch4n,[]);     title('ch4 (Gray) – norm');
    nexttile(tl1); imshow(RGB_all_1);   title('Merged: 1+2+3+4');
    nexttile(tl1); imshow(RGB_14_1);    title('Merged: 1+4');
    nexttile(tl1); imshow(RGB_12_1);    title('Merged: 1+2');
    nexttile(tl1); imshow(RGB_13_1);    title('Merged: 1+3');
    nexttile(tl1); imshow(RGB_23_1);    title('Merged: 2+3');
    sgtitle(tl1, 'Fig1 — Original normalized', 'FontWeight','bold');

    % Section 2: Fig2 (4x3 to include 2+4 and 3+4)
    tl2 = tiledlayout(mainTL, 2, 6, 'Padding','compact','TileSpacing','compact'); tl2.Layout.Tile = 2;
    nexttile(tl2); imshow(ch1c_f,[]);   title('ch1 – LUT');
    nexttile(tl2); imshow(ch2c_f,[]);   title('ch2 – LUT');
    nexttile(tl2); imshow(ch3c_f,[]);   title('ch3 – LUT');
    nexttile(tl2); imshow(ch4c_f,[]);   title('ch4 – LUT');
    nexttile(tl2); imshow(RGB_all_2);   title('Merged: 1+2+3+4 (LUT)');
    nexttile(tl2); imshow(RGB_14_2);    title('Merged: 1+4 (LUT)');
    nexttile(tl2); imshow(RGB_12_2);    title('Merged: 1+2 (LUT)');
    nexttile(tl2); imshow(RGB_13_2);    title('Merged: 1+3 (LUT)');
    nexttile(tl2); imshow(RGB_23_2);    title('Merged: 2+3 (LUT)');
    nexttile(tl2); imshow(RGB_24_2);    title('Merged: 2+4 (LUT)');
    nexttile(tl2); imshow(RGB_34_2);    title('Merged: 3+4 (LUT)');
    sgtitle(tl2, 'Fig2 — Contrast-enhanced (90%→black, top 0.1%→white)', 'FontWeight','bold');

    % Section 3: Fig3 (2x6) Otsu on LUT-applied
    tl3 = tiledlayout(mainTL, 2, 6, 'Padding','compact','TileSpacing','compact'); tl3.Layout.Tile = 3;
    % Otsu thresholds (showing only)
    T1_16 = round(T1*65535); T2_16 = round(T2*65535); T3_16 = round(T3*65535); T4_16 = round(T4*65535);
    nexttile(tl3); imshow(BW1); title(sprintf('ch1 Otsu T=%.4f (~%d)', T1, T1_16));
    nexttile(tl3); imshow(BW2); title(sprintf('ch2 Otsu T=%.4f (~%d)', T2, T2_16));
    nexttile(tl3); imshow(BW3); title(sprintf('ch3 Otsu T=%.4f (~%d)', T3, T3_16));
    nexttile(tl3); imshow(BW4); title(sprintf('ch4 Otsu T=%.4f (~%d)', T4, T4_16));
    nexttile(tl3); imshow(RGB_123);     title('Otsu merge: 1+2+3');
    nexttile(tl3); imshow(RGB_1234);    title('Otsu merge: 1+2+3+4 (4=+50% white)');
    nexttile(tl3); imshow(RGB_14);      title('Otsu merge: 1+4');
    nexttile(tl3); imshow(RGB_12);      title('Otsu merge: 1+2');
    nexttile(tl3); imshow(RGB_13);      title('Otsu merge: 1+3');
    nexttile(tl3); imshow(RGB_23);      title('Otsu merge: 2+3');
    nexttile(tl3); imshow(RGB_24);      title('Otsu merge: 2+4');
    nexttile(tl3); imshow(RGB_34);      title('Otsu merge: 3+4');
    sgtitle(mainTL, sprintf('All views — %s', [fname ext]), 'Interpreter','none','FontSize',12,'FontWeight','bold');

    % ---------- Fig5 (triples) ----------
    RGB_123_2 = makeRGB(ch1c_f, ch2c_f, ch3c_f);
    RGB_124_2 = min(1, makeRGB(ch1c_f, ch2c_f, zeros(size(ch1c_f))) + blendGray*repmat(ch4c_f,[1 1 3]));
    RGB_134_2 = min(1, makeRGB(ch1c_f, zeros(size(ch1c_f)), ch3c_f) + blendGray*repmat(ch4c_f,[1 1 3]));
    RGB_234_2 = min(1, makeRGB(zeros(size(ch1c_f)), ch2c_f, ch3c_f) + blendGray*repmat(ch4c_f,[1 1 3]));

    f5 = figure('Color','w','Name', sprintf('Fig5 — 3-channel overlap (LUT) — %s', [fname ext]), ...
                'Units','normalized', 'Position',[0.03 0.03 0.94 0.92]);
    tl5 = tiledlayout(f5, 4, 5, 'Padding','compact','TileSpacing','compact');
    placeRow3(tl5, 1, RGB_123_2, norm1, norm2, norm3, overlap123, coef123, 'ch1+2+3');
    placeRow3(tl5, 2, RGB_124_2, norm1, norm2, norm4, overlap124, coef124, 'ch1+2+4');
    placeRow3(tl5, 3, RGB_134_2, norm1, norm3, norm4, overlap134, coef134, 'ch1+3+4');
    placeRow3(tl5, 4, RGB_234_2, norm2, norm3, norm4, overlap234, coef234, 'ch2+3+4');
    sgtitle(tl5, 'Fig5 — 3-channel overlap (LUT): col1=merge, col2–4=normalized, col5=overlap+coef','FontWeight','bold');

    % ---------- Fig4 (pairs) ----------
    f4 = figure('Color','w','Name', sprintf('Fig4 — Overlap analysis (LUT) — %s', [fname ext]), ...
                'Units','normalized', 'Position',[0.04 0.04 0.92 0.90]);
    tl4 = tiledlayout(f4, 6, 4, 'Padding','compact','TileSpacing','compact');
    placeRow(tl4, 1, RGB_12_2, norm1, norm2, overlap12, coef12, 'ch1+2');
    placeRow(tl4, 2, RGB_13_2, norm1, norm3, overlap13, coef13, 'ch1+3');
    placeRow(tl4, 3, RGB_23_2, norm2, norm3, overlap23, coef23, 'ch2+3');
    placeRow(tl4, 4, RGB_14_2, norm1, norm4, overlap14, coef14, 'ch1+4');
    placeRow(tl4, 5, RGB_24_2, norm2, norm4, overlap24, coef24, 'ch2+4');
    placeRow(tl4, 6, RGB_34_2, norm3, norm4, overlap34, coef34, 'ch3+4');
    sgtitle(tl4, 'Fig4 — Overlap (LUT-applied): col1=merge, col2–3=normalized, col4=overlap+coef', 'FontWeight','bold');

    
    % ---------- Return coefficients ----------
    coefStruct = struct('file',[fname ext], ...
        'coef12',coef12,'coef13',coef13,'coef23',coef23, ...
        'coef14',coef14,'coef24',coef24,'coef34',coef34, ...
        'coef123',coef123,'coef124',coef124,'coef134',coef134,'coef234',coef234);
end


% ======================= Small utilities ================================
function [ov, coef] = ovMap(A,B)
    ov   = min(A,B);
    coef = sum(ov(:));
end

function [ov, coef] = ovMap3(A,B,C)
    ov   = min(A, min(B,C));
    coef = sum(ov(:));
end

function placeRow(tl, row, mergedImg, nX, nY, ov, coef, label)
    base = (row-1)*4;
    nexttile(tl, base+1); imshow(mergedImg); title(sprintf('%s: merge (LUT)', label));
    nexttile(tl, base+2); imagesc(nX); axis image; colorbar; title(sprintf('%s: norm chX', label));
    nexttile(tl, base+3); imagesc(nY); axis image; colorbar; title(sprintf('%s: norm chY', label));
    nexttile(tl, base+4); imagesc(ov); axis image; colorbar; title(sprintf('%s: overlap (coef=%.4f)', label, coef));
end

function placeRow3(tl, row, mergedImg, nA, nB, nC, ov, coef, label)
    base = (row-1)*5;
    nexttile(tl, base+1); imshow(mergedImg); title(sprintf('%s: merge (LUT)', label));
    nexttile(tl, base+2); imagesc(nA); axis image; colorbar; title(sprintf('%s: norm chA', label));
    nexttile(tl, base+3); imagesc(nB); axis image; colorbar; title(sprintf('%s: norm chB', label));
    nexttile(tl, base+4); imagesc(nC); axis image; colorbar; title(sprintf('%s: norm chC', label));
    nexttile(tl, base+5); imagesc(ov); axis image; colorbar; title(sprintf('%s: overlap (coef=%.4f)', label, coef));
end

function savefigsafe(figHandle, outFigPath)
    try, savefig(figHandle, outFigPath); catch, warning('Could not save FIG: %s', outFigPath); end
end

function S = fillMissingFields(S, keys)
    % Ensure every struct in S has all keys (fill with NaN)
    for i = 1:numel(S)
        for j = 1:numel(keys)
            k = keys{j};
            if ~isfield(S(i),k), S(i).(k) = NaN; end
        end
    end
end
