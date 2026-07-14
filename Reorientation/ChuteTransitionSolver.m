% =========================================================================
% ChuteTransitionSolver.m  —  Pose Transition Network Solver
% =========================================================================
%
% Transition type encoding (locked convention — do not renumber; every
% place that hard-codes 1:N or switches on transType/prim refers back to
% this):
%   1 = Down-Wall
%   2 = Down-Floor
%   3 = Up-Wall
%   4 = Up-Floor
%
function ChuteTransitionSolver()

close all; clc;

THRESH_WALL    = 2.731;
THRESH_FLOOR   = 2.296;
MAX_HOP_ANGLE_DEG = 180;  
ANGLE_COST_DENOM  = 180;
wallTol       = 0.1;
planeTol      = 0.08;
thetaTol_deg  = 2;
planeMergeTol = 0.05;
omegaTol      = 0.005;
quatMatchTol  = 0.05;

[fname, fpath] = uigetfile('*.stl', 'Select STL file');
if isequal(fname, 0), disp('Cancelled.'); return; end
filePath = fullfile(fpath, fname);
[~, partName, ~] = fileparts(filePath);
fprintf('Loading: %s\n', partName);

ans_in = inputdlg({'Roll angle (deg):', 'Pitch angle (deg):'}, ...
    'Chute Orientation', 1, {'0','0'});
if isempty(ans_in), disp('Cancelled.'); return; end
chuteRoll_deg  = str2double(ans_in{1});
chutePitch_deg = str2double(ans_in{2});
fprintf('Roll=%.1f deg  Pitch=%.1f deg\n', chuteRoll_deg, chutePitch_deg);

fprintf('Running geometry pipeline...\n');
partGeometry = fegeometry(filePath);
partMesh     = generateMesh(partGeometry);
fv           = triangulation(partMesh);
partFaces    = fv.ConnectivityList;
partVertices = fv.Points;

partSpan  = max(partVertices) - min(partVertices);
partScale = max(partSpan);
hTol_part = partScale * 0.03;

convexHullFaces = convhull(partVertices, 'Simplify', true);
chullVertexIdx  = unique(convexHullFaces(:));
centroidCoords  = centroidOfPolyhedron(partVertices, partFaces);
fprintf('  Centroid: [%.3f %.3f %.3f]\n', centroidCoords);

[restingPlaneVerts, restingPlaneEqs] = ...
    findFloorPlanes(partVertices, convexHullFaces, planeTol);
fprintf('  Floor planes: %d\n', numel(restingPlaneVerts));

roll  = -deg2rad(chuteRoll_deg);
pitch =  deg2rad(chutePitch_deg);
q_roll  = q_fromAxisAngle([1,0,0], roll);
q_pitch = q_fromAxisAngle([0,1,0], pitch);
q_chute = q_compose(q_roll, q_pitch);
Rchute  = q_toRotm(q_chute);

gravityDir  = (Rchute' * [0;0;-1])';
slideDir    = normaliseVec((Rchute * [1;0;0])');
floorNorm_c = normaliseVec((Rchute * [0;0;1])');
wallNorm_c  = normaliseVec((Rchute * [0;1;0])');
gWorld      = [0;0;-1];

[rawCandidates, planeQuats] = thetaFromHullEdges( ...
    partVertices, convexHullFaces, chullVertexIdx, ...
    restingPlaneVerts, restingPlaneEqs, ...
    centroidCoords, wallTol, planeTol);
fprintf('  Raw candidates: %d\n', numel(rawCandidates));

allPoses = mergeByTheta(rawCandidates, thetaTol_deg, planeMergeTol, partVertices);
fprintf('  After theta merge: %d\n', numel(allPoses));
clear rawCandidates;

gd = gravityDir(:)/norm(gravityDir);
if abs(gd(1)) < 0.9, tmp = [1;0;0]; else, tmp = [0;1;0]; end
e1g = cross(gd,tmp); e1g=e1g/norm(e1g);
e2g = cross(gd,e1g); e2g=e2g/norm(e2g);

numAll = numel(allPoses);
stableMask = false(numAll,1);
for ci = 1:numAll
    pos   = allPoses(ci);
    vW    = pos.verticesWorld;
    cW    = pos.centroidWorld(:)';
    chV   = vW(chullVertexIdx,:);
    pts2D = [chV*e1g, chV*e2g];
    c2D   = [cW*e1g, cW*e2g];
    try
        hIdx = convhull(pts2D(:,1),pts2D(:,2));
        stableMask(ci) = inpolygon(c2D(1),c2D(2),pts2D(hIdx,1),pts2D(hIdx,2));
    catch; stableMask(ci) = false; end
end

stableIdx = find(stableMask);
numStable = numel(stableIdx);
fprintf('  Stable poses: %d\n', numStable);
if numStable == 0, errordlg('No stable poses found.','Error'); return; end

secUnstable = false(numStable,1);
for si = 1:numStable
    pos   = allPoses(stableIdx(si));
    C0    = pos.centroidWorld(:)';
    vW    = (Rchute*(pos.verticesWorld'-C0'))'+C0;
    supIdx = union(pos.floorContactVertIdx, pos.wallContactVertIdx);
    supV   = vW(supIdx,:);
    if size(supV,1) < 3, secUnstable(si)=true; continue; end
    gC = supV*gWorld; mG = min(gC); pp = mG*gWorld';
    sProj = zeros(size(supV));
    for kk=1:size(supV,1)
        p_=supV(kk,:); d_=dot((p_-pp),gWorld');
        sProj(kk,:)=p_-d_*gWorld';
    end
    sProj = uniquetol(sProj,1e-5,'ByRows',true);
    if size(sProj,1)<3, secUnstable(si)=true; continue; end
    if abs(gWorld(1))<0.9, t2=[1;0;0]; else, t2=[0;1;0]; end
    e1b=cross(gWorld,t2); e1b=e1b/norm(e1b);
    e2b=cross(gWorld,e1b); e2b=e2b/norm(e2b);
    s2D=[sProj*e1b, sProj*e2b];
    cW2=[dot(C0,e1b),dot(C0,e2b)];
    try
        hI=convhull(s2D(:,1),s2D(:,2));
        if ~inpolygon(cW2(1),cW2(2),s2D(hI,1),s2D(hI,2))
            secUnstable(si)=true;
        end
    catch; secUnstable(si)=true; end
end

keepMask = ~secUnstable;
stableIdx = stableIdx(keepMask);
numStable = numel(stableIdx);
fprintf('  After secondary stability: %d\n', numStable);
if numStable == 0, errordlg('No poses passed secondary stability.','Error'); return; end

[Qs, omegas, heights, ~, ~] = computeCSA( ...
    allPoses, stableIdx, chullVertexIdx, Rchute, gWorld);

[stableIdx, Qs, omegas, heights, mergedGroups] = mergeByCSAValues( ...
    stableIdx, Qs, omegas, heights, omegaTol, hTol_part, allPoses, Rchute);
numStable = numel(stableIdx);
fprintf('  After CSA merge: %d poses\n', numStable);

[ratioWall, ratioFloor, ~, ~] = computeTransitionRatios( ...
    allPoses, stableIdx, Rchute, slideDir, floorNorm_c, wallNorm_c, planeTol);

for si = 1:numStable
    if ratioFloor(si)   >= THRESH_FLOOR   || ...
       ratioWall(si)    >= THRESH_WALL
        Qs(si) = 0;
    end
end

refQuats    = zeros(numStable,4);
refQuatsAll = cell(numStable,1);
for si = 1:numStable
    pos     = allPoses(stableIdx(si));
    q_align = planeQuats{pos.floorPlaneIdx};
    q_theta = q_fromAxisAngle([0,0,1], pos.theta);
    q_ref   = q_compose(q_align, q_theta);
    [~,mi]  = max(abs(q_ref)); if q_ref(mi)<0, q_ref=-q_ref; end
    refQuats(si,:) = q_ref;

    grpIdx = mergedGroups{si}(:)';
    qList  = zeros(numel(grpIdx),4);
    for gk = 1:numel(grpIdx)
        pg   = allPoses(grpIdx(gk));
        qa   = planeQuats{pg.floorPlaneIdx};
        qt   = q_fromAxisAngle([0,0,1], pg.theta);
        qg   = q_compose(qa,qt);
        [~,mi]=max(abs(qg)); if qg(mi)<0, qg=-qg; end
        qList(gk,:) = qg;
    end
    refQuatsAll{si} = qList;
end

% =========================================================================
% BUILD TRANSITION GRAPH
% =========================================================================
fprintf('Building transition graph...\n');

% edges cols: [src, dst, cost, type, angleChangeDeg]
edges = zeros(0,5);

maxChain  = numStable + 2;
validPose = Qs > 0;

for si = 1:numStable
    for transType = 1:4
        angleCap = MAX_HOP_ANGLE_DEG;   % same 180° single-hop budget for all
                                         % four transition types (Down-Wall,
                                         % Down-Floor, Up-Wall, Up-Floor) —
                                         % a chain whose cumulative rotation
                                         % would exceed this before reaching
                                         % the next stable face is not a
                                         % valid single hop and gets no edge;
                                         % it must be split into two (or
                                         % more) separate transitions.
        [chainNodes, firstCost, tType, angleChangeDeg] = resolveTransitionChain(si, transType, ...
            allPoses, stableIdx, refQuats, refQuatsAll, ...
            Qs, validPose, ...
            quatMatchTol, planeTol, wallTol, maxChain, ...
            restingPlaneVerts, restingPlaneEqs, chullVertexIdx, ...
            centroidCoords, planeQuats, angleCap, ANGLE_COST_DENOM);

        if numel(chainNodes) < 2, continue; end

        % NOTE: the max-transition-angle cap is now applied dynamically in
        % the UI via the "Max Δ°" field, not hard-coded here — all edges
        % are built (regardless of angle) and filtered live at query time
        % so the cap can be adjusted without rebuilding the graph.
        typeLabel = tType * 10;
        for k = 1:numel(chainNodes)-1
            s = chainNodes(k);
            d = chainNodes(k+1);
            if k == 1
                c = firstCost;
                t = transType;
                aDeg = angleChangeDeg;
            else
                c = 0;
                t = typeLabel;
                aDeg = 0;
            end
            dup = any(edges(:,1)==s & edges(:,2)==d & edges(:,4)==t);
            if ~dup
                edges(end+1,:) = [s, d, c, t, aDeg]; %#ok<AGROW>
            end
        end

        tNames = {'Down-Wall','Down-Floor','Up-Wall','Up-Floor'};
        nodeStr = sprintf('%d', chainNodes(1));
        for k=2:numel(chainNodes)
            nodeStr = [nodeStr sprintf('->%d', chainNodes(k))]; %#ok<AGROW>
        end
        fprintf('  %s: %s  (cost=%.3f, angle=%.1f deg)\n', tNames{transType}, nodeStr, firstCost, angleChangeDeg);
    end
end

fprintf('Total edges: %d\n', size(edges,1));

launchUI(numStable, stableIdx, allPoses, Qs, ratioWall, ratioFloor, ...
    THRESH_WALL, THRESH_FLOOR, edges, partName, ...
    chuteRoll_deg, chutePitch_deg, validPose, Rchute, partFaces);

end % ChuteTransitionSolver


% =========================================================================
%  TRANSITION CATEGORY  (Wall vs Floor, direction-agnostic)
% =========================================================================
%  primType: 1=Down-Wall, 2=Down-Floor, 3=Up-Wall, 4=Up-Floor
%  category: 1=Wall (primType 1 or 3), 2=Floor (primType 2 or 4)
%
%  This is the single source of truth for "a Wall transition is a Wall
%  transition whether it's up-chute or down-chute" — every place that
%  needs to treat Up/Down variants as interchangeable goes through this
%  function instead of re-deriving the grouping locally.
function c = catOf(primType)
    if primType == 1 || primType == 3
        c = 1;   % Wall
    else
        c = 2;   % Floor
    end
end

function s = catLabel(primType)
    if catOf(primType) == 1, s = 'Wall'; else, s = 'Floor'; end
end


% =========================================================================
%  INTERACTIVE UI
% =========================================================================
function launchUI(numStable, stableIdx, allPoses, Qs, ratioWall, ratioFloor, ...
    THRESH_WALL, THRESH_FLOOR, edges, partName, rollDeg, pitchDeg, validPose, Rchute, partFaces)

% ── state ──────────────────────────────────────────────────────────────
state.selA     = 0;
state.selB     = 0;
state.mode     = 'idle';
state.pathEdgesS = [];
state.pathEdgesC = [];
state.seqEdgeLists = {};
state.globalSeqMode = false;

% ── layout ─────────────────────────────────────────────────────────────
cols  = ceil(sqrt(numStable));
rows  = ceil(numStable / cols);
bRad  = 0.35;
xPos  = zeros(numStable,1);
yPos  = zeros(numStable,1);
for si = 1:numStable
    c = mod(si-1, cols);
    r = floor((si-1)/cols);
    xPos(si) = c * 2.5 + 1.25;
    yPos(si) = (rows - r - 1) * 2.5 + 1.25;
end

% ── figure ─────────────────────────────────────────────────────────────
hFig = figure('Name', sprintf('Transition Network — %s', partName), ...
    'Color', [0.11 0.12 0.15], ...
    'Units', 'normalized', 'Position', [0 0 1 1], ...
    'WindowState', 'maximized', ...
    'Resize', 'on', 'NumberTitle', 'off');

% ── Close button ───────────────────────────────────────────────────────
uicontrol(hFig, 'Style', 'pushbutton', 'String', '✕  Close', ...
    'Units', 'normalized', 'Position', [0.938, 0.958, 0.058, 0.033], ...
    'BackgroundColor', [0.50 0.13 0.13], 'ForegroundColor', [1.0 0.82 0.82], ...
    'FontSize', 9, 'FontWeight', 'bold', 'Callback', @(~,~) close(hFig));

% ── Clear button ───────────────────────────────────────────────────────
uicontrol(hFig, 'Style', 'pushbutton', 'String', '↺  Clear', ...
    'Units', 'normalized', 'Position', [0.876, 0.958, 0.058, 0.033], ...
    'BackgroundColor', [0.16 0.20 0.28], 'ForegroundColor', [0.70 0.80 1.00], ...
    'FontSize', 9, 'FontWeight', 'bold', 'Callback', @onClear);

% ── Pose viewer A — START (left top) ───────────────────────────────────
axPoseA = axes(hFig, 'Units', 'normalized', 'Position', [0.005, 0.53, 0.185, 0.40], ...
    'Color', [0.09 0.10 0.13], 'XColor','none','YColor','none','ZColor','none');
title(axPoseA, 'START pose', 'Color', [1 0.82 0.14], 'FontSize', 8);

% ── Pose viewer B — END (left bottom) ─────────────────────────────────
axPoseB = axes(hFig, 'Units', 'normalized', 'Position', [0.005, 0.10, 0.185, 0.40], ...
    'Color', [0.09 0.10 0.13], 'XColor','none','YColor','none','ZColor','none');
title(axPoseB, 'END pose', 'Color', [0.18 0.88 0.52], 'FontSize', 8);

% ── Network axes (centre) ──────────────────────────────────────────────
axNet = axes(hFig, 'Units', 'normalized', 'Position', [0.20, 0.19, 0.50, 0.75], ...
    'Color', [0.11 0.12 0.15], 'XColor','none','YColor','none');
hold(axNet,'on');
xlim(axNet, [0, cols*2.5 + 0.5]);
ylim(axNet, [0, rows*2.5 + 0.5]);
axis(axNet, 'equal');
title(axNet, sprintf('%s  |  Roll=%.0f°  Pitch=%.0f°', partName, rollDeg, pitchDeg), ...
    'Color','w','FontSize',11,'FontWeight','bold');

% ── Info panel (right) — scrollable listbox ────────────────────────────
hInfoList = uicontrol(hFig, 'Style', 'listbox', ...
    'Units', 'normalized', 'Position', [0.720, 0.19, 0.275, 0.75], ...
    'BackgroundColor', [0.08 0.09 0.12], ...
    'ForegroundColor', [0.82 0.84 0.92], ...
    'FontName', 'Courier New', 'FontSize', 8.5, ...
    'String', {'  — select poses or run sequence —'}, ...
    'Max', 2, 'Min', 0, ...
    'Value', [], ...
    'HorizontalAlignment', 'left', ...
    'Enable', 'inactive');

% ── Status bar ─────────────────────────────────────────────────────────
axStatus = axes(hFig, 'Units','normalized', 'Position',[0.005, 0.945, 0.865, 0.045], ...
    'Color',[0.07 0.08 0.11], 'XColor','none','YColor','none');
xlim(axStatus,[0,1]); ylim(axStatus,[0,1]);
hStatusTxt = text(axStatus, 0.012, 0.50, ...
    'Click any pose to select START', ...
    'Color',[0.80 0.82 0.88], 'FontSize',9.5, ...
    'VerticalAlignment','middle','Interpreter','none');

% ── Bottom control strip ───────────────────────────────────────────────
uipanel(hFig, 'Units','normalized', 'Position',[0.005, 0.005, 0.990, 0.115], ...
    'BackgroundColor',[0.08 0.09 0.12], 'ForegroundColor',[0.45 0.50 0.60], ...
    'Title','Sequence Constraints', 'FontSize',8, 'FontWeight','bold', 'BorderType','line');

seqOpts    = {'—','Wall','Floor'};
numSeqSlots = 6;
hSeqDD     = gobjects(numSeqSlots,1);

ddW  = 0.072; ddH = 0.50; ddBot = 0.08;
startX = 0.01; stepX = ddW + 0.010;

for k = 1:numSeqSlots
    annotation(hFig,'textbox', ...
        [0.005+startX+(k-1)*stepX*0.990, 0.005+ddBot*0.115+ddH*0.115+0.005, ddW*0.990, 0.013], ...
        'String',sprintf('Slot %d',k),'Color',[0.48 0.52 0.62],'FontSize',6.5, ...
        'EdgeColor','none','BackgroundColor','none','HorizontalAlignment','center');
    hSeqDD(k) = uicontrol(hFig,'Style','popupmenu','String',seqOpts, ...
        'Units','normalized', ...
        'Position',[0.005+startX+(k-1)*stepX*0.990, 0.005+ddBot*0.115, ddW*0.990, ddH*0.115], ...
        'BackgroundColor',[0.15 0.16 0.21],'ForegroundColor',[0.88 0.88 0.94], ...
        'FontSize',8.5,'Callback',@onSeqChange);
end

metX = startX + numSeqSlots*stepX + 0.015;
annotation(hFig,'textbox', ...
    [0.005+metX*0.990, 0.005+ddBot*0.115+ddH*0.115+0.005, 0.078*0.990, 0.013], ...
    'String','Metric','Color',[0.48 0.52 0.62],'FontSize',6.5, ...
    'EdgeColor','none','BackgroundColor','none','HorizontalAlignment','center');
hMetricDD = uicontrol(hFig,'Style','popupmenu','String',{'Lowest Cost','Fewest Hops'}, ...
    'Units','normalized', ...
    'Position',[0.005+metX*0.990, 0.005+ddBot*0.115, 0.078*0.990, ddH*0.115], ...
    'BackgroundColor',[0.15 0.16 0.21],'ForegroundColor',[0.88 0.88 0.94], ...
    'FontSize',8.5,'Callback',@onSeqChange);

genX = metX + 0.085;
annotation(hFig,'textbox', ...
    [0.005+genX*0.990, 0.005+ddBot*0.115+ddH*0.115+0.005, 0.100*0.990, 0.013], ...
    'String','Auto-find','Color',[0.48 0.52 0.62],'FontSize',6.5, ...
    'EdgeColor','none','BackgroundColor','none','HorizontalAlignment','center');
uicontrol(hFig,'Style','pushbutton','String','⚡ Generate Sequences', ...
    'Units','normalized', ...
    'Position',[0.005+genX*0.990, 0.005+ddBot*0.115, 0.100*0.990, ddH*0.115], ...
    'BackgroundColor',[0.14 0.25 0.18],'ForegroundColor',[0.45 0.92 0.60], ...
    'FontSize',8,'FontWeight','bold','Callback',@onGenerateSeqs);

% ── Global / Individual toggle ─────────────────────────────────────────
togX = genX + 0.107;
annotation(hFig,'textbox', ...
    [0.005+togX*0.990, 0.005+ddBot*0.115+ddH*0.115+0.005, 0.095*0.990, 0.013], ...
    'String','Seq Scope','Color',[0.48 0.52 0.62],'FontSize',6.5, ...
    'EdgeColor','none','BackgroundColor','none','HorizontalAlignment','center');
hScopeTog = uicontrol(hFig,'Style','pushbutton', ...
    'String','● Individual', ...
    'Units','normalized', ...
    'Position',[0.005+togX*0.990, 0.005+ddBot*0.115, 0.095*0.990, ddH*0.115], ...
    'BackgroundColor',[0.12 0.22 0.35],'ForegroundColor',[0.45 0.75 1.00], ...
    'FontSize',8,'FontWeight','bold','Callback',@onScopeToggle);

% ── Direction dropdown ─────────────────────────────────────────────────
dirX = togX + 0.102;
annotation(hFig,'textbox', ...
    [0.005+dirX*0.990, 0.005+ddBot*0.115+ddH*0.115+0.005, 0.090*0.990, 0.013], ...
    'String','Direction','Color',[0.48 0.52 0.62],'FontSize',6.5, ...
    'EdgeColor','none','BackgroundColor','none','HorizontalAlignment','center');
hDirDD = uicontrol(hFig,'Style','popupmenu','String',{'Paths → END','Paths FROM START'}, ...
    'Units','normalized', ...
    'Position',[0.005+dirX*0.990, 0.005+ddBot*0.115, 0.090*0.990, ddH*0.115], ...
    'BackgroundColor',[0.15 0.16 0.21],'ForegroundColor',[0.88 0.88 0.94], ...
    'FontSize',8.5,'Callback',@onSeqChange);

% ── Transition type filter (multi-select listbox) ──────────────────────
% Placed to the right of the Direction dropdown
filtX = dirX + 0.097;
annotation(hFig,'textbox', ...
    [0.005+filtX*0.990, 0.005+ddBot*0.115+ddH*0.115+0.005, 0.088*0.990, 0.013], ...
    'String','Allow Types','Color',[0.48 0.52 0.62],'FontSize',6.5, ...
    'EdgeColor','none','BackgroundColor','none','HorizontalAlignment','center');
hTypeFilter = uicontrol(hFig,'Style','listbox', ...
    'String',{'Down-Wall','Down-Floor','Up-Wall','Up-Floor'}, ...
    'Units','normalized', ...
    'Position',[0.005+filtX*0.990, 0.005+ddBot*0.115, 0.088*0.990, ddH*0.115*1.25], ...
    'BackgroundColor',[0.13 0.15 0.20],'ForegroundColor',[0.88 0.92 1.00], ...
    'FontSize',8, ...
    'Max',4,'Min',0, ...          % multi-select
    'Value',[1 2 3 4], ...        % all selected by default
    'Callback',@onSeqChange);

% ── Max transition angle cap (deg) ──────────────────────────────────────
angX = filtX + 0.093;
annotation(hFig,'textbox', ...
    [0.005+angX*0.990, 0.005+ddBot*0.115+ddH*0.115+0.005, 0.062*0.990, 0.013], ...
    'String','Max Δ°','Color',[0.48 0.52 0.62],'FontSize',6.5, ...
    'EdgeColor','none','BackgroundColor','none','HorizontalAlignment','center');
hMaxAngleEdit = uicontrol(hFig,'Style','edit','String','180', ...
    'Units','normalized', ...
    'Position',[0.005+angX*0.990, 0.005+ddBot*0.115, 0.062*0.990, ddH*0.115], ...
    'BackgroundColor',[0.15 0.16 0.21],'ForegroundColor',[0.95 0.85 0.30], ...
    'FontSize',8.5,'FontWeight','bold','Callback',@onSeqChange);

% ── static background edges ────────────────────────────────────────────
for ei = 1:size(edges,1)
    s = edges(ei,1); d = edges(ei,2); t = edges(ei,4); aDeg = edges(ei,5);
    if     t == 1,  col = [0.28 0.52 0.88]; ls = '-';   % Down-Wall  (blue)
    elseif t == 2,  col = [0.28 0.76 0.52]; ls = '-';   % Down-Floor (green)
    elseif t == 3,  col = [0.90 0.55 0.15]; ls = '-';   % Up-Wall    (orange)
    elseif t == 4,  col = [0.90 0.30 0.30]; ls = '-';   % Up-Floor   (red-orange)
    elseif t == 10, col = [0.28 0.52 0.88]; ls = '--';
    elseif t == 20, col = [0.28 0.76 0.52]; ls = '--';
    elseif t == 30, col = [0.90 0.55 0.15]; ls = '--';
    else,           col = [0.90 0.30 0.30]; ls = '--';
    end
    bgLbl = [];
    if aDeg > 0.5, bgLbl = sprintf('%.0f°', aDeg); end
    drawArrow(axNet, xPos(s),yPos(s), xPos(d),yPos(d), bRad, col, 0.5, 0.35, bgLbl, ls);
end

% ── pose bubbles ───────────────────────────────────────────────────────
hBubble    = gobjects(numStable,1);
hBubbleTxt = gobjects(numStable,1);
hBubbleSub = gobjects(numStable,1);

for si = 1:numStable
    pos = allPoses(stableIdx(si));
    [faceC, edgeC] = bubbleColors(si, ratioWall, ratioFloor, ...
        THRESH_WALL, THRESH_FLOOR, Qs);
    th   = linspace(0,2*pi,64);
    xCir = xPos(si) + bRad*cos(th);
    yCir = yPos(si) + bRad*sin(th);
    hBubble(si) = fill(axNet, xCir, yCir, faceC, 'EdgeColor',edgeC,'LineWidth',2.0, ...
        'ButtonDownFcn', @(~,~) onBubbleClick(si));
    hBubbleTxt(si) = text(axNet, xPos(si), yPos(si)+0.08, sprintf('%d',si), ...
        'Color','w','FontSize',11,'FontWeight','bold', ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'ButtonDownFcn', @(~,~) onBubbleClick(si));
    subStr = sprintf('P%d θ%d°', pos.floorPlaneIdx, round(rad2deg(pos.theta)));
    hBubbleSub(si) = text(axNet, xPos(si), yPos(si)-0.15, subStr, ...
        'Color',[0.55 0.58 0.68],'FontSize',6.5, ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'ButtonDownFcn', @(~,~) onBubbleClick(si));
end

pathHandles = {};

% ======================================================================
%  CALLBACKS
% ======================================================================

    % ------------------------------------------------------------------
    %  Get currently allowed transition types from filter listbox
    %  Returns a logical [allowDnWall, allowDnFloor, allowUpWall, allowUpFloor]
    % ------------------------------------------------------------------
    function allowed = getAllowedTypes()
        sel = hTypeFilter.Value;   % indices of selected items (1=Down-Wall,2=Down-Floor,3=Up-Wall,4=Up-Floor)
        allowed = [any(sel==1), any(sel==2), any(sel==3), any(sel==4)];
    end

    % ------------------------------------------------------------------
    %  Filter edges by allowed types
    % ------------------------------------------------------------------
    function eFiltered = filterEdgesByType(edgesIn, allowed)
        keep = true(size(edgesIn,1),1);
        for ei2 = 1:size(edgesIn,1)
            t = edgesIn(ei2,4);
            prim = mod(t-1,10)+1;   % 1=dnWall,2=dnFloor,3=upWall,4=upFloor
            if prim == 1 && ~allowed(1), keep(ei2) = false; end
            if prim == 2 && ~allowed(2), keep(ei2) = false; end
            if prim == 3 && ~allowed(3), keep(ei2) = false; end
            if prim == 4 && ~allowed(4), keep(ei2) = false; end
        end
        eFiltered = edgesIn(keep,:);
        % Re-index: edge src/dst are pose indices (1..numStable), not edge indices,
        % so no remapping needed — just subset the rows.
    end

    % ------------------------------------------------------------------
    %  Get currently set maximum transition angle (deg) from the edit box
    % ------------------------------------------------------------------
    function maxDeg = getMaxAngle()
        v = str2double(hMaxAngleEdit.String);
        if isnan(v) || v <= 0
            v = 180;
            hMaxAngleEdit.String = '180';
        end
        maxDeg = v;
    end

    % ------------------------------------------------------------------
    %  Filter edges by maximum transition angle. Only the "first hop" of
    %  a transition (type < 10) carries a nonzero angleChangeDeg; chained
    %  free hops (type >= 10, cost 0) carry aDeg=0 and are never capped
    %  here, since the cap targets the underlying physical rotation step.
    %  Any edge exceeding the cap is dropped, so a path that would need
    %  it simply becomes unreachable (invalid).
    % ------------------------------------------------------------------
    function eFiltered = filterEdgesByAngle(edgesIn, maxDeg)
        keep = edgesIn(:,5) <= maxDeg;
        eFiltered = edgesIn(keep,:);
    end

    % ------------------------------------------------------------------
    %  SCOPE TOGGLE: Individual <-> Global
    % ------------------------------------------------------------------
    function onScopeToggle(~,~)
        state.globalSeqMode = ~state.globalSeqMode;
        if state.globalSeqMode
            hScopeTog.String = '◉ Global';
            hScopeTog.BackgroundColor = [0.20 0.14 0.35];
            hScopeTog.ForegroundColor = [0.80 0.55 1.00];
        else
            hScopeTog.String = '● Individual';
            hScopeTog.BackgroundColor = [0.12 0.22 0.35];
            hScopeTog.ForegroundColor = [0.45 0.75 1.00];
        end
        if isSeqActive()
            onSeqChange([], []);
        end
    end

    % ------------------------------------------------------------------
    %  CLEAR
    % ------------------------------------------------------------------
    function onClear(~,~)
        clearPathHandles();
        resetAllBubbles();
        clearInfoPanel();
        state.selA = 0; state.selB = 0;
        state.mode = 'idle';
        state.pathEdgesS = []; state.pathEdgesC = [];
        state.seqEdgeLists = {};
        cla(axPoseA); title(axPoseA,'START pose','Color',[1 0.82 0.14],'FontSize',8);
        cla(axPoseB); title(axPoseB,'END pose','Color',[0.18 0.88 0.52],'FontSize',8);
        for k2 = 1:numSeqSlots, hSeqDD(k2).Value = 1; end
        hMaxAngleEdit.String = '180';
        hStatusTxt.String = 'Click any pose to select START';
    end

    % ------------------------------------------------------------------
    %  BUBBLE CLICK
    % ------------------------------------------------------------------
    function onBubbleClick(si)
        if state.selA == 0
            state.selA = si; state.selB = 0;
            state.mode = 'idle';
            clearPathHandles(); resetAllBubbles(); clearInfoPanel();
            cla(axPoseB); title(axPoseB,'END pose','Color',[0.18 0.88 0.52],'FontSize',8);
            highlightBubble(si, [0.82 0.65 0.08], [1 0.88 0.28]);
            renderPoseAxes(axPoseA, allPoses(stableIdx(si)), si, 'START', Rchute, partFaces);
            if isSeqActive()
                hStatusTxt.String = sprintf('START = %d  |  Sequence active — showing reachable END poses', si);
                drawnow; runSeqFromStart();
            else
                hStatusTxt.String = sprintf('START = %d  |  Click another pose to set END, or use sequence controls below', si);
            end

        elseif si == state.selA
            onClear();

        else
            if isSeqActive()
                if hDirDD.Value == 1
                    state.selB = si;
                    clearPathHandles(); resetAllBubbles();
                    highlightBubble(state.selA, [0.82 0.65 0.08], [1 0.88 0.28]);
                    highlightBubble(si, [0.12 0.72 0.32], [0.28 0.92 0.52]);
                    renderPoseAxes(axPoseB, allPoses(stableIdx(si)), si, 'END', Rchute, partFaces);
                    hStatusTxt.String = sprintf('START = %d  END = %d  |  Showing sequence path(s)', state.selA, si);
                    drawnow; runSeqBetween(state.selA, si);
                else
                    state.selA = si; state.selB = 0;
                    clearPathHandles(); resetAllBubbles(); clearInfoPanel();
                    cla(axPoseB); title(axPoseB,'END pose','Color',[0.18 0.88 0.52],'FontSize',8);
                    highlightBubble(si, [0.82 0.65 0.08], [1 0.88 0.28]);
                    renderPoseAxes(axPoseA, allPoses(stableIdx(si)), si, 'START', Rchute, partFaces);
                    drawnow; runSeqFromStart();
                end
            else
                state.selB = si;
                clearPathHandles(); resetAllBubbles();
                highlightBubble(state.selA, [0.82 0.65 0.08], [1 0.88 0.28]);
                highlightBubble(si, [0.12 0.72 0.32], [0.28 0.92 0.52]);
                renderPoseAxes(axPoseB, allPoses(stableIdx(si)), si, 'END', Rchute, partFaces);
                state.mode = 'path';
                hStatusTxt.String = sprintf('START=%d  END=%d  |  Computing...', state.selA, si);
                drawnow; computeAndShowPaths();
            end
        end
    end

    % ------------------------------------------------------------------
    %  SEQUENCE DROPDOWN CHANGED
    % ------------------------------------------------------------------
    function onSeqChange(~,~)
        if state.selA == 0 && ~state.globalSeqMode
            hStatusTxt.String = 'Click any pose to select START first (or enable Global scope)';
            return;
        end
        clearPathHandles(); resetAllBubbles(); clearInfoPanel();
        if state.selA > 0
            highlightBubble(state.selA, [0.82 0.65 0.08], [1 0.88 0.28]);
        end
        if isSeqActive()
            if state.globalSeqMode
                runSeqGlobal();
            elseif state.selB > 0 && hDirDD.Value == 1
                runSeqBetween(state.selA, state.selB);
            else
                runSeqFromStart();
            end
        else
            if state.selB > 0 && ~state.globalSeqMode
                state.mode = 'path';
                computeAndShowPaths();
            else
                if state.globalSeqMode
                    hStatusTxt.String = 'Set sequence slots to run global analysis';
                else
                    hStatusTxt.String = sprintf('START = %d  |  Select END pose or set sequence slots below', state.selA);
                end
            end
        end
    end

    % ------------------------------------------------------------------
    %  GENERATE SEQUENCES (auto-find best)
    % ------------------------------------------------------------------
    function onGenerateSeqs(~,~)
        if state.selA == 0 && ~state.globalSeqMode
            hStatusTxt.String = 'Select a START pose first (or switch to Global scope)';
            return;
        end
        hStatusTxt.String = 'Analysing sequences... please wait';
        drawnow;
        useCost = (hMetricDD.Value == 1);
        anchor  = state.selA;
        allowed = getAllowedTypes();
        maxDeg  = getMaxAngle();
        eF      = filterEdgesByType(edges, allowed);
        eF      = filterEdgesByAngle(eF, maxDeg);

        if state.globalSeqMode
            results = computeOptimalSequencesGlobal(numStable, eF, validPose, useCost, maxDeg);
        else
            results = computeOptimalSequencesIndividual(anchor, numStable, eF, validPose, useCost, hDirDD.Value==1, maxDeg);
        end

        if isempty(results)
            hStatusTxt.String = 'No valid sequences found for this network';
            return;
        end

        bestSeq = results(1).seq;
        for k2 = 1:numSeqSlots
            if k2 <= numel(bestSeq)
                hSeqDD(k2).Value = bestSeq(k2) + 1;
            else
                hSeqDD(k2).Value = 1;
            end
        end
        if isfield(results(1),'direction')
            if strcmp(results(1).direction, '→ END')
                hDirDD.Value = 1;
            else
                hDirDD.Value = 2;
            end
        end

        clearPathHandles(); resetAllBubbles(); clearInfoPanel();
        if state.selA > 0
            highlightBubble(state.selA, [0.82 0.65 0.08], [1 0.88 0.28]);
        end
        drawnow;

        if state.globalSeqMode
            runSeqGlobal();
        else
            runSeqFromStart();
        end

        metStr = 'cost'; if ~useCost, metStr = 'hops'; end
        scopeStr = 'Individual'; if state.globalSeqMode, scopeStr = 'Global'; end
        hStatusTxt.String = sprintf( ...
            '[%s] Best seq: %s  |  Coverage=%d  Avg %s=%.3f  |  Adjust slots or click a pose', ...
            scopeStr, seqToString(bestSeq), results(1).coverage, metStr, results(1).avgMetric);
    end

    % ------------------------------------------------------------------
    %  NORMAL PATH MODE
    % ------------------------------------------------------------------
    function computeAndShowPaths()
        A = state.selA; B = state.selB;
        allowed = getAllowedTypes();
        maxDeg  = getMaxAngle();
        eF = filterEdgesByType(edges, allowed);
        eF = filterEdgesByAngle(eF, maxDeg);

        [~, edgesS] = bfsShortestPath(A, B, numStable, eF);
        [~, edgesC] = dijkstraCheapestPath(A, B, numStable, eF);

        % Map filtered-edge indices back to original edge indices for drawing
        % (eF rows correspond to original edge rows; we need original indices
        %  so drawArrow can use edges(:,4) type etc. — but since we pass eF
        %  into the path functions, the returned indices are into eF, not edges.
        %  Re-map them here.)
        state.pathEdgesS = edgesS;
        state.pathEdgesC = edgesC;

        % Check if both paths are identical
        samePath = isequal(sort(edgesS), sort(edgesC));

        clearPathHandles();
        if samePath && ~isempty(edgesS)
            % Draw single path labelled "Shortest & Lowest Cost"
            for k = 1:numel(edgesS)
                ei = edgesS(k);
                s=eF(ei,1); d=eF(ei,2); r=eF(ei,3); t=eF(ei,4); aDeg=eF(ei,5);
                ls='-'; if t>=10, ls='--'; end
                h = drawArrow(axNet,xPos(s),yPos(s),xPos(d),yPos(d), ...
                    bRad,[0.92 0.72 0.08],2.8,1.0,buildEdgeLabel(r,aDeg),ls);
                pathHandles{end+1}=h;
            end
        else
            for k = 1:numel(edgesS)
                ei = edgesS(k);
                s=eF(ei,1); d=eF(ei,2); r=eF(ei,3); t=eF(ei,4); aDeg=eF(ei,5);
                ls='-'; if t>=10, ls='--'; end
                h = drawArrow(axNet,xPos(s),yPos(s),xPos(d),yPos(d), ...
                    bRad,[1.0 0.78 0.08],2.6,1.0,buildEdgeLabel(r,aDeg),ls);
                pathHandles{end+1}=h;
            end
            for k = 1:numel(edgesC)
                ei = edgesC(k);
                s=eF(ei,1); d=eF(ei,2); r=eF(ei,3); t=eF(ei,4); aDeg=eF(ei,5);
                ls='-'; if t>=10, ls='--'; end
                h = drawArrow(axNet,xPos(s),yPos(s),xPos(d),yPos(d), ...
                    bRad*0.82,[0.18 0.88 0.88],2.6,1.0,buildEdgeLabel(r,aDeg),ls);
                pathHandles{end+1}=h;
            end
        end

        % Pass eF to info panel so it uses the right edge table
        updatePathInfoPanel(eF, samePath);

        if isempty(edgesS)
            hStatusTxt.String = sprintf('No path: %d → %d  |  Click any pose to reset', A, B);
        elseif samePath
            costS = sum(eF(edgesS,3)); hopsS = countRealHops(edgesS,eF);
            hStatusTxt.String = sprintf( ...
                'START=%d → END=%d  |  ■Gold: Shortest & Lowest Cost  %d hop(s) %.3f  |  Click any pose to reset', ...
                A,B,hopsS,costS);
        else
            costS = sum(eF(edgesS,3)); hopsS = countRealHops(edgesS,eF);
            costC = sum(eF(edgesC,3)); hopsC = countRealHops(edgesC,eF);
            hStatusTxt.String = sprintf( ...
                'START=%d → END=%d  |  ■Gold: %d hop(s) %.3f  |  ■Cyan: %d hop(s) %.3f  |  Click any pose to reset', ...
                A,B,hopsS,costS,hopsC,costC);
        end
    end

    % ------------------------------------------------------------------
    %  SEQUENCE SEARCH: from START (individual mode)
    %  Also used for "→ END" direction by searching in the reversed graph
    %  and then flipping results — mirrors FROM START logic exactly.
    % ------------------------------------------------------------------
    function runSeqFromStart()
        state.mode = 'seq';
        seq     = getSequence();
        useCost = (hMetricDD.Value == 1);
        isToEnd = (hDirDD.Value == 1);
        anchor  = state.selA;
        allowed = getAllowedTypes();
        maxDeg  = getMaxAngle();
        eF      = filterEdgesByType(edges, allowed);
        eF      = filterEdgesByAngle(eF, maxDeg);

        if isToEnd
            % "→ END": find all START poses that can reach anchor via seq.
            % Mirror FROM START: reverse the graph, search from anchor to
            % each candidate start in reversed graph, then flip the paths.
            allResults = findAllSeqPathsToEnd(anchor, seq, numStable, eF, useCost, validPose, maxDeg);
        else
            allResults = findAllSeqPathsFromStart(anchor, seq, numStable, eF, useCost, validPose, maxDeg);
        end

        clearPathHandles();
        resetAllBubbles();
        highlightBubble(anchor, [0.82 0.65 0.08], [1 0.88 0.28]);

        pathCols = seqPathColors(numel(allResults));
        state.seqEdgeLists = {};
        for k = 1:numel(allResults)
            col = pathCols(k,:);
            eiList = allResults(k).edgeIdxList;
            state.seqEdgeLists{end+1} = eiList;
            for j = 1:numel(eiList)
                ei = eiList(j);
                s=eF(ei,1); d=eF(ei,2); r=eF(ei,3); t=eF(ei,4); aDeg=eF(ei,5);
                ls='-'; if t>=10, ls='--'; end
                rad_ = bRad*(0.65+0.10*mod(k-1,3));
                h = drawArrow(axNet,xPos(s),yPos(s),xPos(d),yPos(d), ...
                    rad_,col,2.2,0.95,buildEdgeLabel(r,aDeg),ls);
                pathHandles{end+1}=h;
            end
            if isToEnd
                ep = allResults(k).startPose;
                if ep ~= anchor, highlightBubble(ep, col*0.6+[0.06 0.06 0.08], col); end
            else
                ep = allResults(k).endPose;
                if ep ~= anchor, highlightBubble(ep, col*0.6+[0.06 0.06 0.08], col); end
            end
        end

        updateSeqInfoPanel(anchor, seq, allResults, useCost, isToEnd, false, eF);

        metStr='cost'; if ~useCost, metStr='hops'; end
        dirLabel = sprintf('from Pose %d', anchor);
        if isToEnd, dirLabel = sprintf('to Pose %d', anchor); end
        if isempty(allResults)
            hStatusTxt.String = sprintf('No reachable poses %s with this sequence  |  Click pose to change START', dirLabel);
        else
            hStatusTxt.String = sprintf( ...
                '%d pose(s) reachable %s  |  Metric: %s  |  Click END pose for direct path, click START to reset', ...
                numel(allResults), dirLabel, metStr);
        end
    end

    % ------------------------------------------------------------------
    %  SEQUENCE SEARCH: global (all anchors)
    % ------------------------------------------------------------------
    function runSeqGlobal()
        state.mode = 'seq';
        seq     = getSequence();
        useCost = (hMetricDD.Value == 1);
        isToEnd = (hDirDD.Value == 1);
        allowed = getAllowedTypes();
        maxDeg  = getMaxAngle();
        eF      = filterEdgesByType(edges, allowed);
        eF      = filterEdgesByAngle(eF, maxDeg);

        allResults = struct('startPose',{},'endPose',{},'edgeIdxList',{},'hops',{},'cost',{});
        for anchor = 1:numStable
            if ~validPose(anchor), continue; end
            if isToEnd
                res = findAllSeqPathsToEnd(anchor, seq, numStable, eF, useCost, validPose, maxDeg);
                for ri = 1:numel(res)
                    entry.startPose   = res(ri).startPose;
                    entry.endPose     = anchor;
                    entry.edgeIdxList = res(ri).edgeIdxList;
                    entry.hops        = res(ri).hops;
                    entry.cost        = res(ri).cost;
                    allResults(end+1) = entry; %#ok<AGROW>
                end
            else
                res = findAllSeqPathsFromStart(anchor, seq, numStable, eF, useCost, validPose, maxDeg);
                for ri = 1:numel(res)
                    entry.startPose   = anchor;
                    entry.endPose     = res(ri).endPose;
                    entry.edgeIdxList = res(ri).edgeIdxList;
                    entry.hops        = res(ri).hops;
                    entry.cost        = res(ri).cost;
                    allResults(end+1) = entry; %#ok<AGROW>
                end
            end
        end

        clearPathHandles();
        resetAllBubbles();

        pathCols = seqPathColors(numel(allResults));
        state.seqEdgeLists = {};
        for k = 1:numel(allResults)
            col = pathCols(k,:);
            eiList = allResults(k).edgeIdxList;
            state.seqEdgeLists{end+1} = eiList;
            for j = 1:numel(eiList)
                ei = eiList(j);
                s=eF(ei,1); d=eF(ei,2); r=eF(ei,3); t=eF(ei,4); aDeg=eF(ei,5);
                ls='-'; if t>=10, ls='--'; end
                rad_ = bRad*(0.60+0.08*mod(k-1,4));
                h = drawArrow(axNet,xPos(s),yPos(s),xPos(d),yPos(d), ...
                    rad_,col,1.8,0.85,buildEdgeLabel(r,aDeg),ls);
                pathHandles{end+1}=h;
            end
            sp = allResults(k).startPose; ep = allResults(k).endPose;
            highlightBubble(sp, col*0.5+[0.08 0.06 0.04], [1.0 0.82 0.28]);
            highlightBubble(ep, col*0.5+[0.04 0.08 0.06], col);
        end

        updateSeqInfoPanel(0, seq, allResults, useCost, isToEnd, true, eF);

        metStr='cost'; if ~useCost, metStr='hops'; end
        if isempty(allResults)
            hStatusTxt.String = 'No paths found globally with this sequence';
        else
            hStatusTxt.String = sprintf( ...
                '[Global] %d total path(s) found  |  Metric: %s  |  Adjust slots or scope', ...
                numel(allResults), metStr);
        end
    end

    % ------------------------------------------------------------------
    %  SEQUENCE SEARCH: between specific START and END
    % ------------------------------------------------------------------
    function runSeqBetween(A, B)
        state.mode = 'seq';
        seq     = getSequence();
        useCost = (hMetricDD.Value == 1);
        allowed = getAllowedTypes();
        maxDeg  = getMaxAngle();
        eF      = filterEdgesByType(edges, allowed);
        eF      = filterEdgesByAngle(eF, maxDeg);

        res = seqConstrainedSearch(A, B, seq, numStable, eF, useCost, maxDeg);
        clearPathHandles();

        pathCols = seqPathColors(1);
        col = pathCols(1,:);
        state.seqEdgeLists = {};
        if ~isempty(res.edgeIdxList)
            state.seqEdgeLists{1} = res.edgeIdxList;
            for j = 1:numel(res.edgeIdxList)
                ei = res.edgeIdxList(j);
                s=eF(ei,1); d=eF(ei,2); r=eF(ei,3); t=eF(ei,4); aDeg=eF(ei,5);
                ls='-'; if t>=10, ls='--'; end
                h = drawArrow(axNet,xPos(s),yPos(s),xPos(d),yPos(d), ...
                    bRad,col,2.5,1.0,buildEdgeLabel(r,aDeg),ls);
                pathHandles{end+1}=h;
            end
            hStatusTxt.String = sprintf('Seq path %d→%d: %d hop(s) cost=%.3f', ...
                A,B,res.hops,res.cost);
        else
            hStatusTxt.String = sprintf('No sequence path from %d to %d with this constraint', A, B);
        end

        if ~isempty(res.edgeIdxList)
            ar.startPose   = A;
            ar.endPose     = B;
            ar.edgeIdxList = res.edgeIdxList;
            ar.hops        = res.hops;
            ar.cost        = res.cost;
            updateSeqInfoPanel(B, seq, ar, useCost, true, false, eF);
        else
            clearInfoPanel();
        end
    end

    % ------------------------------------------------------------------
    %  INFO PANEL: path mode
    % ------------------------------------------------------------------
    function updatePathInfoPanel(eF, samePath)
        A=state.selA; B=state.selB;
        lines = {};

        lines{end+1} = sprintf('Pose %d  →  Pose %d', A, B);
        lines{end+1} = repmat('─',1,36);

        if isempty(state.pathEdgesS)
            lines{end+1} = '  No path found.';
        elseif samePath
            costS = sum(eF(state.pathEdgesS,3));
            hopsS = countRealHops(state.pathEdgesS,eF);
            lines{end+1} = sprintf('■ Shortest & Lowest Cost   %d hop(s)   r=%.4f', hopsS, costS);
            merged = mergeConsecutiveSameType(state.pathEdgesS, eF);
            for k2=1:numel(merged)
                lines{end+1} = formatMergedEdge(merged(k2), '  ');
            end
            lines{end+1} = repmat('─',1,36);
        else
            if ~isempty(state.pathEdgesS)
                costS = sum(eF(state.pathEdgesS,3));
                hopsS = countRealHops(state.pathEdgesS,eF);
                lines{end+1} = sprintf('■ Shortest   %d hop(s)   r=%.4f', hopsS, costS);
                merged = mergeConsecutiveSameType(state.pathEdgesS, eF);
                for k2=1:numel(merged)
                    lines{end+1} = formatMergedEdge(merged(k2), '  ');
                end
                lines{end+1} = repmat('─',1,36);
            end
            if ~isempty(state.pathEdgesC)
                costC = sum(eF(state.pathEdgesC,3));
                hopsC = countRealHops(state.pathEdgesC,eF);
                lines{end+1} = sprintf('■ Cheapest   %d hop(s)   r=%.4f', hopsC, costC);
                merged = mergeConsecutiveSameType(state.pathEdgesC, eF);
                for k2=1:numel(merged)
                    lines{end+1} = formatMergedEdge(merged(k2), '  ');
                end
                lines{end+1} = repmat('─',1,36);
            end
        end

        for si=[A,B]
            pos=allPoses(stableIdx(si));
            lbl='START'; if si==B, lbl='END'; end
            lines{end+1} = sprintf('[%s] Pose %d', lbl, si);
            lines{end+1} = sprintf('  Plane=%-2d  θ=%.0f°', pos.floorPlaneIdx, rad2deg(pos.theta));
            lines{end+1} = sprintf('  rW=%.3f  rF=%.3f', ratioWall(si),ratioFloor(si));
            lines{end+1} = sprintf('  Qs=%.4f', Qs(si));
        end

        hInfoList.String = lines;
        hInfoList.Value  = [];
        try jList = findjobj(hInfoList); jList.ensureIndexIsVisible(0); catch; end
    end

    % ------------------------------------------------------------------
    %  INFO PANEL: sequence mode
    % ------------------------------------------------------------------
    function updateSeqInfoPanel(anchor, seq, allResults, useCost, isToEnd, isGlobal, eF)
        lines = {};

        seqStr = seqToString(seq);
        metStr='Lowest Cost'; if ~useCost, metStr='Fewest Hops'; end
        dirLbl='→ END'; if ~isToEnd, dirLbl='FROM START'; end

        if isGlobal
            lines{end+1} = '[Global Sequence Analysis]';
        else
            anchor_lbl='END'; if ~isToEnd, anchor_lbl='START'; end
            lines{end+1} = sprintf('%s: Pose %d', anchor_lbl, anchor);
        end
        lines{end+1} = sprintf('Seq: %s', seqStr);
        lines{end+1} = sprintf('Dir: %s     Met: %s', dirLbl, metStr);
        lines{end+1} = repmat('─',1,36);

        if isempty(allResults)
            lines{end+1} = '  No paths found.';
            hInfoList.String = lines;
            hInfoList.Value = [];
            return;
        end

        nr = numel(allResults);
        for k=1:nr
            if isstruct(allResults) && ~iscell(allResults)
                res = allResults(k);
            else
                res = allResults{k};
            end

            sp = 0; ep = 0;
            if isfield(res,'startPose'), sp = res.startPose; end
            if isfield(res,'endPose'),   ep = res.endPose;   end

            if isGlobal
                lines{end+1} = sprintf('■ Path %d:  %d → %d', k, sp, ep);
            elseif isToEnd && sp > 0
                lines{end+1} = sprintf('■ START %d → END %d', sp, anchor);
            elseif ~isToEnd && ep > 0
                lines{end+1} = sprintf('■ START %d → END %d', anchor, ep);
            else
                lines{end+1} = sprintf('■ Path %d', k);
            end
            lhops = countRealHops(res.edgeIdxList, eF);
            lines{end+1} = sprintf('  %d logical hop(s)   cost=%.4f', lhops, res.cost);
            merged = mergeConsecutiveSameType(res.edgeIdxList, eF);
            for j=1:numel(merged)
                lines{end+1} = formatMergedEdge(merged(j), '  ');
            end
            if k < nr
                lines{end+1} = repmat('─',1,36);
            end
        end

        hInfoList.String = lines;
        hInfoList.Value  = [];
        try jList = findjobj(hInfoList); jList.ensureIndexIsVisible(0); catch; end
    end

    % ------------------------------------------------------------------
    %  Helpers
    % ------------------------------------------------------------------
    function active = isSeqActive()
        active = false;
        for k2 = 1:numSeqSlots
            if hSeqDD(k2).Value > 1, active=true; return; end
        end
    end

    function seq = getSequence()
        % Returns a vector of transition CATEGORIES (1=Wall, 2=Floor) —
        % NOT primTypes. The dropdown itself only offers Wall/Floor now,
        % so the popup Value (2=Wall,3=Floor) maps directly: category = v-1.
        seq = [];
        for k2 = 1:numSeqSlots
            v = hSeqDD(k2).Value;
            if v > 1, seq(end+1) = v-1; end %#ok<AGROW>
        end
    end

    function clearPathHandles()
        for k2=1:numel(pathHandles)
            h=pathHandles{k2};
            if iscell(h), for j=1:numel(h), try delete(h{j}); catch; end; end
            else, try delete(h); catch; end
            end
        end
        pathHandles={};
    end

    function highlightBubble(si,fc,ec)
        if si < 1 || si > numStable, return; end
        hBubble(si).FaceColor=fc; hBubble(si).EdgeColor=ec; hBubble(si).LineWidth=3.5;
    end

    function resetAllBubbles()
        for si=1:numStable
            [fc,ec]=bubbleColors(si,ratioWall,ratioFloor, ...
                THRESH_WALL,THRESH_FLOOR,Qs);
            hBubble(si).FaceColor=fc; hBubble(si).EdgeColor=ec; hBubble(si).LineWidth=2.0;
        end
    end

    function clearInfoPanel()
        hInfoList.String = {'  — select poses or run sequence —'};
        hInfoList.Value  = [];
    end

end % launchUI


% =========================================================================
%  formatMergedEdge
% =========================================================================
function s = formatMergedEdge(me, indent)
if nargin < 2, indent = '  '; end
typeStr  = me.typeStr;
hopStr   = '';
if me.hopCount > 1
    hopStr = sprintf('×%d', me.hopCount);
end
costStr  = sprintf('r=%.4f', me.totalCost);
angStr   = '';
if me.totalAngle > 0.5
    angStr = sprintf('  Δ%.0f°', me.totalAngle);
end
s = sprintf('%s%d→%d  [%s%s  %s%s]', indent, me.srcNode, me.dstNode, ...
    typeStr, hopStr, costStr, angStr);
end


% =========================================================================
%  mergeConsecutiveSameType
%  Groups consecutive REAL (cost>0) edges into a single displayed "logical
%  hop" whenever they share the same transition CATEGORY (Wall or Floor),
%  regardless of Up/Down sub-type. E.g. Up-Wall immediately followed by
%  Down-Wall merges into one "Wall ×2" entry, matching the sequence-search
%  semantics in seqConstrainedSearch.
% =========================================================================
function mergedEdges = mergeConsecutiveSameType(edgeIdxList, edges)
mergedEdges = struct('srcNode',{},'dstNode',{},'totalCost',{},'totalAngle',{}, ...
    'typeStr',{},'hopCount',{},'isMerged',{});
if isempty(edgeIdxList), return; end
k=1;
while k<=numel(edgeIdxList)
    ei=edgeIdxList(k); s=edges(ei,1); d=edges(ei,2); r=edges(ei,3);
    t=edges(ei,4); aDeg=edges(ei,5);
    primType=mod(t-1,10)+1; cat0=catOf(primType); isReal=(r>0);
    if isReal
        totalCost=r; totalAngle=aDeg; hopCount=1; lastDst=d; pureType=primType;
        while k+1<=numel(edgeIdxList)
            ei2=edgeIdxList(k+1); t2=edges(ei2,4); r2=edges(ei2,3);
            prim2=mod(t2-1,10)+1;
            if catOf(prim2)==cat0 && r2>0
                totalCost=totalCost+r2; totalAngle=totalAngle+edges(ei2,5);
                hopCount=hopCount+1; lastDst=edges(ei2,2); k=k+1;
                if prim2~=pureType, pureType=-1; end   % mixed Up/Down within the run
            else; break; end
        end
        me.srcNode=s; me.dstNode=lastDst; me.totalCost=totalCost;
        me.totalAngle=totalAngle;
        if hopCount>1 && pureType==-1
            me.typeStr = catLabel(primType);           % mixed run -> generic "Wall"/"Floor"
        else
            me.typeStr = edgeTypeStr(t);                % single/uniform -> exact sub-type
        end
        me.hopCount=hopCount; me.isMerged=(hopCount>1);
    else
        me.srcNode=s; me.dstNode=d; me.totalCost=0; me.totalAngle=aDeg;
        me.typeStr=edgeTypeStr(t); me.hopCount=1; me.isMerged=false;
    end
    mergedEdges(end+1)=me; k=k+1; %#ok<AGROW>
end
end


% =========================================================================
%  countRealHops
%  Counts logical hops by transition CATEGORY (Wall/Floor), so a run that
%  alternates Up-Wall/Down-Wall still counts as a single hop — consistent
%  with mergeConsecutiveSameType and the sequence-search matching logic.
% =========================================================================
function n = countRealHops(edgeIdxList, edges)
n        = 0;
lastCat  = -1;
for k = 1:numel(edgeIdxList)
    ei       = edgeIdxList(k);
    cost     = edges(ei, 3);
    if cost == 0, continue; end
    primType = mod(edges(ei, 4) - 1, 10) + 1;
    cat0     = catOf(primType);
    if cat0 ~= lastCat
        n       = n + 1;
        lastCat = cat0;
    end
end
end

function c = computePathCost(edgeIdxList, edges)
c = 0;
for b = 1:numel(edgeIdxList)
    c = c + edges(edgeIdxList(b), 3);
end
end


% =========================================================================
%  computeOptimalSequencesGlobal
% =========================================================================
function results = computeOptimalSequencesGlobal( ...
        numStable, edges, validPose, useCost, maxDeg)

results = struct('seq',{}, 'direction',{}, 'coverage',{}, ...
                 'totalHops',{}, 'avgHops',{}, 'avgMetric',{}, 'totalMetric',{});

allSeqs = enumSeqs(4);

for si = 1:numel(allSeqs)
    seq = allSeqs{si};

    for dirIdx = 1:2
        isToEnd     = (dirIdx == 1);
        coverage    = 0;
        totalHops   = 0;
        totalMetric = 0;

        for anchor = 1:numStable
            if ~validPose(anchor), continue; end

            if isToEnd
                res = findAllSeqPathsToEnd(anchor, seq, numStable, edges, useCost, validPose, maxDeg);
            else
                res = findAllSeqPathsFromStart(anchor, seq, numStable, edges, useCost, validPose, maxDeg);
            end

            for ri = 1:numel(res)
                lh          = countRealHops(res(ri).edgeIdxList, edges);
                coverage    = coverage    + 1;
                totalHops   = totalHops   + lh;
                if useCost
                    totalMetric = totalMetric + ...
                        computePathCost(res(ri).edgeIdxList, edges);
                else
                    totalMetric = totalMetric + lh;
                end
            end
        end

        if coverage == 0, continue; end

        entry.seq         = seq;
        entry.direction   = dirStr(isToEnd);
        entry.coverage    = coverage;
        entry.totalHops   = totalHops;
        entry.avgHops     = totalHops   / coverage;
        entry.avgMetric   = totalMetric / coverage;
        entry.totalMetric = totalMetric;
        results(end+1) = entry; %#ok<AGROW>
    end
end

if ~isempty(results)
    covs    = [results.coverage]';
    seqLens = cellfun(@numel, {results.seq})';
    avgs    = [results.avgMetric]';
    [~, ord] = sortrows([-covs, seqLens, avgs]);
    results  = results(ord);
    results  = results(1:min(20, numel(results)));
end

    function s = dirStr(toEnd)
        if toEnd, s = '→ END'; else, s = 'FROM START'; end
    end
end


% =========================================================================
%  computeOptimalSequencesIndividual
% =========================================================================
function results = computeOptimalSequencesIndividual( ...
        anchor, numStable, edges, validPose, useCost, isToEnd, maxDeg)

results = struct('seq',{}, 'direction',{}, 'coverage',{}, ...
                 'totalHops',{}, 'avgHops',{}, 'avgMetric',{}, 'totalMetric',{});

if ~validPose(anchor), return; end

allSeqs = enumSeqs(4);

for si = 1:numel(allSeqs)
    seq = allSeqs{si};

    if isToEnd
        res = findAllSeqPathsToEnd(anchor, seq, numStable, edges, useCost, validPose, maxDeg);
    else
        res = findAllSeqPathsFromStart(anchor, seq, numStable, edges, useCost, validPose, maxDeg);
    end
    if isempty(res), continue; end

    coverage    = numel(res);
    totalHops   = 0;
    totalMetric = 0;

    for ri = 1:coverage
        lh          = countRealHops(res(ri).edgeIdxList, edges);
        totalHops   = totalHops   + lh;
        if useCost
            totalMetric = totalMetric + computePathCost(res(ri).edgeIdxList, edges);
        else
            totalMetric = totalMetric + lh;
        end
    end

    entry.seq         = seq;
    entry.direction   = dirStr(isToEnd);
    entry.coverage    = coverage;
    entry.totalHops   = totalHops;
    entry.avgHops     = totalHops   / coverage;
    entry.avgMetric   = totalMetric / coverage;
    entry.totalMetric = totalMetric;
    results(end+1) = entry; %#ok<AGROW>
end

if ~isempty(results)
    covs  = [results.coverage];
    thops = [results.totalHops];
    avgs  = [results.avgMetric];
    [~, ord] = sortrows([-covs(:), thops(:), avgs(:)]);
    results  = results(ord);
    results  = results(1:min(20, numel(results)));
end

    function s = dirStr(toEnd)
        if toEnd, s = '→ END'; else, s = 'FROM START'; end
    end
end


% =========================================================================
%  enumSeqs
%  Enumerates sequences over the transition CATEGORY alphabet {1=Wall,
%  2=Floor} rather than the 4 primTypes, since the sequence UI/search
%  layer now only distinguishes Wall vs Floor.
% =========================================================================
function seqs = enumSeqs(maxLen)
seqs = {};
for len = 1:maxLen
    buildSeq([], len);
end

    function buildSeq(current, remaining)
        if remaining == 0
            seqs{end+1} = current; %#ok<AGROW>
            return;
        end
        for t = 1:2
            if isempty(current) || current(end) ~= t
                buildSeq([current, t], remaining - 1);
            end
        end
    end
end


% =========================================================================
%  renderPoseAxes
% =========================================================================
function renderPoseAxes(ax, pos, si, label, Rchute, partFaces)
cla(ax); hold(ax,'on'); axis(ax,'equal'); grid(ax,'on');
set(ax,'Color',[0.09 0.10 0.13],'XColor',[0.38 0.40 0.48],'YColor',[0.38 0.40 0.48], ...
    'ZColor',[0.38 0.40 0.48],'GridColor',[0.22 0.23 0.28],'GridAlpha',0.3,'FontSize',7);
vertsC=pos.verticesWorld; centC=pos.centroidWorld(:)';
floorIdx=pos.floorContactVertIdx; wallIdx=pos.wallContactVertIdx;
floorZ=min(vertsC(floorIdx,3)); vertsC(:,3)=vertsC(:,3)-floorZ; centC(3)=centC(3)-floorZ;
wallY=max(vertsC(wallIdx,2)); vertsC(:,2)=vertsC(:,2)-wallY; centC(2)=centC(2)-wallY;
vertsW=(Rchute*vertsC')'; centW=(Rchute*centC(:))';

% ── DISPLAY-ONLY mirror + rotate ────────────────────────────────────────
% Flip the rendered geometry about the yz-plane (x -> -x). This is purely
% a viewer transform: it only affects vertsW/centW/partFaces as used in
% THIS function, never the underlying pos data, so nothing downstream
% (transitions, sequence search, ratios, etc.) is touched.
% Mirroring reverses triangle winding (and therefore face-normal
% direction), so the face vertex order is flipped to keep lighting
% correct. The camera azimuth is negated ("...and rotate") to compensate
% for the handedness flip so the part still reads naturally in 3D.
vertsW(:,1) = -vertsW(:,1);
centW(1)    = -centW(1);
if ~isempty(partFaces)
    partFaces = partFaces(:, [1 3 2]);
end
view(ax, -40, 24);

if ~isempty(partFaces)
    trisurf(partFaces,vertsW(:,1),vertsW(:,2),vertsW(:,3),'Parent',ax, ...
        'FaceColor',[1.00 0.92 0.22],'FaceAlpha',0.55,'EdgeColor',[0.38 0.32 0.00],'EdgeAlpha',0.12);
else
    scatter3(ax,vertsW(:,1),vertsW(:,2),vertsW(:,3),8,[0.68 0.68 0.68],'filled','MarkerFaceAlpha',0.5);
end
scatter3(ax,vertsW(floorIdx,1),vertsW(floorIdx,2),vertsW(floorIdx,3),20,[0.18 0.82 0.38],'filled');
scatter3(ax,vertsW(wallIdx,1),vertsW(wallIdx,2),vertsW(wallIdx,3),20,[0.88 0.42 0.18],'filled');
scatter3(ax,centW(1),centW(2),centW(3),50,'w','filled');
lighting(ax,'gouraud'); camlight(ax,'headlight');
if strcmp(label,'START'), titleCol=[1.00 0.82 0.14]; else, titleCol=[0.18 0.88 0.52]; end
title(ax,sprintf('[%s] Pose %d  P%d  θ=%.0f°',label,si,pos.floorPlaneIdx,rad2deg(pos.theta)), ...
    'Color',titleCol,'FontSize',8,'FontWeight','bold');
end


% =========================================================================
%  buildEdgeLabel
% =========================================================================
function lbl = buildEdgeLabel(r, aDeg)
if r==0 && aDeg<=0.5,      lbl='·';
elseif r==0,               lbl=sprintf('Δ%.0f°',aDeg);
elseif aDeg>0.5,           lbl=sprintf('%.2f\nΔ%.0f°',r,aDeg);
else,                      lbl=sprintf('%.2f',r);
end
end


% =========================================================================
%  SEQUENCE PATH SEARCH  — TO END (fixed: mirrors FROM START logic)
%  Reverses the graph, searches from endPose to each candidate startPose,
%  then flips paths back to forward direction.
% =========================================================================
function allResults = findAllSeqPathsToEnd(endPose, seq, numNodes, edges, useCost, validPose, maxDeg)
allResults = struct('startPose',{},'edgeIdxList',{},'hops',{},'cost',{});
% Build reversed edge table
edgesRev = [edges(:,2), edges(:,1), edges(:,3), edges(:,4), edges(:,5)];
for startPose = 1:numNodes
    if startPose == endPose, continue; end
    if ~validPose(startPose), continue; end
    % In reversed graph: search from endPose (now "start") to startPose (now "goal")
    res = seqConstrainedSearch(endPose, startPose, seq, numNodes, edgesRev, useCost, maxDeg);
    if ~isempty(res.edgeIdxList)
        % Map reversed edge indices back to forward edge indices and flip path
        fwdEdgeList = mapReversedEdges(res.edgeIdxList, edges, edgesRev);
        entry.startPose   = startPose;
        entry.edgeIdxList = fwdEdgeList;
        entry.hops        = countRealHops(fwdEdgeList, edges);
        entry.cost        = computePathCost(fwdEdgeList, edges);
        allResults(end+1) = entry; %#ok<AGROW>
    end
end
end


% =========================================================================
%  SEQUENCE PATH SEARCH  — from START
% =========================================================================
function allResults = findAllSeqPathsFromStart(startPose, seq, numNodes, edges, useCost, validPose, maxDeg)
allResults = struct('endPose',{},'edgeIdxList',{},'hops',{},'cost',{});
edgesRev = [edges(:,2), edges(:,1), edges(:,3), edges(:,4), edges(:,5)];
for endPose = 1:numNodes
    if endPose == startPose, continue; end
    if ~validPose(endPose), continue; end
    res = seqConstrainedSearch(endPose, startPose, seq, numNodes, edgesRev, useCost, maxDeg);
    if ~isempty(res.edgeIdxList)
        fwdEdgeList = mapReversedEdges(res.edgeIdxList, edges, edgesRev);
        entry.endPose     = endPose;
        entry.edgeIdxList = fwdEdgeList;
        entry.hops        = countRealHops(fwdEdgeList, edges);
        entry.cost        = computePathCost(fwdEdgeList, edges);
        allResults(end+1) = entry; %#ok<AGROW>
    end
end
end

function fwdList = mapReversedEdges(revIdxList, edges, edgesRev)
fwdList = zeros(1, numel(revIdxList));
for k = 1:numel(revIdxList)
    ri=revIdxList(k); rs=edgesRev(ri,1); rd=edgesRev(ri,2); rt=edgesRev(ri,4);
    matches = find(edges(:,1)==rd & edges(:,2)==rs & edges(:,4)==rt);
    if ~isempty(matches), fwdList(k)=matches(1); else, fwdList(k)=ri; end
end
fwdList = fliplr(fwdList);
end


% =========================================================================
%  seqConstrainedSearch
%  seq is a vector of CATEGORIES (1=Wall, 2=Floor). A slot is satisfied
%  by ANY edge whose catOf(primType) matches — i.e. Up-Wall and Down-Wall
%  are interchangeable within a Wall slot, and likewise for Floor. This
%  applies both when first entering a slot and when extending a run
%  within the same slot (so "up-wall, down-wall, down-floor" satisfies
%  the same [Wall, Floor] sequence as "down-wall, down-wall, down-floor").
% =========================================================================
function result = seqConstrainedSearch(startPose, endPose, seq, numNodes, edges, useCost, maxDeg)
result.startPose   = startPose;
result.endPose     = endPose;
result.edgeIdxList = [];
result.hops        = 0;
result.cost        = 0;

nSeq = numel(seq);
INF  = 1e10;

numStates = numNodes * (nSeq + 1);
dist        = INF * ones(numStates, 1);
prevState   = zeros(numStates, 1);
prevEdge    = zeros(numStates, 1);
runAngle    = zeros(numStates, 1);   % cumulative angle of current same-category run
runStart    = zeros(numStates, 1);   % node where the current run began
unvisited   = true(numStates, 1);

s0          = 0 * numNodes + startPose;
dist(s0)    = 0;
runAngle(s0)= 0;
runStart(s0)= startPose;

for iter = 1:numStates
    tmp = dist; tmp(~unvisited) = INF;
    [dMin, curState] = min(tmp);
    if dMin >= INF, break; end
    unvisited(curState) = false;

    curNode = mod(curState - 1, numNodes) + 1;
    curSP   = floor((curState - 1) / numNodes);

    if curNode == endPose && curNode ~= startPose && (nSeq == 0 || curSP == nSeq)
        eiList = reconstructPath(curState, prevState, prevEdge);
        if ~isempty(eiList) && countRealHops(eiList, edges) > 0
            result.edgeIdxList = eiList;
            result.cost        = computePathCost(eiList, edges);
            result.hops        = countRealHops(eiList, edges);
        end
        return;
    end

    eis = find(edges(:,1) == curNode);

    % --- FREE SLOT SKIP: satisfy the next sequence slot without moving ---
    % Lets the search consume a required slot "for free" when the pose
    % already there doesn't need a physical transition to satisfy it.
    % Zero cost / zero angle, so Dijkstra only ever uses it when it is
    % genuinely at least as good as a real hop — it can't crowd out an
    % actually-needed transition, but it will always beat a pointless
    % loop-back that only exists to burn a slot.
    if nSeq > 0 && curSP < nSeq
        ns  = (curSP + 1) * numNodes + curNode;
        alt = dist(curState);
        if alt < dist(ns)
            dist(ns)      = alt;
            prevState(ns) = curState;
            prevEdge(ns)  = -1;      % sentinel: virtual skip, no real edge
            runAngle(ns)  = 0;
            runStart(ns)  = curNode;
        end
    end

    for k = 1:numel(eis)
        ei       = eis(k);
        nb       = edges(ei, 2);
        eCost    = edges(ei, 3);
        eType    = edges(ei, 4);
        primType = mod(eType - 1, 10) + 1;
        eCat     = catOf(primType);   % Wall/Floor category — direction-agnostic
        eAng     = edges(ei, 5);
        isChain  = (eCost == 0);

        if isChain
            % free continuation of the SAME logical transition — doesn't
            % start a new run and carries no extra angle
            ns  = curSP * numNodes + nb;
            alt = dist(curState);
            if alt < dist(ns)
                dist(ns)      = alt;
                prevState(ns) = curState;
                prevEdge(ns)  = ei;
                runAngle(ns)  = runAngle(curState);
                runStart(ns)  = runStart(curState);
            end

        elseif nSeq == 0
            ns  = nb;
            w   = eCost; if ~useCost, w = 1; end
            alt = dist(curState) + w;
            if alt < dist(ns)
                dist(ns)      = alt;
                prevState(ns) = curState;
                prevEdge(ns)  = ei;
                runAngle(ns)  = eAng;
                runStart(ns)  = curNode;
            end

        else
            % --- advance to a NEW sequence slot: starts a fresh run ---
            % Any Wall edge (Up or Down) satisfies a "Wall" slot; any
            % Floor edge (Up or Down) satisfies a "Floor" slot.
            if curSP < nSeq && eCat == seq(curSP + 1)
                newSP = curSP + 1;
                ns    = newSP * numNodes + nb;
                w     = eCost; if ~useCost, w = 1; end
                alt   = dist(curState) + w;
                if alt < dist(ns)
                    dist(ns)      = alt;
                    prevState(ns) = curState;
                    prevEdge(ns)  = ei;
                    runAngle(ns)  = eAng;      % fresh run, resets to this hop's own angle
                    runStart(ns)  = curNode;   % this hop's own start
                end
            end

            % --- stay in the SAME slot: extend the current same-category run ---
            % An Up-Wall hop can directly follow a Down-Wall hop (or vice
            % versa) within the same "Wall" slot — only the category has
            % to match, not the exact sub-type.
            if curSP > 0 && eCat == seq(curSP)
                newRunAngle = runAngle(curState) + eAng;
                if newRunAngle > maxDeg
                    continue;   % cumulative rotation over the cap — prune
                end
                ns  = curSP * numNodes + nb;
                w   = eCost; if ~useCost, w = 1; end
                alt = dist(curState) + w;
                if alt < dist(ns)
                    dist(ns)      = alt;
                    prevState(ns) = curState;
                    prevEdge(ns)  = ei;
                    runAngle(ns)  = newRunAngle;
                    runStart(ns)  = runStart(curState);
                end
            end
        end
    end
end

end

function edgeIdxList = reconstructPath(goalState, prevState, prevEdge)
edgeIdxList = []; cur = goalState;
while true
    ei = prevEdge(cur);
    if ei == 0
        break;                              % true path start — no predecessor
    elseif ei == -1
        cur = prevState(cur);               % virtual slot skip — no edge to record
    else
        edgeIdxList = [ei, edgeIdxList]; %#ok<AGROW>
        cur = prevState(cur);
    end
end
end


% =========================================================================
%  seqToString
%  seq entries are CATEGORIES now (1=Wall, 2=Floor).
% =========================================================================
function s = seqToString(seq)
names = {'Wall','Floor'};
if isempty(seq), s = '(none)'; return; end
parts = cell(1, numel(seq));
for k = 1:numel(seq), parts{k} = names{seq(k)}; end
s = strjoin(parts, ' → ');
end


% =========================================================================
%  seqPathColors
% =========================================================================
function cols = seqPathColors(n)
base = [1.00 0.75 0.15; 0.18 0.88 0.85; 0.90 0.32 0.42; 0.42 0.85 0.32; ...
        0.68 0.42 0.95; 0.95 0.58 0.18; 0.32 0.62 0.95; 0.85 0.85 0.28];
cols = zeros(n, 3);
for k = 1:n, cols(k,:) = base(mod(k-1, size(base,1))+1, :); end
end


% =========================================================================
%  bubbleColors
% =========================================================================
function [faceC, edgeC] = bubbleColors(si, ratioWall, ratioFloor, ...
    THRESH_WALL, THRESH_FLOOR, Qs)
hasWall    = (ratioWall(si)  >= THRESH_WALL);
hasFloor   = (ratioFloor(si) >= THRESH_FLOOR);
isZeroed   = (Qs(si) == 0);
if isZeroed,              faceC=[0.38 0.11 0.11]; edgeC=[0.72 0.22 0.22];
elseif hasWall&&hasFloor, faceC=[0.48 0.22 0.62]; edgeC=[0.78 0.52 0.92];
elseif hasWall,            faceC=[0.12 0.32 0.60]; edgeC=[0.38 0.62 0.92];
elseif hasFloor,           faceC=[0.10 0.42 0.28]; edgeC=[0.28 0.78 0.52];
else,                      faceC=[0.20 0.22 0.28]; edgeC=[0.52 0.56 0.68];
end
end


% =========================================================================
%  edgeTypeStr
% =========================================================================
function s = edgeTypeStr(t)
switch t
    case 1,  s='dnWall';   case 2,  s='dnFloor';  case 3,  s='upWall';   case 4,  s='upFloor';
    case 10, s='dnWall↓'; case 20, s='dnFloor↓'; case 30, s='upWall↓'; case 40, s='upFloor↓';
    otherwise, s='?';
end
end


% =========================================================================
%  drawArrow
% =========================================================================
function handles = drawArrow(ax, x1,y1, x2,y2, bRad, col, lw, alpha, labelStr, ls)
handles = {}; if nargin < 11, ls = '-'; end
dx=x2-x1; dy=y2-y1; d=sqrt(dx^2+dy^2); if d<1e-6, return; end
ux=dx/d; uy=dy/d; perp=0.12; ox=-uy*perp; oy=ux*perp;
sx=x1+ux*bRad+ox; sy=y1+uy*bRad+oy;
ex=x2-ux*bRad+ox; ey=y2-uy*bRad+oy;
mx=(sx+ex)/2+ox*1.5; my=(sy+ey)/2+oy*1.5;
t=linspace(0,1,40);
bx=(1-t).^2.*sx+2*(1-t).*t.*mx+t.^2.*ex;
by=(1-t).^2.*sy+2*(1-t).*t.*my+t.^2.*ey;
h1=plot(ax,bx,by,ls,'Color',col,'LineWidth',lw);
try h1.Color=[col(:)',alpha]; catch; h1.Color=col*alpha+[0.11 0.12 0.15]*(1-alpha); end
handles{end+1}=h1;
aLen=0.18; aAng=25; ex2=bx(end); ey2=by(end); ex1=bx(end-3); ey1=by(end-3);
ang=atan2(ey2-ey1,ex2-ex1); ang1=ang+deg2rad(180+aAng); ang2=ang+deg2rad(180-aAng);
axH=[ex2,ex2+aLen*cos(ang1),ex2+aLen*cos(ang2),ex2];
ayH=[ey2,ey2+aLen*sin(ang1),ey2+aLen*sin(ang2),ey2];
h2=fill(ax,axH,ayH,col,'EdgeColor',col,'FaceAlpha',alpha);
handles{end+1}=h2;
if ~isempty(labelStr)
    lx=bx(20)+ox*0.5; ly=by(20)+oy*0.5;
    h3=text(ax,lx,ly,labelStr,'Color',col,'FontSize',7,'FontWeight','bold', ...
        'HorizontalAlignment','center','BackgroundColor',[0.09 0.10 0.13],'Margin',1);
    handles{end+1}=h3;
end
end


% =========================================================================
%  PATHFINDING (normal mode)
% =========================================================================
function [path, edgeIdxPath] = bfsShortestPath(src, dst, numNodes, edges)
path=[]; edgeIdxPath=[];
if src==dst, path=src; return; end
INF=1e10; dist=INF*ones(numNodes,1); prev=zeros(numNodes,1); prevEdge=zeros(numNodes,1);
dist(src)=0; unvisited=true(numNodes,1);
for iter=1:numNodes
    tmp=dist; tmp(~unvisited)=INF; [dMin,cur]=min(tmp);
    if dMin>=INF, break; end; if cur==dst, break; end
    unvisited(cur)=false; mask=edges(:,1)==cur; eis=find(mask);
    for k=1:numel(eis)
        ei=eis(k); nb=edges(ei,2); cost=edges(ei,3);
        if ~unvisited(nb), continue; end
        hopW=double(cost>0); alt=dist(cur)+hopW;
        if alt<dist(nb), dist(nb)=alt; prev(nb)=cur; prevEdge(nb)=ei; end
    end
end
if prev(dst)==0, return; end
node=dst; path=dst; edgeIdxPath=[];
while prev(node)~=0
    edgeIdxPath=[prevEdge(node),edgeIdxPath]; node=prev(node); path=[node,path]; %#ok<AGROW>
end
end

function [path, edgeIdxPath] = dijkstraCheapestPath(src, dst, numNodes, edges)
path=[]; edgeIdxPath=[];
if src==dst, path=src; return; end
INF=1e10; dist=INF*ones(numNodes,1); prev=zeros(numNodes,1); prevEdge=zeros(numNodes,1);
dist(src)=0; unvisited=true(numNodes,1);
for iter=1:numNodes
    tmp=dist; tmp(~unvisited)=INF; [dMin,cur]=min(tmp);
    if dMin>=INF, break; end; if cur==dst, break; end
    unvisited(cur)=false; mask=edges(:,1)==cur; eis=find(mask);
    for k=1:numel(eis)
        ei=eis(k); nb=edges(ei,2); cost=edges(ei,3);
        if ~unvisited(nb), continue; end
        alt=dist(cur)+cost;
        if alt<dist(nb), dist(nb)=alt; prev(nb)=cur; prevEdge(nb)=ei; end
    end
end
if prev(dst)==0, return; end
node=dst; path=dst; edgeIdxPath=[];
while prev(node)~=0
    edgeIdxPath=[prevEdge(node),edgeIdxPath]; node=prev(node); path=[node,path]; %#ok<AGROW>
end
end


% =========================================================================
%  resolveTransitionChain
% =========================================================================
function [chainNodes,firstCost,transType,angleChangeDeg]=resolveTransitionChain( ...
    srcSi,transType,allPoses,stableIdx,refQuats,refQuatsAll, ...
    Qs,validPose, ...
    quatMatchTol,planeTol,wallTol,maxChain, ...
    restingPlaneVerts,restingPlaneEqs,chullVertexIdx, ...
    centroidCoords,planeQuats,angleCapDeg,angleCostDenom)

chainNodes=[]; firstCost=0; angleChangeDeg=0;

% NOTE: cost is no longer seeded from ratioWall/ratioFloor here. All four
% transition types (Down-Wall, Down-Floor, Up-Wall, Up-Floor) now use a
% single, uniform cost basis: cumulative rotation angle, scaled by
% angleCostDenom and capped at 1 (see firstCost assignment below, once a
% valid landing pose is found). ratioWall/ratioFloor are still used
% elsewhere (Qs zeroing, bubbleColors) — just not for edge cost.
%
% angleCapDeg is now the SAME 180° single-hop budget for all four
% transition types (see MAX_HOP_ANGLE_DEG at the top of the script). If
% the cumulative rotation across this chain's physical steps would exceed
% it before reaching a valid landing pose, the chain is abandoned (no
% edge built) — e.g. a 180° flip followed by a further 90° rotation is
% two transitions, not one. The main graph-building loop already calls
% this function fresh for every node (valid or not) as a source, so the
% second transition is available as its own edge once the first lands.
%
% NOTE: this function still builds the RAW per-primType (1-4) physical
% edges — the Wall/Floor category collapsing happens one layer up, in the
% sequence-search functions (seqConstrainedSearch, countRealHops,
% mergeConsecutiveSameType) via catOf(). This function is unchanged.

visited=false(numel(stableIdx),1); visited(srcSi)=true;
curSi=srcSi; chain=srcSi; cumAngle=0;

for chainStep=1:maxChain
    [vertsC,centC]=initPose(curSi,allPoses,stableIdx); pos0=allPoses(stableIdx(curSi));

    switch transType
        case 1
            contactIdx=pos0.wallContactVertIdx(:);
            [vertsR,centR,~,~,~,q_rot]=wallRotation_q(vertsC,centC,contactIdx);
        case 2
            contactIdx=pos0.floorContactVertIdx(:);
            [vertsR,centR,~,~,~,q_rot]=floorRotation_q(vertsC,centC,contactIdx);
        case 3
            contactIdx=pos0.wallContactVertIdx(:);
            [vertsR,centR,~,~,~,q_rot]=wallRotationUp_q(vertsC,centC,contactIdx);
        case 4
            contactIdx=pos0.floorContactVertIdx(:);
            [vertsR,centR,~,~,~,q_rot]=floorRotationUp_q(vertsC,centC,contactIdx);
    end

    if isempty(q_rot), chainNodes=[]; return; end

    q_id=[1 0 0 0];
    stepAngleDeg=rad2deg(q_geodesic(q_id,q_rot));
    cumAngle=cumAngle+stepAngleDeg;
    if chainStep==1, angleChangeDeg=stepAngleDeg; end

    if cumAngle > angleCapDeg
        chainNodes=[]; return;
    end

    q_acc=refQuats(curSi,:);
    q_acc=q_compose(q_acc,q_rot);
    [~,mi]=max(abs(q_acc)); if q_acc(mi)<0, q_acc=-q_acc; end

    if isempty(vertsR), [vertsR,centR,~,~]=reseat(vertsC,centC); end
    [matchSi,~]=matchQuatComposed(q_acc,refQuatsAll,quatMatchTol);
    if matchSi==0
        matchSi=secondarySettle(vertsR,centR, ...
            restingPlaneVerts,restingPlaneEqs,chullVertexIdx, ...
            centroidCoords,planeQuats,refQuatsAll, ...
            planeTol,wallTol,quatMatchTol*3);
    end
    if matchSi==0||matchSi==srcSi||visited(matchSi), chainNodes=[]; return; end

    chain(end+1)=matchSi; %#ok<AGROW>
    if validPose(matchSi)
        chainNodes=chain;
        % Uniform angle-based cost for ALL transition types, capped at 1:
        %   1 deg of cumulative rotation = 1/angleCostDenom cost.
        firstCost=min(cumAngle/angleCostDenom, 1);
        return;
    end
    visited(matchSi)=true; curSi=matchSi;
end
chainNodes=[];
end


% =========================================================================
%  secondarySettle
% =========================================================================
function matchSi=secondarySettle(vertsIn,centIn, ...
    restingPlaneVerts,restingPlaneEqs,chullVertexIdx, ...
    centroidCoords,planeQuats,refQuatsAll, ...
    planeTol,wallTol,matchTol)
matchSi=0; [vertsS,centS,~,~]=reseat(vertsIn,centIn);
zMin=min(vertsS(:,3)); onFloor=find(abs(vertsS(:,3)-zMin)<planeTol);
if numel(onFloor)<2, return; end
bestPlane=0; bestOverlap=0;
for pi_=1:numel(restingPlaneVerts)
    ov=numel(intersect(onFloor,restingPlaneVerts{pi_}));
    if ov>bestOverlap, bestOverlap=ov; bestPlane=pi_; end
end
if bestPlane==0, bestPlane=1; end
q_align=planeQuats{bestPlane}; R_floor=q_toRotm(q_align); C0=centS(:)';
vertsRot=(R_floor*(vertsS'-C0'))'+C0;
thetaSamples=linspace(0,2*pi,361); thetaSamples=thetaSamples(1:end-1);
bestDist=inf;
for ti=1:numel(thetaSamples)
    theta=thetaSamples(ti); allVW=rotatePtsAroundZ(vertsRot,C0,theta);
    fZc=min(allVW(:,3)); onF=find(abs(allVW(:,3)-fZc)<planeTol);
    if numel(onF)<2, continue; end
    projXY=allVW(chullVertexIdx,1:2); yWall=max(projXY(:,2));
    onWL=find(abs(projXY(:,2)-yWall)<wallTol);
    if numel(onWL)<2, continue; end
    wX=projXY(onWL,1); cX=C0(1);
    if cX<min(wX)-wallTol||cX>max(wX)+wallTol, continue; end
    q_theta=q_fromAxisAngle([0,0,1],theta); q_cand=q_compose(q_align,q_theta);
    [~,mi]=max(abs(q_cand)); if q_cand(mi)<0, q_cand=-q_cand; end
    [mSi,dist]=matchQuatComposed(q_cand,refQuatsAll,matchTol);
    if mSi>0&&dist<bestDist, bestDist=dist; matchSi=mSi; end
end
end


% =========================================================================
%  WALL / FLOOR ROTATION  (down-chute and up-chute variants)
% =========================================================================
function [vertsR,centR,phi,pivotOut,axisOut,q_rot]=wallRotation_q(vertsC,centC,wallIdx)
vertsR=[]; centR=[]; phi=0; pivotOut=[0,0,0]; axisOut=[0,0,1]; q_rot=[];
if isempty(wallIdx), return; end
wallV=vertsC(wallIdx,:); [~,pivLoc]=max(wallV(:,1)); pivot=wallV(pivLoc,:); pivotOut=pivot;
hullIdx3D=p4_hullIdx(vertsC); if numel(hullIdx3D)<3, return; end
mTol=max(range(vertsC(:,1)),range(vertsC(:,2)))*1e-3;
hXY_raw=vertsC(hullIdx3D,1:2); hXY=uniquetol(hXY_raw,mTol,'ByRows',true,'DataScale',1);
if size(hXY,1)<3, return; end
try ord2D=convhull(hXY(:,1),hXY(:,2)); catch; return; end
ord2D=ord2D(1:end-1); pts2D=hXY(ord2D,:); M=size(pts2D,1); if M<2, return; end
pivXY=pivot(1:2); dists=vecnorm(pts2D-pivXY,2,2); [~,pivLoc2D]=min(dists);
prevLoc=mod(pivLoc2D-2,M)+1; nextLoc=mod(pivLoc2D,M)+1;
prevXY=pts2D(prevLoc,:); nextXY=pts2D(nextLoc,:);
if nextXY(1)>=prevXY(1), neighXY=nextXY; else, neighXY=prevXY; end
dX=neighXY(1)-pivXY(1); dY=neighXY(2)-pivXY(2);
phi=atan2(-dY,dX); if abs(phi)<1e-8, return; end
q_rot=q_fromAxisAngle([0,0,1],phi); axisOut=[0,0,1];
[vertsR,centR]=q_rotateCloud(q_rot,vertsC,pivot,centC);
end

function [vertsR,centR,phi,pivotOut,axisOut,q_rot]=floorRotation_q(vertsC,centC,floorIdx)
vertsR=[]; centR=[]; phi=0; pivotOut=[0,0,0]; axisOut=[0,1,0]; q_rot=[];
if isempty(floorIdx), return; end
floorV=vertsC(floorIdx,:); [~,pivLoc]=max(floorV(:,1)); pivot=floorV(pivLoc,:); pivotOut=pivot;
hullIdx3D=p4_hullIdx(vertsC); if numel(hullIdx3D)<3, return; end
mTol=max(range(vertsC(:,1)),range(vertsC(:,3)))*1e-3;
hXZ_raw=vertsC(hullIdx3D,[1 3]); hXZ=uniquetol(hXZ_raw,mTol,'ByRows',true,'DataScale',1);
if size(hXZ,1)<3, return; end
try ord2D=convhull(hXZ(:,1),hXZ(:,2)); catch; return; end
ord2D=ord2D(1:end-1); pts2D=hXZ(ord2D,:); M=size(pts2D,1); if M<2, return; end
pivXZ=pivot([1 3]); dists=vecnorm(pts2D-pivXZ,2,2); [~,pivLoc2D]=min(dists);
prevLoc=mod(pivLoc2D-2,M)+1; nextLoc=mod(pivLoc2D,M)+1;
prevXZ=pts2D(prevLoc,:); nextXZ=pts2D(nextLoc,:);
if nextXZ(1)>=prevXZ(1), neighXZ=nextXZ; else, neighXZ=prevXZ; end
dX=neighXZ(1)-pivXZ(1); dZ=neighXZ(2)-pivXZ(2);
phi=atan2(-dZ,dX); if abs(phi)<1e-8, return; end
q_rot=q_fromAxisAngle([0,1,0],-phi); axisOut=[0,1,0];
[vertsR,centR]=q_rotateCloud(q_rot,vertsC,pivot,centC);
end

% ---- Up-chute variants -------------------------------------------------
% Mirror of wallRotation_q / floorRotation_q: pivots on the hull vertex
% furthest AGAINST the direction of travel (min-X instead of max-X), since
% up-chute tipping happens over the opposite edge. The prev/next neighbor
% comparison direction is flipped to match (the "next" edge to rotate onto
% is now on the other side of the pivot).
function [vertsR,centR,phi,pivotOut,axisOut,q_rot]=wallRotationUp_q(vertsC,centC,wallIdx)
vertsR=[]; centR=[]; phi=0; pivotOut=[0,0,0]; axisOut=[0,0,1]; q_rot=[];
if isempty(wallIdx), return; end
wallV=vertsC(wallIdx,:); [~,pivLoc]=min(wallV(:,1)); pivot=wallV(pivLoc,:); pivotOut=pivot;
hullIdx3D=p4_hullIdx(vertsC); if numel(hullIdx3D)<3, return; end
mTol=max(range(vertsC(:,1)),range(vertsC(:,2)))*1e-3;
hXY_raw=vertsC(hullIdx3D,1:2); hXY=uniquetol(hXY_raw,mTol,'ByRows',true,'DataScale',1);
if size(hXY,1)<3, return; end
try ord2D=convhull(hXY(:,1),hXY(:,2)); catch; return; end
ord2D=ord2D(1:end-1); pts2D=hXY(ord2D,:); M=size(pts2D,1); if M<2, return; end
pivXY=pivot(1:2); dists=vecnorm(pts2D-pivXY,2,2); [~,pivLoc2D]=min(dists);
prevLoc=mod(pivLoc2D-2,M)+1; nextLoc=mod(pivLoc2D,M)+1;
prevXY=pts2D(prevLoc,:); nextXY=pts2D(nextLoc,:);
if nextXY(1)<=prevXY(1), neighXY=nextXY; else, neighXY=prevXY; end
dX=neighXY(1)-pivXY(1); dY=neighXY(2)-pivXY(2);
phi=-atan2(-dY,dX);
if abs(phi)<1e-8, return; end
q_rot=q_fromAxisAngle([0,0,1],phi); axisOut=[0,0,1];
[vertsR,centR]=q_rotateCloud(q_rot,vertsC,pivot,centC);
end

function [vertsR,centR,phi,pivotOut,axisOut,q_rot]=floorRotationUp_q(vertsC,centC,floorIdx)
vertsR=[]; centR=[]; phi=0; pivotOut=[0,0,0]; axisOut=[0,1,0]; q_rot=[];
if isempty(floorIdx), return; end
floorV=vertsC(floorIdx,:); [~,pivLoc]=min(floorV(:,1)); pivot=floorV(pivLoc,:); pivotOut=pivot;
hullIdx3D=p4_hullIdx(vertsC); if numel(hullIdx3D)<3, return; end
mTol=max(range(vertsC(:,1)),range(vertsC(:,3)))*1e-3;
hXZ_raw=vertsC(hullIdx3D,[1 3]); hXZ=uniquetol(hXZ_raw,mTol,'ByRows',true,'DataScale',1);
if size(hXZ,1)<3, return; end
try ord2D=convhull(hXZ(:,1),hXZ(:,2)); catch; return; end
ord2D=ord2D(1:end-1); pts2D=hXZ(ord2D,:); M=size(pts2D,1); if M<2, return; end
pivXZ=pivot([1 3]); dists=vecnorm(pts2D-pivXZ,2,2); [~,pivLoc2D]=min(dists);
prevLoc=mod(pivLoc2D-2,M)+1; nextLoc=mod(pivLoc2D,M)+1;
prevXZ=pts2D(prevLoc,:); nextXZ=pts2D(nextLoc,:);
if nextXZ(1)<=prevXZ(1), neighXZ=nextXZ; else, neighXZ=prevXZ; end
dX=neighXZ(1)-pivXZ(1); dZ=neighXZ(2)-pivXZ(2);
phi=atan2(-dZ,dX);
if abs(phi)<1e-8, return; end
q_rot=q_fromAxisAngle([0,1,0],phi);
axisOut=[0,1,0];
[vertsR,centR]=q_rotateCloud(q_rot,vertsC,pivot,centC);
end


% =========================================================================
%  INIT POSE / RESEAT / QUAT MATCH
% =========================================================================
function [vertsW,centW]=initPose(si,allPoses,stableIdx)
pos=allPoses(stableIdx(si)); vertsW=pos.verticesWorld; centW=pos.centroidWorld(:)';
fCI=pos.floorContactVertIdx(:); wCI=pos.wallContactVertIdx(:);
fZ=min(vertsW(fCI,3)); vertsW(:,3)=vertsW(:,3)-fZ; centW(3)=centW(3)-fZ;
wY=max(vertsW(wCI,2)); vertsW(:,2)=vertsW(:,2)-wY; centW(2)=centW(2)-wY;
end

function [vertsOut,centOut,floorIdx,wallIdx]=reseat(vertsC,centC)
span=max(max(vertsC,[],1)-min(vertsC,[],1)); cTol=max(0.05,span*0.01);
hullIdx=p4_hullIdx(vertsC);
zVals=vertsC(hullIdx,3); zMin=min(zVals);
vertsC(:,3)=vertsC(:,3)-zMin; centC(3)=centC(3)-zMin;
yVals=vertsC(hullIdx,2); yMax=max(yVals);
if ~isnan(yMax), vertsC(:,2)=vertsC(:,2)-yMax; centC(2)=centC(2)-yMax; end
zVals2=vertsC(hullIdx,3); yVals2=vertsC(hullIdx,2);
floorIdx=hullIdx(abs(zVals2-0)<=cTol); wallIdx=hullIdx(abs(yVals2-0)<=cTol);
vertsOut=vertsC; centOut=centC;
end

function [matchSi,bestDist]=matchQuatComposed(q_query,refQuatsAll,tol)
matchSi=0; bestDist=inf;
for sj=1:numel(refQuatsAll)
    qList=refQuatsAll{sj}; if isempty(qList), continue; end
    for kk=1:size(qList,1)
        d=q_geodesic(q_query,qList(kk,:)); if d<bestDist, bestDist=d; matchSi=sj; end
    end
end
if bestDist>tol, matchSi=0; end
end

function idx=p4_hullIdx(vertsC)
try hf=convhull(vertsC,'Simplify',true); idx=unique(hf(:));
catch; idx=(1:size(vertsC,1))'; end
end


% =========================================================================
%  QUATERNION PRIMITIVES
% =========================================================================
function q=q_fromAxisAngle(axis,phi)
axis=axis(:)'/norm(axis(:)); s=sin(phi/2);
q=[cos(phi/2),s*axis(1),s*axis(2),s*axis(3)];
end
function q=q_compose(q1,q2)
w1=q1(1);x1=q1(2);y1=q1(3);z1=q1(4); w2=q2(1);x2=q2(2);y2=q2(3);z2=q2(4);
q=[w2*w1-x2*x1-y2*y1-z2*z1,w2*x1+x2*w1+y2*z1-z2*y1, ...
   w2*y1-x2*z1+y2*w1+z2*x1,w2*z1+x2*y1-y2*x1+z2*w1];
end
function [Vout,cout]=q_rotateCloud(q,V,pivot,cent)
Rmat=q_toRotm(q); p=pivot(:)';
Vout=(Rmat*(V-p)')'+p; cout=(Rmat*(cent(:)'-p)')'+p;
end
function R=q_toRotm(q)
w=q(1);x=q(2);y=q(3);z=q(4);
R=[1-2*(y^2+z^2),2*(x*y-z*w),2*(x*z+y*w);
   2*(x*y+z*w),1-2*(x^2+z^2),2*(y*z-x*w);
   2*(x*z-y*w),2*(y*z+x*w),1-2*(x^2+y^2)];
end
function q=q_fromRotm(R)
tr=R(1,1)+R(2,2)+R(3,3);
if tr>0
    s=0.5/sqrt(tr+1); q=[0.25/s,(R(3,2)-R(2,3))*s,(R(1,3)-R(3,1))*s,(R(2,1)-R(1,2))*s];
elseif (R(1,1)>R(2,2))&&(R(1,1)>R(3,3))
    s=2*sqrt(1+R(1,1)-R(2,2)-R(3,3));
    q=[(R(3,2)-R(2,3))/s,0.25*s,(R(1,2)+R(2,1))/s,(R(1,3)+R(3,1))/s];
elseif R(2,2)>R(3,3)
    s=2*sqrt(1+R(2,2)-R(1,1)-R(3,3));
    q=[(R(1,3)-R(3,1))/s,(R(1,2)+R(2,1))/s,0.25*s,(R(2,3)+R(3,2))/s];
else
    s=2*sqrt(1+R(3,3)-R(1,1)-R(2,2));
    q=[(R(2,1)-R(1,2))/s,(R(1,3)+R(3,1))/s,(R(2,3)+R(3,2))/s,0.25*s];
end
q=q/norm(q);
end
function d=q_geodesic(q1,q2)
dp=abs(dot(q1(:),q2(:))); dp=min(dp,1.0); d=2*acos(dp);
end


% =========================================================================
%  GEOMETRY PIPELINE
% =========================================================================
function c=centroidOfPolyhedron(vertex,faces)
v1=vertex(faces(:,1),:); v2=vertex(faces(:,2),:); v3=vertex(faces(:,3),:);
vec1=v2-v1; vec2=v3-v1; ta=0.5*cross(vec1,vec2);
A=sqrt(ta(:,1).^2+ta(:,2).^2+ta(:,3).^2); tot=sum(A); ct=(v1+v2+v3)/3;
c=[sum(A.*ct(:,1)),sum(A.*ct(:,2)),sum(A.*ct(:,3))]/tot;
end

function [planeVerts,planeEqs]=findFloorPlanes(partVertices,convexHullFaces,planeTol)
planeVerts={}; planeEqs=[];
for fi=1:size(convexHullFaces,1)
    v1=partVertices(convexHullFaces(fi,1),:); v2=partVertices(convexHullFaces(fi,2),:);
    v3=partVertices(convexHullFaces(fi,3),:);
    nv=cross(v2-v1,v3-v1); if norm(nv)<1e-10, continue; end
    D=dot(nv,v1); res=partVertices*nv'-D; onP=find(abs(res)<planeTol);
    if isempty(onP), continue; end
    sP=sort(onP(:))'; dup=false;
    for k=1:numel(planeVerts), if isequal(sort(planeVerts{k}),sP), dup=true; break; end; end
    if ~dup, planeVerts{end+1}=sP; planeEqs(end+1,:)=[nv,D]; end %#ok<AGROW>
end
end

function [rawCandidates,planeQuats]=thetaFromHullEdges( ...
    partVertices,convexHullFaces,chullVertexIdx, ...
    restingPlaneVerts,restingPlaneEqs,centroidCoords,wallTol,planeTol)
rawCandidates=struct('floorPlaneIdx',{},'wallSide',{},'wallContactVertIdx',{}, ...
    'floorContactVertIdx',{},'theta',{},'verticesWorld',{},'centroidWorld',{}, ...
    'floorZ',{},'centroidHeightAboveFloor',{});
numP=numel(restingPlaneVerts); planeQuats=repmat({[1,0,0,0]},1,numP);
edgeTol=wallTol; neighbourDeg=5; neighbourN=20; fallbackSamples=720;
for planeIdx=1:numP
    planeEq=restingPlaneEqs(planeIdx,:); planeNorm=planeEq(1:3)/norm(planeEq(1:3));
    supportIdx=restingPlaneVerts{planeIdx}; if numel(supportIdx)<2, continue; end
    faceN=planeNorm(:); targetN=[0;0;-1]; dotFT=dot(faceN,targetN);
    if abs(dotFT-1)<1e-8, q_align=[1,0,0,0];
    elseif abs(dotFT+1)<1e-8
        perp=null(faceN'); ax=perp(:,1)/norm(perp(:,1)); q_align=q_fromAxisAngle(ax',pi);
    else
        ax=cross(faceN,targetN); ax=ax/norm(ax); ang=acos(max(-1,min(1,dotFT)));
        q_align=q_fromAxisAngle(ax',ang);
    end
    planeQuats{planeIdx}=q_align;
    R_floor=q_toRotm(q_align); C0=centroidCoords(:)';
    vertsRot=(R_floor*(partVertices'-C0'))'+C0; centRot=C0;
    suppRot=vertsRot(supportIdx,:); floorZ=mean(suppRot(:,3));
    cHAF=centRot(3)-floorZ; if cHAF<1e-6, continue; end
    hullVerts2D=vertsRot(chullVertexIdx,1:2);
    try ord2D=convhull(hullVerts2D(:,1),hullVerts2D(:,2),'Simplify',false); catch; ord2D=[]; end
    exactThetas=[];
    if ~isempty(ord2D)
        edgePts=hullVerts2D(ord2D(1:end-1),:); M=size(edgePts,1);
        for ei=1:M
            p0=edgePts(ei,:); p1=edgePts(mod(ei,M)+1,:);
            dXY=p1-p0; eL=norm(dXY); if eL<1e-10, continue; end
            dXYn=dXY/eL; dp=hullVerts2D-p0;
            dPerp=abs(dp(:,1)*dXYn(2)-dp(:,2)*dXYn(1));
            if sum(dPerp<=edgeTol)<2, continue; end
            tExact=atan2(dXYn(1),dXYn(2));
            for tD=[0,pi], exactThetas(end+1)=mod(tExact+tD,2*pi); end %#ok<AGROW>
        end
    end
    thetaCands=exactThetas(:)'; dN=deg2rad(neighbourDeg);
    for et=exactThetas
        neigh=et+linspace(-dN,dN,2*neighbourN+1); thetaCands=[thetaCands,neigh]; %#ok<AGROW>
    end
    fb=linspace(0,2*pi,fallbackSamples+1); fb=fb(1:end-1);
    thetaCands=[thetaCands,fb]; thetaCands=mod(thetaCands,2*pi);
    thetaCands=sort(unique(round(thetaCands/deg2rad(0.01))*deg2rad(0.01)));
    thetaWallKeys={};
    for ti=1:numel(thetaCands)
        theta=thetaCands(ti); allVW=rotatePtsAroundZ(vertsRot,centRot,theta);
        fZc=min(allVW(:,3)); fResid=abs(allVW(:,3)-fZc); onFloor=find(fResid<planeTol);
        if numel(onFloor)<2, continue; end
        fPts2D=allVW(onFloor,1:2); cXY=[centRot(1),centRot(2)];
        if numel(onFloor)>=3
            try hI=convhull(fPts2D(:,1),fPts2D(:,2));
                floorOK=inpolygon(cXY(1),cXY(2),fPts2D(hI,1),fPts2D(hI,2));
            catch; floorOK=false; end
        else
            floorOK=(cXY(1)>=min(fPts2D(:,1))-wallTol)&&(cXY(1)<=max(fPts2D(:,1))+wallTol);
        end
        if ~floorOK, continue; end
        projPts=allVW(chullVertexIdx,1:2); yWall=max(projPts(:,2));
        onWL=find(abs(projPts(:,2)-yWall)<wallTol);
        if numel(onWL)>=2
            onWall=chullVertexIdx(onWL); wPts=allVW(onWall,1:2);
            wPts_u=uniquetol(wPts,wallTol,'ByRows',true);
            if size(wPts_u,1)>=2
                wX=wPts_u(:,1); cX=centRot(1);
                if cX>=min(wX)-wallTol&&cX<=max(wX)+wallTol
                    wc=sort(onWall(:)'); fc=sort(onFloor(:)');
                    key=sprintf('W_%s__F_%s',mat2str(wc),mat2str(fc));
                    if ~any(strcmp(thetaWallKeys,key))
                        thetaWallKeys{end+1}=key; %#ok<AGROW>
                        pos.floorPlaneIdx=planeIdx; pos.wallSide=1;
                        pos.wallContactVertIdx=wc; pos.floorContactVertIdx=fc;
                        pos.theta=theta; pos.verticesWorld=allVW;
                        pos.centroidWorld=centRot; pos.floorZ=floorZ;
                        pos.centroidHeightAboveFloor=cHAF;
                        rawCandidates(end+1)=pos; %#ok<AGROW>
                    end
                end
            end
        end
    end
end
end

function pts_out=rotatePtsAroundZ(pts_in,pivot,theta)
s=sin(theta); c=cos(theta); R=[c,-s,0;s,c,0;0,0,1];
pts_out=(R*(pts_in-pivot)')'+pivot;
end

function merged=mergeByTheta(rawCandidates,thetaTol_deg,planeMergeTol,partVertices) %#ok<INUSD>
merged=struct('floorPlaneIdx',{},'wallSide',{},'wallContactVertIdx',{}, ...
    'floorContactVertIdx',{},'theta',{},'verticesWorld',{},'centroidWorld',{}, ...
    'floorZ',{},'centroidHeightAboveFloor',{});
if isempty(rawCandidates), return; end
floorIds=[rawCandidates.floorPlaneIdx]; wallSides=[rawCandidates.wallSide];
thetas=[rawCandidates.theta];
groups=unique([floorIds(:),wallSides(:)],'rows');
for gi=1:size(groups,1)
    fp=groups(gi,1); ws=groups(gi,2);
    sel=find(floorIds==fp&wallSides==ws); if isempty(sel), continue; end
    groupThetas=thetas(sel(:))'; tTol=deg2rad(thetaTol_deg); used=false(size(sel));
    for ii=1:numel(sel)
        if used(ii), continue; end
        t0=rawCandidates(sel(ii)).theta; dt=angle(exp(1i*(groupThetas-t0)));
        nearby=abs(dt)<tTol; clIdx=sel(nearby); used(nearby)=true;
        nW=arrayfun(@(i)numel(rawCandidates(i).wallContactVertIdx),clIdx);
        [~,best]=max(nW); rep=rawCandidates(clIdx(best));
        aW=[]; aF=[];
        for kk=1:numel(clIdx)
            aW=union(aW,rawCandidates(clIdx(kk)).wallContactVertIdx);
            aF=union(aF,rawCandidates(clIdx(kk)).floorContactVertIdx);
        end
        rep.wallContactVertIdx=aW; rep.floorContactVertIdx=aF;
        rep.theta=atan2(mean(sin([rawCandidates(clIdx).theta])),mean(cos([rawCandidates(clIdx).theta])));
        if rep.theta<0, rep.theta=rep.theta+2*pi; end
        merged(end+1)=rep; %#ok<AGROW>
    end
end
end

function [Qs,omegas,heights,projCentroids,hullProjs]=computeCSA( ...
    allPoses,stableIdx,chullVertexIdx,Rchute,gWorld) %#ok<INUSD>
numS=numel(stableIdx); Qs=zeros(numS,1); omegas=zeros(numS,1);
heights=zeros(numS,1); projCentroids=NaN(numS,3); hullProjs=cell(numS,1);
gWorld=gWorld(:)/norm(gWorld(:));
for si=1:numS
    pos=allPoses(stableIdx(si)); C0=pos.centroidWorld(:)';
    vW=(Rchute*(pos.verticesWorld'-C0'))'+C0; cW=C0;
    supIdx=union(pos.floorContactVertIdx,pos.wallContactVertIdx);
    supV=vW(supIdx,:); if size(supV,1)<3, continue; end
    fV=vW(pos.floorContactVertIdx,:); fM=mean(fV,1);
    if size(fV,1)>=3, [~,~,V]=svd(fV-fM); fN=V(:,3)'; else, fN=-gWorld(:)'; end
    if dot(fN,gWorld')>0, fN=-fN; end; fN=fN/norm(fN);
    rd=[0,0,-1]; dn=dot(rd,fN); if abs(dn)<1e-10, continue; end
    t=dot(fM-cW,fN)/dn; if t<0, continue; end
    h=t; if h<1e-10, continue; end
    heights(si)=h; projCentroids(si,:)=cW+t*rd;
    gC=supV*gWorld; mG=min(gC); pp=mG*gWorld';
    sProj=zeros(size(supV));
    for kk=1:size(supV,1)
        p_=supV(kk,:); d_=dot((p_-pp),gWorld'); sProj(kk,:)=p_-d_*gWorld';
    end
    sProj=uniquetol(sProj,1e-5,'ByRows',true); if size(sProj,1)<3, continue; end
    if abs(gWorld(1))<0.9, t2=[1;0;0]; else, t2=[0;1;0]; end
    e1b=cross(gWorld,t2); e1b=e1b/norm(e1b); e2b=cross(gWorld,e1b); e2b=e2b/norm(e2b);
    s2D=[sProj*e1b,sProj*e2b];
    try hI=convhull(s2D(:,1),s2D(:,2)); catch; continue; end
    hull3D=sProj(hI(1:end-1),:); nH=size(hull3D,1); if nH<3, continue; end
    Vs=hull3D-cW; norms_Vs=vecnorm(Vs,2,2); if any(norms_Vs<1e-10), continue; end
    Vs=Vs./norms_Vs; omega=0;
    for k=1:nH-2
        a_=Vs(1,:); b_=Vs(k+1,:); c_=Vs(k+2,:);
        nv=abs(dot(a_,cross(b_,c_))); dv=1+dot(a_,b_)+dot(b_,c_)+dot(a_,c_);
        if abs(dv)<1e-14, continue; end
        omega=omega+2*atan2(nv,dv);
    end
    omegas(si)=omega; Qs(si)=omega/h; hullProjs{si}=hull3D;
end
end

function [stIdxOut,QsOut,omOut,htOut,mGroups]=mergeByCSAValues( ...
    stableIdx,Qs,omegas,heights,omegaTol,hTol,allPoses,Rchute) %#ok<INUSD>
n=numel(stableIdx); clouds=cell(n,1);
for si=1:n
    pos=allPoses(stableIdx(si)); C0=pos.centroidWorld(:)';
    vW=(Rchute*(pos.verticesWorld'-C0'))'+C0; clouds{si}=vW-mean(vW,1);
end
used=false(n,1); stIdxOut=[]; QsOut=[]; omOut=[]; htOut=[]; mGroups={};
for ii=1:n
    if used(ii), continue; end
    gM=false(n,1); gM(ii)=true;
    for jj=ii+1:n
        if used(jj), continue; end
        if cloudsMatch(clouds{ii},clouds{jj},hTol), gM(jj)=true; end
    end
    gIdx=find(gM); used(gM)=true;
    [~,bL]=max(Qs(gIdx)); bG=gIdx(bL);
    stIdxOut(end+1)=stableIdx(bG); QsOut(end+1)=sum(Qs(gIdx)); %#ok<AGROW>
    omOut(end+1)=omegas(bG); htOut(end+1)=heights(bG); %#ok<AGROW>
    mGroups{end+1}=stableIdx(gIdx(:)'); %#ok<AGROW>
end
stIdxOut=stIdxOut(:); QsOut=QsOut(:); omOut=omOut(:); htOut=htOut(:);
if numel(stableIdx)~=numel(stIdxOut)
    fprintf('  CSA merge: %d -> %d\n',numel(stableIdx),numel(stIdxOut));
end
end

function tf=cloudsMatch(A,B,tol)
if isempty(A)||isempty(B), tf=false; return; end
tf=allPtsMatched(A,B,tol)&&allPtsMatched(B,A,tol);
end
function tf=allPtsMatched(src,tgt,tol)
tf=true;
for k=1:size(src,1), if min(vecnorm(tgt-src(k,:),2,2))>tol, tf=false; return; end; end
end

function [ratioWall,ratioFloor,transitions,momentArmGeo]=computeTransitionRatios( ...
    allPoses,stableIdx,Rchute,slideDir,floorNorm_c,wallNorm_c,planeTol) %#ok<INUSL>
% NOTE: ratioWall/ratioFloor are still computed here for the stability
% gating (Qs zeroing) and bubble-color logic — they are just no longer
% used as edge cost inside resolveTransitionChain, which now uses a
% uniform angle-based cost for every transition type.
numS=numel(stableIdx);
ratioWall=zeros(numS,1); ratioFloor=zeros(numS,1);
transitions=false(numS,1);
eGeo=struct('pivotWall_c',[],'pivotFloor_c',[],'centC',[], ...
    'l_A_wall',0,'l_w',0,'l_A_floor',0,'l_f',0);
momentArmGeo=repmat(eGeo,numS,1);
for si=1:numS
    pos=allPoses(stableIdx(si)); vC=pos.verticesWorld; cC=pos.centroidWorld(:)';
    fCI=pos.floorContactVertIdx; wCI=pos.wallContactVertIdx;
    fZ=min(vC(fCI,3)); vC(:,3)=vC(:,3)-fZ; cC(3)=cC(3)-fZ;
    wY=max(vC(wCI,2)); vC(:,2)=vC(:,2)-wY; cC(2)=cC(2)-wY;
    momentArmGeo(si).centC=cC; fV=vC(fCI,:); wV=vC(wCI,:);
    if ~isempty(wV)
        [~,iW]=max(wV(:,1)); pW=wV(iW,:);
        l_A_w=abs(cC(2)); l_w=abs(cC(1)-pW(1));
        momentArmGeo(si).pivotWall_c=pW; momentArmGeo(si).l_A_wall=l_A_w; momentArmGeo(si).l_w=l_w;
        if l_w>planeTol, ratioWall(si)=l_A_w/l_w; end
    end
    if ~isempty(fV)
        [~,iF]=max(fV(:,1)); pF=fV(iF,:);
        l_A_f=abs(cC(3)); l_f=abs(cC(1)-pF(1));
        momentArmGeo(si).pivotFloor_c=pF; momentArmGeo(si).l_A_floor=l_A_f; momentArmGeo(si).l_f=l_f;
        if l_f>planeTol, ratioFloor(si)=l_A_f/l_f; end
    end
end
end

function v=normaliseVec(v)
n=norm(v); if n>0, v=v/n; end
end
