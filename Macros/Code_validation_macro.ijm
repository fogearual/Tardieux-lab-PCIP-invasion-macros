// ======================================================================
//  Counting control - automated counts vs manual counts
// ======================================================================
//
//  Same folder layout and the SAME detection settings as the main
//  pipeline, so the counts produced here are identical to it.
//
//  Sampling, independently for every experiment:
//    - conditions are shuffled, then one random .czi is drawn from each
//      in turn; a condition is drawn again only once every other one has
//      met its share (target / number of conditions)
//    - drawing stops once all three targets are met; the image being
//      processed is always finished
//    - a class is recorded only while its own target is unmet, so the
//      sampled sets are nested: cells in the first images, green in more,
//      red in all of them
//
//  Output, in <selected folder>/Counting code control/:
//    <experiment>/   one 6-slice stack per sampled image
//       1 blue raw          count the cells here
//       2 green raw         count the parasites here
//       3 red raw
//       4 merge G+R raw     count the double positives here
//       5 merge G+R annotated (green and red)
//       6 blue annotated
//    Counting control <name>.xls   auto counts, empty manual columns,
//                                  error columns computed on entry
//
//  Count on slices 1 to 4 and only then compare with 5 and 6: seeing the
//  outlines before counting biases the count and voids the comparison.
// ======================================================================


// ======================================================================
//  1. SAMPLING
// ======================================================================
targetCells = 500;   // stop recording cells once this many are reached
targetGreen = 300;   // idem for green parasites
targetRed   = 150;   // idem for double positives

// The one constraint on the draw: at most this many images from a single
// condition, so no condition dominates the control. Beyond that the draw
// is purely random over the whole experiment.
maxPerCondition = 4;

randomSeed  = 0;     // 0 = new draw every run; any other value = reproducible

outFolderName = "Counting code control";


// ======================================================================
//  2. DETECTION  -  the settings to adjust first
// ======================================================================

// --- GREEN: parasites ---------------------------------------------------
// Threshold = Triangle threshold x factor. Raising the factor demands more
// fluorescence, so less background is picked up. The value used is logged.
greenThresholdFactor = 4.25;
greenMinArea         = 7.25;   // um2; a single tachyzoite is 6-15 um2
greenMaxArea         = 62.5;   // um2; above this, cluster or debris
greenWatershed       = true;   // split touching parasites; set false if a
                               // single elongated parasite is cut in two

// --- RED: colocalisation marker -----------------------------------------
redThresholdFactor = 2.8;      // same logic as green
minOverlapPercent  = 92;       // % of a green object that must be red-positive
                               // to count as double positive

// --- BLUE: nuclei -------------------------------------------------------
// An object larger than nucleiDoubletFactor x the field median area is
// treated as several overlapping cells. The median is recomputed for every
// image, so the criterion adapts to the field.
splitBigNuclei      = true;
nucleiDoubletFactor = 1.7;
nucleiMaxSplit      = 4;       // max cells assigned to one object


// ======================================================================
//  3. DISPLAY
// ======================================================================
displayBoost = 1.25;   // brightens red and green on screen only;
                       // 1.0 = unchanged. Never affects any measurement.
labelScale   = 1.0;    // object number size
recapScale   = 1.0;    // size of the summary text at the bottom of images
strokeScale  = 1.65;   // outline thickness

// true keeps one window open per image; with 8 conditions x 50 images
// this exhausts memory. Keep false for batch runs.
keepWindowsOpen = false;


// ======================================================================
//  4. CHANNEL -> COLOUR
// ======================================================================
// Channel colours are read from the CZI metadata. The fixed numbers below
// are used only if that fails, or if autoDetectChannels is set to false.
autoDetectChannels = true;

// Channel number as produced by "Split Channels". 0 = channel absent.
chDIC   = 1;
chRED   = 2;
chGREEN = 3;
chBLUE  = 4;


// ======================================================================
//  5. PREPROCESSING  -  change only if you know why
// ======================================================================
// Multiply and Enhance Contrast do not change detection: the automatic
// threshold is recomputed on the rescaled histogram. Sensitivity is
// controlled by the factors in section 2.
nucleiMultiply        = 2.5;
nucleiBgRolling       = 80;             // Subtract Background radius, px
nucleiBlurSigma       = 1.5;
nucleiThresholdMethod = "Triangle";     // "Triangle", "Otsu", "Li", "Yen"...
nucleiWatershed       = true;
nucleiMinArea         = 105;            // um2
nucleiMaxArea         = "Infinity";
nucleiCirc            = "0.20-1.00";

greenMultiply         = 2.0;
greenBgRolling        = 30;
greenBlurSigma        = 0.8;
greenThresholdMethod  = "Triangle";
greenCirc             = "0.10-1.00";

redMultiply           = 2.0;
redBgRolling          = 30;
redBlurSigma          = 0.8;
redThresholdMethod    = "Triangle";

enhanceSaturated      = 0.10;


// ----------------------------------------------------------------------
// Globals shared with the functions. In the ImageJ macro language a
// function only sees its parameters and variables declared with "var".
// ----------------------------------------------------------------------
var vExp   = newArray(0);      // one entry per sampled image
var vCond  = newArray(0);
var vImg   = newArray(0);
var vCells = newArray(0);      // -1 when the class was not recorded
var vGreen = newArray(0);
var vRed   = newArray(0);

var mAcc = newArray(0);        // recursive scan
var cAcc = newArray(0);
var pAcc = newArray(0);

var dcDIC = 0;                 // detected channels
var dcRED = 0;
var dcGREEN = 0;
var dcBLUE = 0;

var condsAll = newArray(0);
var rootName = "";
var singleManip = 0;


// ======================================================================
//  PATHS
// ======================================================================
sep = File.separator;
if (randomSeed != 0) random("seed", randomSeed);

rootFolder = getDirectory("Select the folder containing the experiments");
scanDir(rootFolder, "", "");

print("\\Clear");
if (pAcc.length == 0) {
    print("Contents of " + rootFolder + " :");
    dbg = getFileList(rootFolder);
    for (a = 0; a < dbg.length; a++) print("   " + dbg[a]);
    exit("No .czi file found under " + rootFolder);
}

rootName   = cleanName(rootFolder);
outputRoot = rootFolder + outFolderName + sep;
File.makeDirectory(outputRoot);
tmpZip = outputRoot + "_tmp_nuclei.zip";

// ---- experiments, alphabetical ----
manipList = newArray(0);
for (t = 0; t < pAcc.length; t++) {
    known = 0;
    for (u = 0; u < manipList.length; u++) if (manipList[u] == mAcc[t]) known = 1;
    if (known == 0) manipList = Array.concat(manipList, mAcc[t]);
}
Array.sort(manipList);
nM = manipList.length;
singleManip = 0;
if (nM == 1) { if (manipList[0] == "-") singleManip = 1; }

// ---- conditions, controls first then increasing concentration ----
for (t = 0; t < pAcc.length; t++) {
    known = 0;
    for (u = 0; u < condsAll.length; u++) if (condsAll[u] == cAcc[t]) known = 1;
    if (known == 0) condsAll = Array.concat(condsAll, cAcc[t]);
}
nC = condsAll.length;
sk = newArray(nC);
for (c = 0; c < nC; c++) {
    num = firstNumber(condsAll[c]);
    if (num < 0) num = 99999;
    if (isControlName(condsAll[c]) == 1) sk[c] = num;
    else                                 sk[c] = 1000000000 + num;
}
for (a = 0; a < nC - 1; a++) {
    mn = a;
    for (b = a + 1; b < nC; b++) if (sk[b] < sk[mn]) mn = b;
    if (mn != a) {
        sw = condsAll[a]; condsAll[a] = condsAll[mn]; condsAll[mn] = sw;
        sv = sk[a];       sk[a]       = sk[mn];       sk[mn]       = sv;
    }
}

print(pAcc.length + " image(s) available, " + nM + " experiment(s), " + nC + " condition(s)");
print("targets per experiment: " + targetCells + " cells, " + targetGreen +
      " green, " + targetRed + " double positives");
print("");


// ======================================================================
//  SAMPLING AND COUNTING
// ======================================================================
setBatchMode(true);

for (m = 0; m < nM; m++) {

    manip = manipList[m];
    expFolder = outputRoot + manip + sep;
    if (singleManip == 1) expFolder = outputRoot;
    File.makeDirectory(expFolder);

    print("=== " + manip + " ===");

    // Every image of this experiment, in random order. The only constraint
    // is maxPerCondition: beyond that the draw is purely random.
    pool = newArray(0);
    for (t = 0; t < pAcc.length; t++) {
        if (mAcc[t] != manip) continue;
        pool = Array.concat(pool, t);
    }
    if (pool.length == 0) continue;
    for (a = pool.length - 1; a > 0; a--) {
        b = floor(random() * (a + 1));
        sw = pool[a]; pool[a] = pool[b]; pool[b] = sw;
    }

    drawnC = newArray(nC);       // images already drawn from each condition
    totCells = 0; totGreen = 0; totRed = 0;
    nTaken = 0;

    for (k = 0; k < pool.length; k++) {

        needCells = 0; if (totCells < targetCells) needCells = 1;
        needGreen = 0; if (totGreen < targetGreen) needGreen = 1;
        needRed   = 0; if (totRed   < targetRed)   needRed   = 1;
        if (needCells + needGreen + needRed == 0) k = pool.length;
        if (k >= pool.length) continue;

        t = pool[k];
        ci = -1;
        for (u = 0; u < nC; u++) if (condsAll[u] == cAcc[t]) ci = u;
        if (ci < 0) continue;
        if (drawnC[ci] >= maxPerCondition) continue;
        drawnC[ci] = drawnC[ci] + 1;

        path = pAcc[t];
        base = File.getName(path);
        base = substring(base, 0, lengthOf(base) - 4);
        print("[" + manip + " / " + cAcc[t] + "] " + base);

        processImage(path, base, expFolder);

        // A class is recorded only while its own target is unmet, so the
        // sampled sets are nested.
        recCells = -1; recGreen = -1; recRed = -1;
        if (needCells == 1) { recCells = pCells; totCells = totCells + pCells; }
        if (needGreen == 1) { recGreen = pGreen; totGreen = totGreen + pGreen; }
        if (needRed   == 1) { recRed   = pRed;   totRed   = totRed + pRed; }

        vExp   = Array.concat(vExp, manip);
        vCond  = Array.concat(vCond, cAcc[t]);
        vImg   = Array.concat(vImg, base);
        vCells = Array.concat(vCells, recCells);
        vGreen = Array.concat(vGreen, recGreen);
        vRed   = Array.concat(vRed, recRed);
        nTaken = nTaken + 1;
    }

    print("   " + nTaken + " image(s) sampled -> " + totCells + " cells, " +
          totGreen + " green, " + totRed + " double positives");
    if (totCells < targetCells)
        print("   NOTE: cells target (" + targetCells + ") not reached");
    if (totGreen < targetGreen)
        print("   NOTE: green target (" + targetGreen + ") not reached");
    if (totRed < targetRed)
        print("   NOTE: double positive target (" + targetRed + ") not reached." +
              " At the high doses there are almost none left; raise" +
              " maxPerCondition or lower targetRed.");
    print("");
}

setBatchMode(false);
if (File.exists(tmpZip)) File.delete(tmpZip);

writeControlBook(outputRoot + "Counting control " + rootName + ".xls");

print("Done. " + vImg.length + " image(s) sampled.");
print("1. count by hand on slices 1 to 4 of the stacks in " + outputRoot);
print("2. type your counts in the Manual columns of the .xls");
print("3. the error columns compute themselves");


// ======================================================================
//  FUNCTIONS - sampling
// ======================================================================

// Recursive scan: the last two folders above a file give the experiment
// and the condition. With a single level, the experiment is "-".
function scanDir(dir, lvl1, lvl2) {
    lst = getFileList(dir);
    Array.sort(lst);
    for (q = 0; q < lst.length; q++) {
        nm = lst[q];
        full = dir + nm;
        if (File.isDirectory(full)) {
            cn = cleanName(nm);
            if (cn == "Results") continue;
            if (cn == outFolderName) continue;
            if (lvl1 == "")      scanDir(full, cn, "");
            else if (lvl2 == "") scanDir(full, lvl1, cn);
            else                 scanDir(full, lvl2, cn);
        } else if (endsWith(toLowerCase(nm), ".czi")) {
            if (lvl2 == "") { mm = "-";  cc = lvl1; }
            else            { mm = lvl1; cc = lvl2; }
            if (cc == "") cc = "-";
            mAcc = Array.concat(mAcc, mm);
            cAcc = Array.concat(cAcc, cc);
            pAcc = Array.concat(pAcc, full);
        }
    }
}



// Channel colour read from the LUT Bio-Formats applies from the CZI
// metadata. Grey = DIC; otherwise the dominant component gives the colour.
// Composite LUTs (cyan, magenta, yellow) go to the colour still free.
function detectChannels(t) {
    dcDIC = 0; dcRED = 0; dcGREEN = 0; dcBLUE = 0;
    selectWindow(t);
    getDimensions(dw, dh, dnc, dns, dnf);
    for (c = 1; c <= dnc; c++) {
        if (dnc > 1) Stack.setChannel(c);
        getLut(lr, lg, lb);
        last = lr.length - 1;
        rr = lr[last];
        gg = lg[last];
        bb = lb[last];
        if (rr == gg && gg == bb) {
            if (dcDIC == 0) dcDIC = c;
        } else if (bb > rr && bb > gg) {
            if (dcBLUE == 0) dcBLUE = c;
        } else if (gg > rr && gg > bb) {
            if (dcGREEN == 0) dcGREEN = c;
        } else if (rr > gg && rr > bb) {
            if (dcRED == 0) dcRED = c;
        } else if (gg == bb && gg > rr) {          // cyan
            if (dcBLUE == 0) dcBLUE = c;
            else if (dcGREEN == 0) dcGREEN = c;
        } else if (rr == bb && rr > gg) {          // magenta
            if (dcRED == 0) dcRED = c;
            else if (dcBLUE == 0) dcBLUE = c;
        } else if (rr == gg && rr > bb) {          // jaune
            if (dcRED == 0) dcRED = c;
            else if (dcGREEN == 0) dcGREEN = c;
        }
    }
}

function isControlName(name) {
    lc = toLowerCase(name);
    if (indexOf(lc, "control") >= 0) return 1;
    if (indexOf(lc, "controle") >= 0) return 1;
    if (indexOf(lc, "ctrl") >= 0) return 1;
    if (indexOf(lc, "temoin") >= 0) return 1;
    return 0;
}

// First number contained in a string, -1 if none.
function firstNumber(s) {
    digits = "0123456789";
    buf = "";
    started = 0;
    stop = 0;
    for (i2 = 0; i2 < lengthOf(s); i2++) {
        if (stop == 1) continue;
        c = substring(s, i2, i2 + 1);
        if (indexOf(digits, c) >= 0) { buf = buf + c; started = 1; }
        else if (started == 1 && (c == "." || c == ",")) { buf = buf + "."; }
        else if (started == 1) stop = 1;
    }
    if (buf == "") return -1;
    if (endsWith(buf, ".")) buf = substring(buf, 0, lengthOf(buf) - 1);
    return parseFloat(buf);
}

// Folder name without trailing separator
function cleanName(s) {
    while (endsWith(s, "/") || endsWith(s, "\\"))
        s = substring(s, 0, lengthOf(s) - 1);
    return File.getName(s);
}


// ======================================================================
//  FUNCTIONS - one image
// ======================================================================
// Counts are returned in pCells / pGreen / pRed, and the 6-slice stack is
// written to outFolder. The segmentation below is identical to the main
// pipeline, so the numbers match it exactly.
var pCells = 0;
var pGreen = 0;
var pRed   = 0;
var tmpZip = "";

function processImage(path, base, outFolder) {

    pCells = 0; pGreen = 0; pRed = 0;
    roiManager("Reset");
    run("Clear Results");

    // ---- Import, max Z projection, channel split ----
    run("Bio-Formats Importer",
        "open=[" + path + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
    original = getTitle();
    run("Z Project...", "projection=[Max Intensity]");
    projected = getTitle();
    selectWindow(original);
    close();

    selectWindow(projected);
    getDimensions(dimW, dimH, nCh, nSl, nFr);

    // ---- Which channel is which colour? ----
    // Read from the LUT Bio-Formats applies to each channel from the CZI
    // metadata: grey = DIC, otherwise the dominant component gives the
    // colour. If this fails, fall back to the fixed numbers of section 4.
    uDIC = chDIC; uRED = chRED; uGREEN = chGREEN; uBLUE = chBLUE;
    if (autoDetectChannels) {
        detectChannels(projected);
        okDet = 1;
        if (dcRED == 0)   okDet = 0;
        if (dcGREEN == 0) okDet = 0;
        if (dcBLUE == 0)  okDet = 0;
        if (okDet == 1) {
            uDIC = dcDIC; uRED = dcRED; uGREEN = dcGREEN; uBLUE = dcBLUE;
            print("   channels detected: DIC=" + uDIC + " RED=" + uRED +
                  " GREEN=" + uGREEN + " BLUE=" + uBLUE);
        } else {
            print("   channel detection failed (LUTs not coloured)" +
                  " -> falling back to the fixed numbers in section 4");
        }
    }

    run("Split Channels");

    // Name from C1-/C2-/...: do NOT use getList("image.titles"), which also
    // lists windows left open by previously processed files.
    if (uDIC   > 0 && isOpen("C" + uDIC   + "-" + projected)) { selectWindow("C" + uDIC   + "-" + projected); rename(base + "_DIC");   }
    if (uRED   > 0 && isOpen("C" + uRED   + "-" + projected)) { selectWindow("C" + uRED   + "-" + projected); rename(base + "_RED");   }
    if (uGREEN > 0 && isOpen("C" + uGREEN + "-" + projected)) { selectWindow("C" + uGREEN + "-" + projected); rename(base + "_GREEN"); }
    if (uBLUE  > 0 && isOpen("C" + uBLUE  + "-" + projected)) { selectWindow("C" + uBLUE  + "-" + projected); rename(base + "_BLUE");  }

    blueName  = base + "_BLUE";
    greenName = base + "_GREEN";
    redName   = base + "_RED";
    dicName   = base + "_DIC";

    if (!isOpen(blueName) || !isOpen(greenName) || !isOpen(redName)) {
        print("   ERROR: missing channels -> image skipped");
        continue;
    }
    if (isOpen(dicName)) { selectWindow(dicName); close(); }

    // ---- Annotation sizes, scaled to the image ----
    selectWindow(blueName);
    imgW = getWidth();
    labelFS = maxOf(8, round(imgW / 100 * labelScale));
    recapFS = maxOf(8, round(imgW /  90 * recapScale));
    lineW   = maxOf(2, round(imgW / 600 * strokeScale));

    // ==================================================================
    //  NUCLEI
    // ==================================================================
    nucleusMaskTitle = "Nucleus_mask_" + base;

    selectWindow(blueName);
    run("Duplicate...", "title=[" + nucleusMaskTitle + "]");
    run("Enhance Contrast", "saturated=" + enhanceSaturated);
    run("Multiply...", "value=" + nucleiMultiply);
    run("Subtract Background...", "rolling=" + nucleiBgRolling);
    run("Gaussian Blur...", "sigma=" + nucleiBlurSigma);
    setAutoThreshold(nucleiThresholdMethod + " dark");
    run("Convert to Mask");
    run("8-bit");
    run("Fill Holes");
    if (nucleiWatershed) run("Watershed");

    roiManager("Reset");
    selectWindow(nucleusMaskTitle);
    run("Analyze Particles...",
        "size=" + nucleiMinArea + "-" + nucleiMaxArea +
        " circularity=" + nucleiCirc + " show=Nothing add");

    nucleiROIs  = roiManager("count");
    nucleiCount = nucleiROIs;
    doubletIdx  = newArray(0);

    // Clusters: an object clearly larger than the field median counts as
    // several overlapping cells.
    if (splitBigNuclei && nucleiROIs > 3) {
        // getValue reads the measurement directly, bypassing the Results
        // table. That table is shared with the rest of ImageJ, so another
        // macro or a renamed window could empty it and getResult then
        // failed with "No results found".
        nAreas = newArray(nucleiROIs);
        for (n = 0; n < nucleiROIs; n++) {
            selectWindow(blueName);
            roiManager("Select", n);
            nAreas[n] = getValue("Area");
        }
        selectWindow(blueName);
        run("Select None");

        sortedA = Array.copy(nAreas);
        Array.sort(sortedA);
        medA = sortedA[floor(nucleiROIs / 2)];

        totalCells = 0;
        for (n = 0; n < nucleiROIs; n++) {
            est = 1;
            if (medA > 0 && nAreas[n] >= nucleiDoubletFactor * medA) {
                est = round(nAreas[n] / medA);
                if (est < 2) est = 2;
                if (est > nucleiMaxSplit) est = nucleiMaxSplit;
                doubletIdx = Array.concat(doubletIdx, n);
            }
            totalCells = totalCells + est;
        }
        nucleiCount = totalCells;
        print("   nuclei " + nucleiROIs + " objects, median " + d2s(medA, 1) +
              " um2, " + doubletIdx.length + " clusters -> " + nucleiCount + " cells");
    } else {
        print("   nuclei = " + nucleiCount);
    }

    nucleiZip = tmpZip;
    if (nucleiROIs > 0) roiManager("Save", nucleiZip);

    selectWindow(nucleusMaskTitle);
    close();

    // ------------------------------------------------------------------
// ==================================================================
    //  GREEN PARASITES
    // ==================================================================
    parasiteMaskTitle = "Parasite_mask_" + base;

    selectWindow(greenName);
    run("Duplicate...", "title=[" + parasiteMaskTitle + "]");
    run("Enhance Contrast", "saturated=" + enhanceSaturated);
    run("Multiply...", "value=" + greenMultiply);
    run("Subtract Background...", "rolling=" + greenBgRolling);
    run("Gaussian Blur...", "sigma=" + greenBlurSigma);

    run("Select None");
    setAutoThreshold(greenThresholdMethod + " dark");
    getThreshold(grLo, grHi);
    grLoAuto = grLo;
    grLo = round(grLo * greenThresholdFactor);
    setThreshold(grLo, grHi);
    print("   green threshold " + grLo + " (auto " + grLoAuto + " x " + greenThresholdFactor + ")");
    run("Convert to Mask");
    run("8-bit");
    run("Fill Holes");
    if (greenWatershed) run("Watershed");

    roiManager("Reset");
    selectWindow(parasiteMaskTitle);
    run("Analyze Particles...",
        "size=" + greenMinArea + "-" + greenMaxArea +
        " circularity=" + greenCirc + " show=Nothing add");
    greenCount = roiManager("count");
    print("   green parasites = " + greenCount);

    selectWindow(parasiteMaskTitle);
    close();

    // ------------------------------------------------------------------
// ==================================================================
    //  RED / GREEN COLOCALISATION
    //  Green ROIs stay loaded: indices 0 to greenCount-1
    // ==================================================================
    redMaskTitle = "Red_mask_" + base;

    selectWindow(redName);
    run("Duplicate...", "title=[" + redMaskTitle + "]");
    run("Enhance Contrast", "saturated=" + enhanceSaturated);
    run("Multiply...", "value=" + redMultiply);
    run("Subtract Background...", "rolling=" + redBgRolling);
    run("Gaussian Blur...", "sigma=" + redBlurSigma);
    setAutoThreshold(redThresholdMethod + " dark");
    getThreshold(redLo, redHi);
    redLoAuto = redLo;
    redLo = round(redLo * redThresholdFactor);
    setThreshold(redLo, redHi);
    print("   red threshold " + redLo + " (auto " + redLoAuto + " x " + redThresholdFactor + ")");
    run("Convert to Mask");
    run("8-bit");
    run("Fill Holes");

    colocIndices = newArray(0);
    n20 = 0; n35 = 0; n55 = 0;
    for (r = 0; r < greenCount; r++) {
        selectWindow(redMaskTitle);
        roiManager("Select", r);
        // Mean of the red mask inside the green ROI: 255 = fully positive,
        // 0 = not at all. Read without the Results table.
        coveragePercent = (getValue("Mean") / 255) * 100;
        if (coveragePercent >= 20) n20++;
        if (coveragePercent >= 35) n35++;
        if (coveragePercent >= 55) n55++;
        if (coveragePercent >= minOverlapPercent)
            colocIndices = Array.concat(colocIndices, r);
    }
    colocCount = colocIndices.length;
    // If the three counters are close, the overlap threshold hardly
    // matters; if they diverge, many objects are borderline.
    print("   double positives = " + colocCount + " / " + greenCount +
          "   (>=20% : " + n20 + ", >=35% : " + n35 + ", >=55% : " + n55 + ")");
    selectWindow(redMaskTitle);
    run("Select None");

    pCells = nucleiCount;
    pGreen = greenCount;
    pRed   = colocCount;

    // ------------------------------------------------------------------
    //  6-slice stack: raw first, annotated last
    // ------------------------------------------------------------------
    greenIdx = makeSeq(greenCount);

    // 1-3 raw single channels
    u1 = rawPanel(blueName,  1.0,          base + "_u1_blue");
    u2 = rawPanel(greenName, 1.0,          base + "_u2_green");
    u3 = rawPanel(redName,   displayBoost, base + "_u3_red");

    // 4 raw green+red merge
    dR = base + "_dR";
    dG = base + "_dG";
    boostCopy(redName,   dR, displayBoost);
    boostCopy(greenName, dG, displayBoost);
    mU = base + "_mU";
    run("Merge Channels...", "c1=[" + dR + "] c2=[" + dG + "] create keep");
    rename(mU);
    run("RGB Color");
    u4 = base + "_u4_merge";
    flattenClean(u4);
    if (isOpen(mU)) { selectWindow(mU); close(); }
    if (isOpen(dR)) { selectWindow(dR); close(); }
    if (isOpen(dG)) { selectWindow(dG); close(); }

    // 5 annotated merge: green outlines plus red double positives
    eR = base + "_eR";
    eG = base + "_eG";
    boostCopy(redName,   eR, displayBoost);
    boostCopy(greenName, eG, displayBoost);
    mA = base + "_mA";
    run("Merge Channels...", "c1=[" + eR + "] c2=[" + eG + "] create keep");
    rename(mA);
    run("RGB Color");
    q5 = base + "_q5_merge_annotated";
    flattenClean(q5);
    drawROIs(greenIdx,     0, 255,   0, lineW);
    drawROIs(colocIndices, 255,   0, 0, lineW);
    labelROIs(greenIdx,     0, 255,   0, labelFS, 0);
    labelROIs(colocIndices, 255,   0, 0, labelFS, labelFS + 2);
    drawRecap("Green = " + greenCount + "   Green+Red = " + colocCount,
              255, 255, 255, recapFS);
    if (isOpen(mA)) { selectWindow(mA); close(); }
    if (isOpen(eR)) { selectWindow(eR); close(); }
    if (isOpen(eG)) { selectWindow(eG); close(); }
    if (isOpen(redMaskTitle)) { selectWindow(redMaskTitle); close(); }

    // 6 annotated blue: the nuclei ROIs have to be reloaded, the green
    // segmentation emptied the ROI manager.
    roiManager("Reset");
    nucleiIdx = newArray(0);
    if (nucleiROIs > 0 && File.exists(nucleiZip)) {
        roiManager("Open", nucleiZip);
        nucleiIdx = makeSeq(roiManager("count"));
    }
    q6 = rawPanel(blueName, 1.0, base + "_q6_blue_annotated");
    drawROIs(nucleiIdx,   0, 255, 255, lineW);
    drawROIs(doubletIdx, 255,   0, 255, lineW);
    labelROIs(nucleiIdx, 0, 255, 255, labelFS, 0);
    drawRecap("Cells = " + nucleiCount + "   (" + nucleiROIs + " objects, " +
              doubletIdx.length + " clusters)", 0, 255, 255, recapFS);

    stackTitle = base + "_control";
    run("Concatenate...",
        "title=[" + stackTitle + "]" +
        " image1=[" + u1 + "] image2=[" + u2 + "] image3=[" + u3 + "]" +
        " image4=[" + u4 + "] image5=[" + q5 + "] image6=[" + q6 + "]");
    selectWindow(stackTitle);
    setSlice(1); setMetadata("Label", "1 - blue, raw: count the cells");
    setSlice(2); setMetadata("Label", "2 - green, raw: count the parasites");
    setSlice(3); setMetadata("Label", "3 - red, raw");
    setSlice(4); setMetadata("Label", "4 - merge G+R, raw: count the double positives");
    setSlice(5); setMetadata("Label", "5 - merge G+R annotated by the code");
    setSlice(6); setMetadata("Label", "6 - blue annotated by the code");
    setSlice(1);

    leftovers = getList("image.titles");
    for (w = leftovers.length - 1; w >= 0; w--) {
        t2 = leftovers[w];
        if (endsWith(t2, "_control") || endsWith(t2, "_control.tif")) continue;
        if (indexOf(t2, base) >= 0 || startsWith(t2, "Composite") ||
            startsWith(t2, "MAX_")  || startsWith(t2, "C1-") ||
            startsWith(t2, "C2-")   || startsWith(t2, "C3-") ||
            startsWith(t2, "C4-")) {
            selectWindow(t2);
            close();
        }
    }

    selectWindow(stackTitle);
    saveAs("Tiff", outFolder + stackTitle + ".tif");
    close();

    if (isOpen(blueName))  { selectWindow(blueName);  close(); }
    if (isOpen(greenName)) { selectWindow(greenName); close(); }
    if (isOpen(redName))   { selectWindow(redName);   close(); }
    roiManager("Reset");
    run("Clear Results");
}

// Flattened RGB copy of one channel, no annotation.
function rawPanel(src, boost, outTitle) {
    tmp = outTitle + "_tmp";
    boostCopy(src, tmp, boost);
    flattenClean(outTitle);
    if (isOpen(tmp)) { selectWindow(tmp); close(); }
    return outTitle;
}


// ======================================================================
//  WORKBOOK
// ======================================================================
// Three groups of three columns: Auto, Manual, Error %. The Manual cells
// are yellow and left empty - you fill them in, and the error columns
// compute themselves. Grey cells fall outside the sampling: the target
// for that class was already met, nothing to count there.
// The bottom block compares TOTALS rather than averaging percentages: a
// field with 2 objects must not weigh as much as one with 40.
function writeControlBook(outPath) {

    nb = vImg.length;
    if (nb == 0) return;

    b1 = "<Borders>" +
        "<Border ss:Position=\"Bottom\" ss:LineStyle=\"Continuous\" ss:Weight=\"1\"/>" +
        "<Border ss:Position=\"Left\" ss:LineStyle=\"Continuous\" ss:Weight=\"1\"/>" +
        "<Border ss:Position=\"Right\" ss:LineStyle=\"Continuous\" ss:Weight=\"1\"/>" +
        "<Border ss:Position=\"Top\" ss:LineStyle=\"Continuous\" ss:Weight=\"1\"/></Borders>";
    b2 = "<Borders>" +
        "<Border ss:Position=\"Bottom\" ss:LineStyle=\"Continuous\" ss:Weight=\"2\"/>" +
        "<Border ss:Position=\"Left\" ss:LineStyle=\"Continuous\" ss:Weight=\"2\"/>" +
        "<Border ss:Position=\"Right\" ss:LineStyle=\"Continuous\" ss:Weight=\"2\"/>" +
        "<Border ss:Position=\"Top\" ss:LineStyle=\"Continuous\" ss:Weight=\"2\"/></Borders>";
    al = "<Alignment ss:Horizontal=\"Center\" ss:Vertical=\"Center\" ss:WrapText=\"1\"/>";

    x = "<?xml version=\"1.0\"?>\n<?mso-application progid=\"Excel.Sheet\"?>\n" +
        "<Workbook xmlns=\"urn:schemas-microsoft-com:office:spreadsheet\"" +
        " xmlns:ss=\"urn:schemas-microsoft-com:office:spreadsheet\">\n<Styles>\n" +
        "<Style ss:ID=\"h\">"  + al + b1 + "<Font ss:Bold=\"1\"/>" +
            "<Interior ss:Color=\"#E8E8E8\" ss:Pattern=\"Solid\"/></Style>\n" +
        "<Style ss:ID=\"hg\">" + al + b2 + "<Font ss:Bold=\"1\" ss:Size=\"12\"/>" +
            "<Interior ss:Color=\"#C6D9F0\" ss:Pattern=\"Solid\"/></Style>\n" +
        "<Style ss:ID=\"c\">"  + al + b1 + "</Style>\n" +
        "<Style ss:ID=\"a\">"  + al + b1 +
            "<Interior ss:Color=\"#F2F2F2\" ss:Pattern=\"Solid\"/></Style>\n" +
        "<Style ss:ID=\"in\">" + al + b2 + "<Font ss:Bold=\"1\" ss:Color=\"#1F4E79\"/>" +
            "<Interior ss:Color=\"#FFF2CC\" ss:Pattern=\"Solid\"/></Style>\n" +
        "<Style ss:ID=\"na\">" + al + b1 +
            "<Interior ss:Color=\"#D9D9D9\" ss:Pattern=\"Solid\"/></Style>\n" +
        "<Style ss:ID=\"p\">"  + al + b1 + "<NumberFormat ss:Format=\"0.00\"/></Style>\n" +
        "<Style ss:ID=\"t\">"  + al + b2 + "<Font ss:Bold=\"1\"/>" +
            "<Interior ss:Color=\"#DCE6F1\" ss:Pattern=\"Solid\"/></Style>\n" +
        "<Style ss:ID=\"tp\">" + al + b2 + "<Font ss:Bold=\"1\"/>" +
            "<Interior ss:Color=\"#DCE6F1\" ss:Pattern=\"Solid\"/>" +
            "<NumberFormat ss:Format=\"0.00\"/></Style>\n</Styles>\n" +
        "<Worksheet ss:Name=\"Counting control\">\n<Table>\n" +
        "<Column ss:Width=\"150\"/><Column ss:Width=\"130\"/><Column ss:Width=\"150\"/>" +
        "<Column ss:Width=\"62\"/><Column ss:Width=\"70\"/><Column ss:Width=\"70\"/>" +
        "<Column ss:Width=\"62\"/><Column ss:Width=\"70\"/><Column ss:Width=\"70\"/>" +
        "<Column ss:Width=\"62\"/><Column ss:Width=\"70\"/><Column ss:Width=\"70\"/>\n";

    // ---------- header, two rows ----------
    x = x + "<Row ss:Height=\"26\">" +
        "<Cell ss:StyleID=\"h\" ss:MergeDown=\"1\"><Data ss:Type=\"String\">Experiment</Data></Cell>" +
        "<Cell ss:StyleID=\"h\" ss:MergeDown=\"1\"><Data ss:Type=\"String\">Condition</Data></Cell>" +
        "<Cell ss:StyleID=\"h\" ss:MergeDown=\"1\"><Data ss:Type=\"String\">Image</Data></Cell>" +
        "<Cell ss:StyleID=\"hg\" ss:MergeAcross=\"2\"><Data ss:Type=\"String\">CELLS</Data></Cell>" +
        "<Cell ss:StyleID=\"hg\" ss:MergeAcross=\"2\"><Data ss:Type=\"String\">GREEN PARASITES</Data></Cell>" +
        "<Cell ss:StyleID=\"hg\" ss:MergeAcross=\"2\"><Data ss:Type=\"String\">DOUBLE POSITIVES</Data></Cell>" +
        "</Row>\n";
    x = x + "<Row ss:Height=\"24\">" +
        "<Cell ss:Index=\"4\" ss:StyleID=\"h\"><Data ss:Type=\"String\">Auto</Data></Cell>" +
        subHead() + subHeadFull() + subHeadFull() + "</Row>\n";

    firstRow = 3;
    rw = 2;

    // ---------- one row per sampled image ----------
    for (m = 0; m < manipList.length; m++) {

        nExp = 0;
        for (t = 0; t < nb; t++) if (vExp[t] == manipList[m]) nExp = nExp + 1;
        if (nExp == 0) continue;

        firstOfExp = 1;
        for (c = 0; c < condsAll.length; c++) {

            nCd = 0;
            for (t = 0; t < nb; t++) {
                if (vExp[t] != manipList[m]) continue;
                if (vCond[t] != condsAll[c]) continue;
                nCd = nCd + 1;
            }
            if (nCd == 0) continue;

            firstOfCond = 1;
            for (t = 0; t < nb; t++) {
                if (vExp[t] != manipList[m]) continue;
                if (vCond[t] != condsAll[c]) continue;

                // Experiment merged over its rows, condition merged over
                // its images; the following rows start at the image column.
                x = x + "<Row>";
                if (firstOfExp == 1) {
                    x = x + "<Cell ss:StyleID=\"t\" ss:MergeDown=\"" + (nExp - 1) +
                        "\"><Data ss:Type=\"String\">" + manipList[m] + "</Data></Cell>";
                    x = x + "<Cell ss:StyleID=\"t\" ss:MergeDown=\"" + (nCd - 1) +
                        "\"><Data ss:Type=\"String\">" + condsAll[c] + "</Data></Cell>";
                    x = x + sCell(vImg[t]);
                    firstOfExp = 0;
                    firstOfCond = 0;
                } else if (firstOfCond == 1) {
                    x = x + "<Cell ss:Index=\"2\" ss:StyleID=\"t\" ss:MergeDown=\"" +
                        (nCd - 1) + "\"><Data ss:Type=\"String\">" + condsAll[c] +
                        "</Data></Cell>";
                    x = x + sCell(vImg[t]);
                    firstOfCond = 0;
                } else {
                    x = x + sCellIdx(vImg[t], 3);
                }
                x = x + grp(vCells[t], 4) + grp(vGreen[t], 7) + grp(vRed[t], 10);
                x = x + "</Row>\n";
                rw = rw + 1;
            }
        }
    }
    lastRow = rw;

    // ---------- totals ----------
    x = x + "<Row></Row>\n";
    x = x + "<Row ss:Height=\"28\">" + tCell("TOTAL") + tCell("") + tCell("") +
        totGrp(4, firstRow, lastRow) + totGrp(7, firstRow, lastRow) +
        totGrp(10, firstRow, lastRow) + "</Row>\n";

    x = x + "</Table>\n</Worksheet>\n</Workbook>\n";
    File.saveString(x, outPath);
}

function subHead() {
    return "<Cell ss:StyleID=\"h\"><Data ss:Type=\"String\">Manual</Data></Cell>" +
           "<Cell ss:StyleID=\"h\"><Data ss:Type=\"String\">Error %</Data></Cell>";
}
function subHeadFull() {
    return "<Cell ss:StyleID=\"h\"><Data ss:Type=\"String\">Auto</Data></Cell>" + subHead();
}

// Three cells for one class. col = column of the Auto cell.
// A negative auto value means the class was outside the sampling.
function grp(autoVal, col) {
    if (autoVal < 0) {
        return "<Cell ss:StyleID=\"na\"></Cell><Cell ss:StyleID=\"na\"></Cell>" +
               "<Cell ss:StyleID=\"na\"></Cell>";
    }
    cA = col; cM = col + 1;
    s = "<Cell ss:StyleID=\"a\"><Data ss:Type=\"Number\">" + autoVal + "</Data></Cell>";
    s = s + "<Cell ss:StyleID=\"in\"></Cell>";
    s = s + "<Cell ss:StyleID=\"p\" ss:Formula=\"=IF(OR(RC" + cM + "=&quot;&quot;,RC" + cM +
            "=0),&quot;&quot;,(RC" + cA + "-RC" + cM + ")/RC" + cM + "*100)\"></Cell>";
    return s;
}

// Totals for one class, and the bias computed on those totals.
function totGrp(col, r1, r2) {
    cA = col; cM = col + 1;
    s = "<Cell ss:StyleID=\"t\" ss:Formula=\"=SUM(R" + r1 + "C" + cA + ":R" + r2 + "C" + cA + ")\"></Cell>";
    s = s + "<Cell ss:StyleID=\"t\" ss:Formula=\"=SUM(R" + r1 + "C" + cM + ":R" + r2 + "C" + cM + ")\"></Cell>";
    s = s + "<Cell ss:StyleID=\"tp\" ss:Formula=\"=IF(RC" + cM + "=0,&quot;&quot;,(RC" + cA +
            "-RC" + cM + ")/RC" + cM + "*100)\"></Cell>";
    return s;
}

function tCell(s) { return "<Cell ss:StyleID=\"t\"><Data ss:Type=\"String\">" + s + "</Data></Cell>"; }
function sCell(s) { return "<Cell ss:StyleID=\"c\"><Data ss:Type=\"String\">" + s + "</Data></Cell>"; }
function sCellIdx(s, idx) {
    return "<Cell ss:Index=\"" + idx + "\" ss:StyleID=\"c\"><Data ss:Type=\"String\">" +
           s + "</Data></Cell>";
}


// ======================================================================
//  FUNCTIONS - image annotation
// ======================================================================

// Duplicate a channel and stretch its display range only.
// Pixel values are untouched, so measurements are unaffected.
function boostCopy(srcTitle, newTitle, boost) {
    selectWindow(srcTitle);
    run("Duplicate...", "title=[" + newTitle + "]");
    if (boost > 1) {
        getMinAndMax(mn, mx);
        setMinAndMax(mn, mn + (mx - mn) / boost);
    }
}

function makeSeq(nn) {
    a = newArray(nn);
    for (q = 0; q < nn; q++) a[q] = q;
    return a;
}

// Flatten to RGB with no overlay and no "Show All": annotations are then
// written into the pixels, which is reliable in batch mode.
function flattenClean(newTitle) {
    run("Select None");
    run("Remove Overlay");
    roiManager("Show None");
    run("Flatten");
    rename(newTitle);
}

function drawROIs(indices, colR, colG, colB, lineWidth) {
    if (indices.length == 0) return;
    run("Line Width...", "line=" + lineWidth);
    setForegroundColor(colR, colG, colB);
    for (q = 0; q < indices.length; q++) {
        roiManager("Select", indices[q]);
        run("Draw", "slice");
    }
    run("Select None");
}

// Numbers 1..N next to each ROI. The black halo keeps them readable on any
// background; without it, red vanishes on the red channel.
function labelROIs(indices, colR, colG, colB, fontSize, dy) {
    nLab = indices.length;
    if (nLab == 0) return;

    imgWl = getWidth();
    imgHl = getHeight();
    dx = maxOf(4, round(fontSize / 4));
    lx = newArray(nLab);
    ly = newArray(nLab);

    for (q = 0; q < nLab; q++) {
        roiManager("Select", indices[q]);
        getSelectionBounds(bx, by, bw, bh);
        px = bx + bw + dx;
        py = by + bh / 2 + dy;
        if (px > imgWl - fontSize * 1.6) px = bx - round(fontSize * 1.3);
        px = minOf(maxOf(px, 2), imgWl - round(fontSize * 1.5));
        py = minOf(maxOf(py, fontSize + 2), imgHl - 2);
        lx[q] = px;
        ly[q] = py;
    }

    run("Select None");
    setFont("SansSerif", fontSize, "bold");   // order: name, SIZE, style
    setJustification("left");
    halo = maxOf(1, round(fontSize / 12));

    setColor(0, 0, 0);
    for (q = 0; q < nLab; q++) {
        txt = "" + (q + 1);
        drawString(txt, lx[q] - halo, ly[q]);
        drawString(txt, lx[q] + halo, ly[q]);
        drawString(txt, lx[q], ly[q] - halo);
        drawString(txt, lx[q], ly[q] + halo);
    }
    setColor(colR, colG, colB);
    for (q = 0; q < nLab; q++) drawString("" + (q + 1), lx[q], ly[q]);
}

function drawRecap(txt, colR, colG, colB, fontSize) {
    run("Select None");
    setFont("SansSerif", fontSize, "bold");
    setJustification("left");
    halo = maxOf(1, round(fontSize / 14));
    px = 20;
    py = getHeight() - round(fontSize * 0.7);

    setColor(0, 0, 0);
    drawString(txt, px - halo, py);
    drawString(txt, px + halo, py);
    drawString(txt, px, py - halo);
    drawString(txt, px, py + halo);

    setColor(colR, colG, colB);
    drawString(txt, px, py);
}
