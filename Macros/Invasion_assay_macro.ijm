// ======================================================================
//  Toxoplasma invasion assay - automated quantification from .czi images
// ======================================================================
//
//  Folder layout (detected automatically):
//     <folder>/<experiment>/<condition>/*.czi     several experiments
//     <folder>/<condition>/*.czi                  a single experiment
//
//  For every image, a 5-slice QC stack:
//     1  BLUE        nuclei, cyan outlines + numbers (clusters in magenta)
//     2  GREEN       parasites, green outlines + numbers
//     3  RED         double positives only, red outlines + numbers
//     4  MERGE R+G   double positives
//     5  MERGE FULL  the three channels overlaid
//
//  Output, in <selected folder>/Results/:
//     Stacks_<name>/                QC stacks
//     Data/                         nucleus ROIs + one CSV per condition
//     Raw data.csv                  all experiments, raw counts
//     Final results <name>.xls      summary per experiment (top) and
//                                   per-image detail (bottom), live formulas
//     Histogram <name>.png          one panel per experiment
//     Histogram pooled <name>.png   all experiments pooled
//
//  Reported measure: green-only parasites per 100 host cells, i.e.
//  (green - double positive) / nuclei x 100.
//
//  Requires Fiji/ImageJ with Bio-Formats. Calibrated images assumed:
//  all area settings below are in um^2.
// ======================================================================

// ----------------------------------------------------------------------
// Globals shared with the functions. In the ImageJ macro language a
// function only sees its parameters and variables declared with "var".
// ----------------------------------------------------------------------
var summaryManip = newArray(0);
var summaryCond  = newArray(0);
var summaryNames = newArray(0);
var summaryCells = newArray(0);
var summaryGreen = newArray(0);
var summaryColoc = newArray(0);

var manipList = newArray(0);   // experiments, alphabetical
var condsAll  = newArray(0);   // conditions, final display order
var aggCells  = newArray(0);   // per experiment x condition, index m*nC+c
var aggGreen  = newArray(0);
var aggColoc  = newArray(0);
var aggN      = newArray(0);

var gOf    = newArray(0);      // group index of each condition
var gCount = 0;

var ordStart = newArray(0);    // image indices grouped by experiment x condition
var ordIdx   = newArray(0);

var mAcc = newArray(0);        // recursive scan
var cAcc = newArray(0);
var pAcc = newArray(0);

var dcDIC = 0;                 // detected channels
var dcRED = 0;
var dcGREEN = 0;
var dcBLUE = 0;

var rootName = "";
var singleManip = 0;           // 1 when no experiment level was found
var exclList = newArray(0);    // excluded conditions, see section 2


// ======================================================================
//  1. INPUT
// ======================================================================
// "ASK"    prompt at start
// "AUTO"   one folder; depth detected (experiment/condition or condition)
// "FOLDER" a single folder of .czi files
// "FILE"   a single .czi file
inputMode = "ASK";


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

// --- Conditions to leave out of every calculation ------------------------
// One entry per excluded condition, written "<experiment folder>/<condition
// folder>", or just "<condition folder>" when there is no experiment level.
// Excluded images are still processed and still produce QC stacks; they are
// marked in the workbook and shown as a red cross in the histograms.
// Example: newArray("Invasion PCIP - 11-08-26/10 uM PCIP")
excludedConditions = newArray();

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


// ======================================================================
//  PATHS
// ======================================================================
sep = File.separator;

if (inputMode == "ASK") {
    Dialog.create("Image source");
    Dialog.addChoice("Process:", newArray(
        "A folder (experiments and/or conditions in subfolders)",
        "A single folder of CZI files",
        "A single CZI file"),
        "A folder (experiments and/or conditions in subfolders)");
    Dialog.show();
    ch = Dialog.getChoice();
    if (startsWith(ch, "A single file"))       inputMode = "FILE";
    else if (startsWith(ch, "A single folder"))  inputMode = "FOLDER";
    else                                          inputMode = "AUTO";
}

manList = newArray(0);
condList = newArray(0);
pathList = newArray(0);

if (inputMode == "FILE") {
    p = File.openDialog("Select a .czi file");
    rootFolder = File.getParent(p) + sep;
    cnFile = cleanName(rootFolder);
    manList  = Array.concat(manList,  "-");
    condList = Array.concat(condList, cnFile);
    pathList = Array.concat(pathList, p);
}
else if (inputMode == "FOLDER") {
    rootFolder = getDirectory("Select the folder containing the CZI files");
    cond = cleanName(rootFolder);
    fl = getFileList(rootFolder);
    Array.sort(fl);
    for (a = 0; a < fl.length; a++) {
        if (endsWith(toLowerCase(fl[a]), ".czi")) {
            manList  = Array.concat(manList,  "-");
            condList = Array.concat(condList, cond);
            pathList = Array.concat(pathList, rootFolder + fl[a]);
        }
    }
}
else {
    rootFolder = getDirectory("Select the folder containing the experiments");
    // Recursive scan: the last two folders above a file give the
    // experiment and the condition. Depth-tolerant.
    scanDir(rootFolder, "", "");
    manList  = mAcc;
    condList = cAcc;
    pathList = pAcc;
}

print("\\Clear");

if (pathList.length == 0) {
    print("Contents of " + rootFolder + " :");
    dbg = getFileList(rootFolder);
    for (a = 0; a < dbg.length; a++) print("   " + dbg[a]);
    exit("No .czi file found under " + rootFolder +
         "\nSee the Log window for the folder contents.");
}

rootName   = cleanName(rootFolder);
outputRoot = rootFolder + "Results" + sep;
stacksRoot = outputRoot + "Stacks_" + rootName + sep;
dataRoot   = outputRoot + "Data" + sep;
File.makeDirectory(outputRoot);
File.makeDirectory(stacksRoot);
File.makeDirectory(dataRoot);

exclList = excludedConditions;
print(pathList.length + " image(s) to process, output in " + outputRoot);
for (a = 0; a < exclList.length; a++) print("excluded from calculations: " + exclList[a]);


// ======================================================================
//  PROCESSING
// ======================================================================
setBatchMode(true);

for (i = 0; i < pathList.length; i++) {

    path  = pathList[i];
    manip = manList[i];
    cond  = condList[i];
    base  = File.getName(path);
    base  = substring(base, 0, lengthOf(base) - 4);

    if (manip == "-") {
        stackFolder = stacksRoot + cond + sep;
        dataFolder  = dataRoot + cond + sep;
    } else {
        File.makeDirectory(stacksRoot + manip + sep);
        File.makeDirectory(dataRoot + manip + sep);
        stackFolder = stacksRoot + manip + sep + cond + sep;
        dataFolder  = dataRoot + manip + sep + cond + sep;
    }
    File.makeDirectory(stackFolder);
    File.makeDirectory(dataFolder);

    print("[" + (i+1) + "/" + pathList.length + "] " + manip + " / " + cond + " / " + base);

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

    nucleiZip = dataFolder + base + "_nuclei_ROIs.zip";
    if (nucleiROIs > 0) roiManager("Save", nucleiZip);

    selectWindow(nucleusMaskTitle);
    close();

    // ---- Panel BLUE ----
    nucleiIdx = makeSeq(nucleiROIs);
    blueDup = "BLUE_dup_" + base;
    selectWindow(blueName);
    run("Duplicate...", "title=[" + blueDup + "]");
    blueQC = base + "_BLUE_QC";
    flattenClean(blueQC);

    drawROIs(nucleiIdx,   0, 255, 255, lineW);
    drawROIs(doubletIdx, 255,   0, 255, lineW);   // amas : a verifier a l'oeil
    labelROIs(nucleiIdx, 0, 255, 255, labelFS, 0);
    drawRecap("Number of cells = " + nucleiCount + "   (" + nucleiROIs +
              " objets, " + doubletIdx.length + " amas)", 0, 255, 255, recapFS);

    if (isOpen(blueDup)) { selectWindow(blueDup); close(); }

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

    // ---- Panel GREEN ----
    greenIdx = makeSeq(greenCount);
    greenDup = "GREEN_dup_" + base;
    selectWindow(greenName);
    run("Duplicate...", "title=[" + greenDup + "]");
    greenQC = base + "_GREEN_QC";
    flattenClean(greenQC);

    drawROIs(greenIdx, 0, 255, 0, lineW);
    labelROIs(greenIdx, 0, 255, 0, labelFS, 0);
    drawRecap("Green parasites = " + greenCount, 0, 255, 0, recapFS);

    if (isOpen(greenDup)) { selectWindow(greenDup); close(); }

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

    // ---- Panel RED ----
    redDup = "RED_dup_" + base;
    boostCopy(redName, redDup, displayBoost);
    redQC = base + "_RED_QC";
    flattenClean(redQC);

    drawROIs(colocIndices, 255, 0, 0, lineW);
    labelROIs(colocIndices, 255, 0, 0, labelFS, 0);
    drawRecap("Double positive (red+green) = " + colocCount, 255, 0, 0, recapFS);

    if (isOpen(redDup)) { selectWindow(redDup); close(); }

    // ---- Brightened display copies (measurements unchanged) ----
    redDisp   = base + "_RED_disp";
    greenDisp = base + "_GREEN_disp";
    boostCopy(redName,   redDisp,   displayBoost);
    boostCopy(greenName, greenDisp, displayBoost);

    // ---- Panel MERGE R+G ----
    mergeRG = base + "_Merge_RG";
    run("Merge Channels...",
        "c1=[" + redDisp + "] c2=[" + greenDisp + "] create keep");
    rename(mergeRG);
    run("RGB Color");
    mergeRGQC = base + "_Merge_RG_QC";
    flattenClean(mergeRGQC);

    drawROIs(colocIndices, 255, 0, 0, lineW);
    labelROIs(colocIndices, 255, 0, 0, labelFS, 0);
    drawRecap("Green+Red double positive = " + colocCount + " / " + greenCount,
              255, 255, 255, recapFS);

    if (isOpen(mergeRG))      { selectWindow(mergeRG);      close(); }
    if (isOpen(redMaskTitle)) { selectWindow(redMaskTitle); close(); }

    // ---- Panel MERGE COMPLET ----
    nucleiOffset = greenCount;
    if (nucleiROIs > 0 && File.exists(nucleiZip)) roiManager("Open", nucleiZip);
    nucleiInMerge = roiManager("count") - nucleiOffset;
    if (nucleiInMerge < 0) nucleiInMerge = 0;
    nucleiMergeIdx = newArray(nucleiInMerge);
    for (n = 0; n < nucleiInMerge; n++) nucleiMergeIdx[n] = nucleiOffset + n;

    fullMerge = base + "_Merge_Full";
    run("Merge Channels...",
        "c1=[" + redDisp + "] c2=[" + greenDisp + "] c3=[" + blueName + "] create keep");
    rename(fullMerge);
    run("RGB Color");
    fullMergeQC = base + "_Merge_Full_QC";
    flattenClean(fullMergeQC);

    drawROIs(nucleiMergeIdx, 0, 255, 255, lineW);
    drawROIs(greenIdx,       0, 255,   0, lineW);
    drawROIs(colocIndices, 255,   0,   0, lineW);   // rouge par-dessus le vert

    labelROIs(nucleiMergeIdx, 0, 255, 255, labelFS, 0);
    labelROIs(greenIdx,       0, 255,   0, labelFS, 0);
    labelROIs(colocIndices, 255,   0,   0, labelFS, labelFS + 2);

    drawRecap("Cells = " + nucleiCount + "   Green = " + greenCount +
              "   Red+Green = " + colocCount, 255, 255, 255, recapFS);

    if (isOpen(fullMerge)) { selectWindow(fullMerge); close(); }

    // ---- Assemblage de la pile et sauvegarde ----
    stackTitle = base + "_QC_Stack";
    run("Concatenate...",
        "title=[" + stackTitle + "]" +
        " image1=[" + blueQC + "]" +
        " image2=[" + greenQC + "]" +
        " image3=[" + redQC + "]" +
        " image4=[" + mergeRGQC + "]" +
        " image5=[" + fullMergeQC + "]");

    selectWindow(stackTitle);
    setSlice(1); setMetadata("Label", "BLUE - nuclei");
    setSlice(2); setMetadata("Label", "GREEN - parasites");
    setSlice(3); setMetadata("Label", "RED - double positive only");
    setSlice(4); setMetadata("Label", "MERGE RED+GREEN");
    setSlice(5); setMetadata("Label", "MERGE FULL");
    setSlice(1);

    // Close everything left for this image except the QC stack.
    // "RGB Color" on a composite leaves a second window with the same title.
    leftovers = getList("image.titles");
    for (w = leftovers.length - 1; w >= 0; w--) {
        t = leftovers[w];
        if (endsWith(t, "_QC_Stack") || endsWith(t, "_QC_Stack.tif")) continue;
        if (indexOf(t, base) >= 0 || startsWith(t, "Composite") ||
            startsWith(t, "MAX_")  || startsWith(t, "C1-") ||
            startsWith(t, "C2-")   || startsWith(t, "C3-") ||
            startsWith(t, "C4-")) {
            selectWindow(t);
            close();
        }
    }

    selectWindow(stackTitle);
    saveAs("Tiff", stackFolder + stackTitle + ".tif");
    if (keepWindowsOpen) {
        setBatchMode("show");
        setBatchMode(true);
    } else {
        close();
    }

    if (isOpen(redDisp))   { selectWindow(redDisp);   close(); }
    if (isOpen(greenDisp)) { selectWindow(greenDisp); close(); }
    if (isOpen(blueName))  { selectWindow(blueName);  close(); }
    if (isOpen(greenName)) { selectWindow(greenName); close(); }
    if (isOpen(redName))   { selectWindow(redName);   close(); }
    roiManager("Reset");
    run("Clear Results");

    summaryManip = Array.concat(summaryManip, manip);
    summaryCond  = Array.concat(summaryCond,  cond);
    summaryNames = Array.concat(summaryNames, base);
    summaryCells = Array.concat(summaryCells, nucleiCount);
    summaryGreen = Array.concat(summaryGreen, greenCount);
    summaryColoc = Array.concat(summaryColoc, colocCount);
}

setBatchMode(false);


// ======================================================================
//  RESULT ASSEMBLY
// ======================================================================
nRows = summaryNames.length;

// ---- experiments, alphabetical ----
for (t = 0; t < nRows; t++) {
    known = 0;
    for (u = 0; u < manipList.length; u++) if (manipList[u] == summaryManip[t]) known = 1;
    if (known == 0) manipList = Array.concat(manipList, summaryManip[t]);
}
Array.sort(manipList);
nM = manipList.length;
singleManip = 0;
if (nM == 1) { if (manipList[0] == "-") singleManip = 1; }

// ---- union of conditions ----
for (t = 0; t < nRows; t++) {
    known = 0;
    for (u = 0; u < condsAll.length; u++) if (condsAll[u] == summaryCond[t]) known = 1;
    if (known == 0) condsAll = Array.concat(condsAll, summaryCond[t]);
}
nC = condsAll.length;

// ---- per experiment x condition totals ----
aggCells = newArray(nM * nC);
aggGreen = newArray(nM * nC);
aggColoc = newArray(nM * nC);
aggN     = newArray(nM * nC);
for (t = 0; t < nRows; t++) {
    mi = -1;
    for (u = 0; u < nM; u++) if (manipList[u] == summaryManip[t]) mi = u;
    ci = -1;
    for (u = 0; u < nC; u++) if (condsAll[u] == summaryCond[t]) ci = u;
    if (mi < 0) continue;
    if (ci < 0) continue;
    if (isExcluded(summaryManip[t], summaryCond[t]) == 1) continue;
    k = mi * nC + ci;
    aggCells[k] = aggCells[k] + summaryCells[t];
    aggGreen[k] = aggGreen[k] + summaryGreen[t];
    aggColoc[k] = aggColoc[k] + summaryColoc[t];
    aggN[k]     = aggN[k] + 1;
}

// ---- display order of the conditions ----
// Controls first, then increasing concentration. Whether a condition is
// present in every experiment does NOT enter the sort: the detail tables
// are aligned with grey cells, which must not disturb the dose order.
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
        for (m = 0; m < nM; m++) {
            i1 = m * nC + a;
            i2 = m * nC + mn;
            v1 = aggCells[i1]; aggCells[i1] = aggCells[i2]; aggCells[i2] = v1;
            v2 = aggGreen[i1]; aggGreen[i1] = aggGreen[i2]; aggGreen[i2] = v2;
            v3 = aggColoc[i1]; aggColoc[i1] = aggColoc[i2]; aggColoc[i2] = v3;
            v4 = aggN[i1];     aggN[i1]     = aggN[i2];     aggN[i2]     = v4;
        }
    }
}

// ---- condition groups ----
// Same concentration = same group (controls form a group of their own).
// Once sorted, group members are contiguous, so a pooled row can sum them.
gOf = newArray(nC);
gCount = 0;
for (c = 0; c < nC; c++) {
    if (c == 0) {
        gOf[0] = 0;
        gCount = 1;
    } else {
        if (sameGroup(condsAll[c - 1], condsAll[c]) == 1) {
            gOf[c] = gOf[c - 1];
        } else {
            gOf[c] = gCount;
            gCount = gCount + 1;
        }
    }
}


// ---- image index, grouped by experiment then condition ----
// Avoids rescanning the whole table for every cell written.
ordStart = newArray(nM * nC);
ordIdx   = newArray(nRows);
pos = 0;
for (m = 0; m < nM; m++) {
    for (c = 0; c < nC; c++) {
        ordStart[m * nC + c] = pos;
        for (t = 0; t < nRows; t++) {
            if (summaryManip[t] != manipList[m]) continue;
            if (summaryCond[t]  != condsAll[c])  continue;
            ordIdx[pos] = t;
            pos = pos + 1;
        }
    }
}


// ======================================================================
//  OUTPUT FILES
// ======================================================================
csv = "Experiment,Condition,Image,Cells,Green,Green_and_Red,Green_only,Green_only_per_100_cells\n";
for (m = 0; m < nM; m++) {
    for (c = 0; c < nC; c++) {
        for (t = 0; t < nRows; t++) {
            if (summaryManip[t] != manipList[m]) continue;
            if (summaryCond[t]  != condsAll[c])  continue;
            seuls = summaryGreen[t] - summaryColoc[t];
            if (summaryCells[t] > 0) pv = seuls * 100.0 / summaryCells[t]; else pv = 0;
            csv = csv + summaryManip[t] + "," + summaryCond[t] + "," + summaryNames[t] +
                  "," + summaryCells[t] + "," + summaryGreen[t] + "," + summaryColoc[t] +
                  "," + seuls + "," + num2(pv) + "\n";
        }
    }
}
File.saveString(csv, outputRoot + "Raw data.csv");

// one CSV per condition, next to its ROIs
for (m = 0; m < nM; m++) {
    for (c = 0; c < nC; c++) {
        if (aggN[m * nC + c] == 0) continue;
        sub = "Image,Cells,Green,Green_and_Red,Green_only,Green_only_per_100_cells\n";
        for (t = 0; t < nRows; t++) {
            if (summaryManip[t] != manipList[m]) continue;
            if (summaryCond[t]  != condsAll[c])  continue;
            seuls = summaryGreen[t] - summaryColoc[t];
            if (summaryCells[t] > 0) pv = seuls * 100.0 / summaryCells[t]; else pv = 0;
            sub = sub + summaryNames[t] + "," + summaryCells[t] + "," + summaryGreen[t] +
                  "," + summaryColoc[t] + "," + seuls + "," + num2(pv) + "\n";
        }
        if (singleManip == 1) subDir = dataRoot + condsAll[c] + sep;
        else                  subDir = dataRoot + manipList[m] + sep + condsAll[c] + sep;
        File.saveString(sub, subDir + "Raw data - " + condsAll[c] + ".csv");
    }
}

writeWorkbook(outputRoot + "Final results " + rootName + ".xls");

if (nC > 0) drawHistogram(outputRoot + "Histogramme " + rootName + ".png");
if (nC > 0 && singleManip == 0 && nM > 1)
    drawGlobalHistogram(outputRoot + "Histogramme global " + rootName + ".png");

print("Done. " + nRows + " image(s), " + nM + " experiment(s), " + nC + " condition(s).");
if (exclList.length > 0)
    print(exclList.length + " condition(s) excluded from the calculations.");
print("Output: " + outputRoot);


// ======================================================================
//  FUNCTIONS - assembly and output
// ======================================================================

// Recursive scan: the last two folders above a file give the experiment
// and the condition. With a single level, the experiment is "-" and that
// level is the condition.
function scanDir(dir, lvl1, lvl2) {
    lst = getFileList(dir);
    Array.sort(lst);
    for (q = 0; q < lst.length; q++) {
        nm = lst[q];
        full = dir + nm;
        if (File.isDirectory(full)) {
            cn = cleanName(nm);
            if (cn == "Results") continue;
            if (cn == "Validation") continue;
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

// Is this experiment x condition excluded from the calculations?
function isExcluded(manip, cond) {
    for (q = 0; q < exclList.length; q++) {
        if (exclList[q] == cond) return 1;
        if (exclList[q] == manip + "/" + cond) return 1;
    }
    return 0;
}

// Do two conditions belong to the same group?
// Two controls: yes. Otherwise, same leading number in the name: yes.
function sameGroup(a, b) {
    ca = isControlName(a);
    cb = isControlName(b);
    if (ca == 1 && cb == 1) return 1;
    if (ca != cb) return 0;
    na = firstNumber(a);
    nb = firstNumber(b);
    if (na < 0) return 0;
    if (nb < 0) return 0;
    if (na == nb) return 1;
    return 0;
}

// Group label: longest common prefix of its members, which gives
// "Control" for Control 1 / Control 2, and "0.25 uM PCIP" for
// 0.25 uM PCIP / 0.25 uM PCIP_slide 1 / _slide 2.
function groupLabel(g) {
    lab = "";
    started = 0;
    for (c = 0; c < condsAll.length; c++) {
        if (gOf[c] != g) continue;
        if (started == 0) { lab = condsAll[c]; started = 1; }
        else              lab = commonPrefix(lab, condsAll[c]);
    }
    lab = trimEnd(lab);
    if (lengthOf(lab) < 2) {
        for (c = 0; c < condsAll.length; c++) {
            if (gOf[c] == g) return condsAll[c];
        }
    }
    return lab;
}

function commonPrefix(a, b) {
    n = minOf(lengthOf(a), lengthOf(b));
    k = 0;
    stop = 0;
    for (q = 0; q < n; q++) {
        if (stop == 1) continue;
        if (substring(a, q, q + 1) == substring(b, q, q + 1)) k = q + 1;
        else stop = 1;
    }
    return substring(a, 0, k);
}

function trimEnd(s) {
    while (lengthOf(s) > 0) {
        c2 = substring(s, lengthOf(s) - 1, lengthOf(s));
        if (c2 == " " || c2 == "_" || c2 == "-") s = substring(s, 0, lengthOf(s) - 1);
        else return s;
    }
    return s;
}

// Sort key: controls first, then the leading number (0.05, 0.1, 1, 10)
function condSortKey(name) {
    num = firstNumber(name);
    if (isControlName(name) == 1) {
        if (num < 0) return -1000000;
        return -1000000 + num;
    }
    if (num < 0) return 1000000;
    return num;
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

// Number as text, 2 decimals, DOT separator. d2s follows the system
// locale and would emit a comma, breaking both the CSV and the XML.
function num2(v) {
    return replace(d2s(v, 2), ",", ".");
}

function cleanName(s) {
    while (endsWith(s, "/") || endsWith(s, "\\"))
        s = substring(s, 0, lengthOf(s) - 1);
    return File.getName(s);
}


// ======================================================================
//  WORKBOOK
// ======================================================================
// Mise en page :
//   top     summary, one stacked block per experiment
//   bottom  per-image detail, experiments side by side
//
// The detail is built first into a separate string, which yields the row
// numbers of the TOTAL lines; the summary then points at them by formula,
// so every value stays live if a count is edited.
function writeWorkbook(outPath) {

    nMw = manipList.length;
    nCw = condsAll.length;
    BW  = 8;          // 7 colonnes par manip + 1 colonne de separation

    // ---------- summary geometry ----------
    // Stacked blocks, one experiment after another, merged left column.
    // Conditions are walked group by group (same concentration). As soon as
    // a group holds two or more conditions, a pooled row sums them just
    // below. For the controls that pooled row is the 100 % reference -
    // otherwise the reference control would be picked arbitrarily.
    sRowOf   = newArray(nMw * nCw);
    grpRowOf = newArray(nMw * gCount);
    gFirst   = newArray(nMw * gCount);
    gLast    = newArray(nMw * gCount);
    ctrlS    = newArray(nMw);
    nPresM   = newArray(nMw);

    sr = 1;                                   // row 1 = headers
    for (m = 0; m < nMw; m++) {
        tot = 0;
        for (g = 0; g < gCount; g++) {
            k2 = 0;
            isCtrlG = 0;
            for (c = 0; c < nCw; c++) {
                if (gOf[c] != g) continue;
                if (aggN[m * nCw + c] == 0) continue;
                if (isControlName(condsAll[c]) == 1) isCtrlG = 1;
                sr = sr + 1;
                sRowOf[m * nCw + c] = sr;
                if (k2 == 0) gFirst[m * gCount + g] = sr;
                gLast[m * gCount + g] = sr;
                k2 = k2 + 1;
            }
            if (k2 == 0) continue;
            tot = tot + k2;
            if (k2 > 1) {
                sr = sr + 1;
                grpRowOf[m * gCount + g] = sr;
                tot = tot + 1;
            } else {
                grpRowOf[m * gCount + g] = 0;
            }
            if (isCtrlG == 1) {
                if (k2 > 1) ctrlS[m] = grpRowOf[m * gCount + g];
                else        ctrlS[m] = gLast[m * gCount + g];
            }
        }
        nPresM[m] = tot;
        sr = sr + 1;                          // separator between experiments
    }
    sr = sr + 1;                              // blank row
    summaryHeight = sr;

    // ---------- detail geometry ----------
    // Side-by-side blocks: experiment m uses columns m*8+1 to m*8+7. The
    // block of a given condition starts on the SAME row in every experiment
    // so they can be compared horizontally. If one experiment has fewer
    // images its extra rows stay empty; the TOTAL sum ignores them.
    maxK   = newArray(nCw);
    totRel = newArray(nCw);
    for (c = 0; c < nCw; c++) {
        mk = 0;
        for (m = 0; m < nMw; m++) {
            if (aggN[m * nCw + c] > mk) mk = aggN[m * nCw + c];
        }
        maxK[c] = mk;
    }

    rr = 0;
    if (singleManip == 0) rr = rr + 1;        // ligne des titres de manip
    rr = rr + 1;                              // ligne des en-tetes
    for (c = 0; c < nCw; c++) {
        if (maxK[c] == 0) continue;
        rr = rr + maxK[c];
        rr = rr + 1;
        totRel[c] = rr;                       // ligne TOTAL, relative
        rr = rr + 1;                          // separation
    }

    // ---------- write ----------
    nCols = nMw * BW;
    if (nCols < 8) nCols = 8;
    x = xmlHead(nCols);

    // ---------- SUMMARY ----------
    x = x + "<Row ss:Height=\"34\">" +
        hCell("Experiment") + hCell("Condition") + hCell("Cells") + hCell("Green") +
        hCell("Green + Red") + hCell("Green only") +
        hCell("Verts seuls / 100 cells") + hCell("% of control") + "</Row>\n";

    for (m = 0; m < nMw; m++) {

        nPres = nPresM[m];
        if (nPres == 0) {
            x = x + "<Row></Row>\n";
            continue;
        }

        labM = manipList[m];
        if (singleManip == 0) labM = "EXP " + (m + 1) + "&#10;" + manipList[m];

        first = 1;
        baseD = m * BW + 1;               // column of the experiment's DETAIL block

        for (g = 0; g < gCount; g++) {

            k2 = 0;
            for (c = 0; c < nCw; c++) {
                if (gOf[c] != g) continue;
                if (aggN[m * nCw + c] == 0) continue;
                tr = summaryHeight + totRel[c];
                rowStr = sumRowRef(first, labM, nPres - 1, condsAll[c],
                                   baseD, tr, ctrlS[m]);
                x = x + rowStr;
                first = 0;
                k2 = k2 + 1;
            }

            if (k2 > 1) {
                gl = groupLabel(g) + " (pooled)";
                rowStr = sumRowGroup(first, labM, nPres - 1, gl,
                                     gFirst[m * gCount + g],
                                     gLast[m * gCount + g], ctrlS[m]);
                x = x + rowStr;
                first = 0;
            }
        }

        x = x + "<Row ss:Height=\"10\"></Row>\n";
    }
    x = x + "<Row></Row>\n";

    // ---------- DETAIL, exp. cote a cote ----------
    if (singleManip == 0) {
        x = x + "<Row ss:Height=\"30\">";
        for (m = 0; m < nMw; m++) {
            x = x + "<Cell ss:Index=\"" + (m * BW + 1) + "\" ss:StyleID=\"hm\"" +
                " ss:MergeAcross=\"6\"><Data ss:Type=\"String\">" +
                "EXP " + (m + 1) + "  -  " + manipList[m] + "</Data></Cell>";
        }
        x = x + "</Row>\n";
    }

    x = x + "<Row ss:Height=\"34\">";
    for (m = 0; m < nMw; m++) {
        x = x + "<Cell ss:Index=\"" + (m * BW + 1) + "\" ss:StyleID=\"h\">" +
            "<Data ss:Type=\"String\">Condition</Data></Cell>" +
            hCell("Image") + hCell("Cells") + hCell("Green") +
            hCell("Green + Red") + hCell("Green only") +
            hCell("Verts seuls / 100 cells");
    }
    x = x + "</Row>\n";

    for (c = 0; c < nCw; c++) {
        if (maxK[c] == 0) continue;

        // --- image rows ---
        for (j = 0; j < maxK[c]; j++) {
            x = x + "<Row>";
            for (m = 0; m < nMw; m++) {
                baseD = m * BW + 1;
                kk = aggN[m * nCw + c];

                if (kk == 0) {
                    x = x + "<Cell ss:Index=\"" + baseD + "\" ss:StyleID=\"na\"></Cell>";
                    for (z = 0; z < 6; z++) x = x + "<Cell ss:StyleID=\"na\"></Cell>";
                } else if (j < kk) {
                    t = ordIdx[ordStart[m * nCw + c] + j];
                    if (j == 0) {
                        lab = condsAll[c];
                        if (singleManip == 0) lab = lab + "&#10;Manip " + (m + 1);
                        // MergeDown = maxK: the label spans the image rows
                        // AND the TOTAL row, even when this experiment has
                        // fewer images than the largest one.
                        x = x + "<Cell ss:Index=\"" + baseD + "\" ss:StyleID=\"t\"" +
                            " ss:MergeDown=\"" + maxK[c] + "\"><Data ss:Type=\"String\">" +
                            lab + "</Data></Cell>" + sCell(summaryNames[t]);
                    } else {
                        x = x + sCellIdx(summaryNames[t], baseD + 1);
                    }
                    x = x + nCell(summaryCells[t]) + nCell(summaryGreen[t]) +
                        nCell(summaryColoc[t]) +
                        fCell("c", "=RC[-2]-RC[-1]") +
                        fCell("p", "=IF(RC[-4]=0,0,RC[-1]/RC[-4]*100)");
                } else {
                    // this experiment has fewer images: padding rows
                    x = x + "<Cell ss:Index=\"" + (baseD + 1) + "\" ss:StyleID=\"c\"></Cell>";
                    for (z = 0; z < 5; z++) x = x + "<Cell ss:StyleID=\"c\"></Cell>";
                }
            }
            x = x + "</Row>\n";
        }

        // --- TOTAL row ---
        x = x + "<Row ss:Height=\"28\">";
        for (m = 0; m < nMw; m++) {
            baseD = m * BW + 1;
            kk = aggN[m * nCw + c];
            if (kk == 0) {
                x = x + "<Cell ss:Index=\"" + baseD + "\" ss:StyleID=\"na\"></Cell>";
                for (z = 0; z < 6; z++) x = x + "<Cell ss:StyleID=\"na\"></Cell>";
            } else {
                x = x + tCellIdx("TOTAL", baseD + 1) +
                    fCell("t", "=SUM(R[-" + maxK[c] + "]C:R[-1]C)") +
                    fCell("t", "=SUM(R[-" + maxK[c] + "]C:R[-1]C)") +
                    fCell("t", "=SUM(R[-" + maxK[c] + "]C:R[-1]C)") +
                    fCell("t", "=RC[-2]-RC[-1]") +
                    fCell("tp", "=IF(RC[-4]=0,0,RC[-1]/RC[-4]*100)");
            }
        }
        x = x + "</Row>\n";

        x = x + "<Row ss:Height=\"8\"></Row>\n";
    }

    x = x + "</Table>\n</Worksheet>\n</Workbook>\n";
    File.saveString(x, outPath);
}

// A summary row pointing at the matching TOTAL row of the detail.
function sumRowRef(isFirst, labM, mergeN, condLabel, baseD, tr, ctrlRow) {
    r = "<Row ss:Height=\"26\">";
    if (isFirst == 1) r = r + mCell(labM, mergeN) + tCell(condLabel);
    else              r = r + tCellIdx(condLabel, 2);
    r = r + fCell("t",  "=R" + tr + "C" + (baseD + 2)) +
            fCell("t",  "=R" + tr + "C" + (baseD + 3)) +
            fCell("t",  "=R" + tr + "C" + (baseD + 4)) +
            fCell("t",  "=RC[-2]-RC[-1]") +
            fCell("tp", "=IF(RC[-4]=0,0,RC[-1]/RC[-4]*100)");
    if (ctrlRow > 0)
        r = r + fCell("tp", "=IF(R" + ctrlRow + "C7=0,0,RC[-1]/R" + ctrlRow + "C7*100)");
    else
        r = r + "<Cell ss:StyleID=\"na\"></Cell>";
    return r + "</Row>\n";
}

// Pooled row: sum of the group member rows just above. Live formulas,
// like the rest of the table.
function sumRowGroup(isFirst, labM, mergeN, grpLab, r1, r2, ctrlRow) {
    r = "<Row ss:Height=\"26\">";
    if (isFirst == 1) r = r + mCell(labM, mergeN) + tCell(grpLab);
    else              r = r + tCellIdx(grpLab, 2);
    r = r + fCell("t",  "=SUM(R" + r1 + "C3:R" + r2 + "C3)") +
            fCell("t",  "=SUM(R" + r1 + "C4:R" + r2 + "C4)") +
            fCell("t",  "=SUM(R" + r1 + "C5:R" + r2 + "C5)") +
            fCell("t",  "=RC[-2]-RC[-1]") +
            fCell("tp", "=IF(RC[-4]=0,0,RC[-1]/RC[-4]*100)") +
            fCell("tp", "=IF(R" + ctrlRow + "C7=0,0,RC[-1]/R" + ctrlRow + "C7*100)");
    return r + "</Row>\n";
}

function xmlHead(nCols) {
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

    h = "<?xml version=\"1.0\"?>\n<?mso-application progid=\"Excel.Sheet\"?>\n" +
        "<Workbook xmlns=\"urn:schemas-microsoft-com:office:spreadsheet\"" +
        " xmlns:ss=\"urn:schemas-microsoft-com:office:spreadsheet\">\n<Styles>\n" +
        "<Style ss:ID=\"h\">"  + al + b1 + "<Font ss:Bold=\"1\"/>" +
            "<Interior ss:Color=\"#E8E8E8\" ss:Pattern=\"Solid\"/></Style>\n" +
        "<Style ss:ID=\"hm\">" + al + b2 + "<Font ss:Bold=\"1\" ss:Size=\"12\"/>" +
            "<Interior ss:Color=\"#C6D9F0\" ss:Pattern=\"Solid\"/></Style>\n" +
        "<Style ss:ID=\"c\">"  + al + b1 + "</Style>\n" +
        "<Style ss:ID=\"p\">"  + al + b1 + "<NumberFormat ss:Format=\"0.00\"/></Style>\n" +
        "<Style ss:ID=\"na\">" + al + b1 +
            "<Interior ss:Color=\"#D9D9D9\" ss:Pattern=\"Solid\"/></Style>\n" +
        "<Style ss:ID=\"t\">"  + al + b2 + "<Font ss:Bold=\"1\"/>" +
            "<Interior ss:Color=\"#DCE6F1\" ss:Pattern=\"Solid\"/></Style>\n" +
        "<Style ss:ID=\"tp\">" + al + b2 + "<Font ss:Bold=\"1\"/>" +
            "<Interior ss:Color=\"#DCE6F1\" ss:Pattern=\"Solid\"/>" +
            "<NumberFormat ss:Format=\"0.00\"/></Style>\n</Styles>\n" +
        "<Worksheet ss:Name=\"Resultats\">\n<Table>\n";

    // Column 8 of each block carries "% of control" in the summary and acts
    // as the separator in the detail, so it cannot be narrow.
    for (q = 0; q < nCols; q++) {
        w = 80;
        r8 = q % 8;
        if (r8 == 0) w = 135;      // Manip (summary) / Condition (detail)
        if (r8 == 1) w = 130;      // Condition (summary) / Image (detail)
        if (r8 == 5) w = 110;
        if (r8 == 6) w = 125;
        if (r8 == 7) w = 95;
        h = h + "<Column ss:Width=\"" + w + "\"/>";
    }
    return h + "\n";
}

function hCell(s) { return "<Cell ss:StyleID=\"h\"><Data ss:Type=\"String\">" + s + "</Data></Cell>"; }
function sCell(s) { return "<Cell ss:StyleID=\"c\"><Data ss:Type=\"String\">" + s + "</Data></Cell>"; }
function nCell(v) { return "<Cell ss:StyleID=\"c\"><Data ss:Type=\"Number\">" + v + "</Data></Cell>"; }
function tCell(s) { return "<Cell ss:StyleID=\"t\"><Data ss:Type=\"String\">" + s + "</Data></Cell>"; }
function mCell(s, mergeDown) {
    return "<Cell ss:StyleID=\"t\" ss:MergeDown=\"" + mergeDown +
           "\"><Data ss:Type=\"String\">" + s + "</Data></Cell>";
}
function sCellIdx(s, idx) {
    return "<Cell ss:Index=\"" + idx + "\" ss:StyleID=\"c\"><Data ss:Type=\"String\">" +
           s + "</Data></Cell>";
}
function tCellIdx(s, idx) {
    return "<Cell ss:Index=\"" + idx + "\" ss:StyleID=\"t\"><Data ss:Type=\"String\">" +
           s + "</Data></Cell>";
}
function fCell(st, f) {
    return "<Cell ss:StyleID=\"" + st + "\" ss:Formula=\"" + f +
           "\"><Data ss:Type=\"Number\"></Data></Cell>";
}


// ======================================================================
//  HISTOGRAMS
// ======================================================================
// Two separate images:
//   Histogram <name>.png         one panel per experiment, every condition
//                                shown separately
//   Histogram pooled <name>.png  all experiments together, conditions of
//                                equal concentration pooled
// Both use the same vertical scale, computed on the bars they contain.
// ======================================================================

// Records to plot: one per (panel, bar).
var pName  = newArray(0);
var pPct   = newArray(0);
var pSD    = newArray(0);
var pCells = newArray(0);
var pNf    = newArray(0);
var pNm    = newArray(0);      // number of experiments; 0 outside pooled panel
var pStart = newArray(0);
var pCount = newArray(0);
var pTitle = newArray(0);
var pGlob  = newArray(0);
var wl1 = "";
var wl2 = "";

function resetRecords(nPan) {
    pName  = newArray(0);
    pPct   = newArray(0);
    pSD    = newArray(0);
    pCells = newArray(0);
    pNf    = newArray(0);
    pNm    = newArray(0);
    pStart = newArray(nPan);
    pCount = newArray(nPan);
    pTitle = newArray(nPan);
    pGlob  = newArray(nPan);
}

// ---------------- one panel per experiment ----------------
function drawHistogram(outPath) {
    nMh  = manipList.length;
    nCh2 = condsAll.length;
    if (nMh == 0) return;
    if (nCh2 == 0) return;

    resetRecords(nMh);

    for (m = 0; m < nMh; m++) {
        pStart[m] = pName.length;
        pGlob[m]  = 0;
        if (singleManip == 0) pTitle[m] = "EXP " + (m + 1) + "  -  " + manipList[m];
        else                  pTitle[m] = manipList[m];
        cnt = 0;

        for (c = 0; c < nCh2; c++) {
            if (aggN[m * nCh2 + c] == 0) continue;

            kk = m * nCh2 + c;
            cc = aggCells[kk];
            gg = aggGreen[kk] - aggColoc[kk];
            nn = aggN[kk];
            pct = 0;
            if (cc > 0) pct = gg * 100.0 / cc;

            // field-to-field spread within this experiment
            vals = newArray(0);
            for (t = 0; t < summaryNames.length; t++) {
                if (summaryManip[t] != manipList[m]) continue;
                if (summaryCond[t] != condsAll[c]) continue;
                if (summaryCells[t] <= 0) continue;
                vv = (summaryGreen[t] - summaryColoc[t]) * 100.0 / summaryCells[t];
                vals = Array.concat(vals, vv);
            }
            sd = 0;
            if (vals.length > 1) {
                Array.getStatistics(vals, vmin, vmax, vmean, vstd);
                sd = vstd;
            }

            pName  = Array.concat(pName, condsAll[c]);
            pPct   = Array.concat(pPct, pct);
            pSD    = Array.concat(pSD, sd);
            pCells = Array.concat(pCells, cc);
            pNf    = Array.concat(pNf, nn);
            pNm    = Array.concat(pNm, 0);
            cnt = cnt + 1;
        }
        pCount[m] = cnt;
    }

    renderPanels(outPath, "Parasites intracellulaires / 100 cells",
                 "SD = ecart-type des valeurs champ par champ, au sein d'une meme manip");
}

// ---------------- single pooled panel ----------------
// Bar height: ratio computed on the totals of all experiments.
// Error bar: standard deviation BETWEEN EXPERIMENTS - the spread that
// matters for a figure, and the reason independent replicates are needed.
function drawGlobalHistogram(outPath) {
    nMh = manipList.length;
    if (nMh == 0) return;
    if (gCount == 0) return;

    resetRecords(1);
    pStart[0] = 0;
    pGlob[0]  = 1;
    pTitle[0] = "ALL EXPERIMENTS  (n = " + nMh + ")";
    cnt = 0;

    for (g = 0; g < gCount; g++) {

        cc = 0; gg = 0; nn = 0; nmc = 0;
        mv = newArray(0);

        for (mm = 0; mm < nMh; mm++) {
            xc = 0; xg = 0; xn = 0;
            for (c = 0; c < condsAll.length; c++) {
                if (gOf[c] != g) continue;
                kk = mm * condsAll.length + c;
                xc = xc + aggCells[kk];
                xg = xg + aggGreen[kk] - aggColoc[kk];
                xn = xn + aggN[kk];
            }
            if (xn == 0) continue;
            nmc = nmc + 1;
            cc = cc + xc;
            gg = gg + xg;
            nn = nn + xn;
            if (xc > 0) mv = Array.concat(mv, xg * 100.0 / xc);
        }
        if (nmc == 0) continue;

        pct = 0;
        if (cc > 0) pct = gg * 100.0 / cc;
        sd = 0;
        if (mv.length > 1) {
            Array.getStatistics(mv, vmin, vmax, vmean, vstd);
            sd = vstd;
        }

        gl = groupLabel(g);
        pName  = Array.concat(pName, gl);
        pPct   = Array.concat(pPct, pct);
        pSD    = Array.concat(pSD, sd);
        pCells = Array.concat(pCells, cc);
        pNf    = Array.concat(pNf, nn);
        pNm    = Array.concat(pNm, nmc);
        cnt = cnt + 1;
    }
    pCount[0] = cnt;

    renderPanels(outPath, "Parasites intracellulaires / 100 cells",
                 "SD = ecart-type ENTRE MANIPS.  Conditions de meme concentration regroupees.");
}

// ---------------- drawing ----------------
function renderPanels(outPath) {

    nPan = pStart.length;
    nBars = pName.length;
    if (nBars == 0) return;

    // vertical scale, integer step chosen from 1 / 2 / 5 x 10^n
    maxV = 0;
    for (u = 0; u < nBars; u++) {
        top = pPct[u] + pSD[u];
        if (top > maxV) maxV = top;
    }
    if (maxV <= 0) maxV = 1;
    maxV = maxV * 1.15;
    step = 1;
    mant = 1;
    while (maxV / step > 5) {
        if      (mant == 1) { mant = 2; step = step * 2; }
        else if (mant == 2) { mant = 5; step = (step / 2) * 5; }
        else                { mant = 1; step = step * 2; }
    }
    nTicks = 0;
    while (nTicks * step < maxV) nTicks++;
    maxV = nTicks * step;

    slot  = 118;    // wide: condition names are long
    axisW = 62;
    padR  = 22;
    gap   = 40;
    mt    = 46;
    mb    = 150;
    ph    = 380;

    totalW = 0;
    for (m = 0; m < nPan; m++) {
        pw = axisW + slot * pCount[m] + padR;
        if (pCount[m] == 0) pw = 0;
        totalW = totalW + pw;
        if (m < nPan - 1) totalW = totalW + gap;
    }
    W = totalW + 40;
    if (W < 700) W = 700;
    H = mt + ph + mb;

    newImage("Histogramme", "RGB white", W, H, 1);

    setJustification("center");

    xCur = 20;
    for (m = 0; m < nPan; m++) {

        if (pCount[m] == 0) continue;

        plotL = xCur + axisW;
        plotW = slot * pCount[m];
        yTop  = mt;
        yBot  = mt + ph;

        setColor(0, 0, 0);
        setFont("SansSerif", 14, "bold");
        setJustification("center");
        drawString(pTitle[m], plotL + plotW / 2, mt - 22);

        setLineWidth(2);
        drawLine(plotL, yTop, plotL, yBot);
        drawLine(plotL, yBot, plotL + plotW, yBot);

        setLineWidth(1);
        setFont("SansSerif", 12, "");
        setJustification("right");
        for (g = 0; g <= nTicks; g++) {
            yv = step * g;
            yy = yBot - round(ph * g / nTicks);
            setColor(205, 205, 205);
            drawLine(plotL + 1, yy, plotL + plotW, yy);
            setColor(0, 0, 0);
            drawLine(plotL - 6, yy, plotL, yy);
            drawString("" + yv, plotL - 11, yy + 5);
        }

        setJustification("center");
        bw = round(slot * 0.52);
        for (q = 0; q < pCount[m]; q++) {
            idx = pStart[m] + q;
            bx  = round(plotL + slot * q + (slot - bw) / 2);
            bh  = round(ph * pPct[idx] / maxV);

            if (pGlob[m] == 1) setColor(105, 105, 105);
            else               setColor(150, 150, 150);
            fillRect(bx, yBot - bh, bw, bh);
            setColor(70, 70, 70);
            setLineWidth(1);
            drawRect(bx, yBot - bh, bw, bh);

            cx = bx + round(bw / 2);
            if (pSD[idx] > 0) {
                hiV = pPct[idx] + pSD[idx];
                loV = pPct[idx] - pSD[idx];
                if (loV < 0) loV = 0;
                yHi = yBot - round(ph * hiV / maxV);
                yLo = yBot - round(ph * loV / maxV);
                setColor(40, 40, 40);
                setLineWidth(2);
                drawLine(cx, yHi, cx, yLo);
                drawLine(cx - 8, yHi, cx + 8, yHi);
                drawLine(cx - 8, yLo, cx + 8, yLo);
                yVal = yHi - 8;
            } else {
                yVal = yBot - bh - 8;
            }

            setColor(0, 0, 0);
            setFont("SansSerif", 13, "bold");
            drawString(d2s(pPct[idx], 2), cx, yVal);

            // condition name, wrapped onto two lines when too long
            setFont("SansSerif", 12, "bold");
            wrapLabel(pName[idx], slot - 10);
            drawString(wl1, cx, yBot + 22);
            if (wl2 != "") drawString(wl2, cx, yBot + 38);

            setFont("SansSerif", 11, "");
            setColor(60, 60, 60);
            drawString("SD = " + d2s(pSD[idx], 2), cx, yBot + 60);
            setColor(105, 105, 105);
            if (pGlob[m] == 1) {
                drawString("n = " + pNm[idx] + " exp.", cx, yBot + 78);
                drawString(pCells[idx] + " cells", cx, yBot + 94);
            } else {
                drawString("n = " + pCells[idx] + " cells", cx, yBot + 78);
                drawString("sur " + pNf[idx] + " fields", cx, yBot + 94);
            }
        }

        xCur = plotL + plotW + padR + gap;
    }

    saveAs("PNG", outPath);
    close();
}

// Split a label onto two lines when wider than maxW; result in wl1/wl2.
// The current font must be the drawing font: getStringWidth depends on it.
function wrapLabel(txt, maxW) {
    wl1 = txt;
    wl2 = "";
    if (getStringWidth(txt) <= maxW) return;
    parts = split(txt, " ");
    if (parts.length < 2) return;
    k = 1;
    for (q = 1; q < parts.length; q++) {
        cand = "";
        for (z = 0; z < q; z++) {
            if (z > 0) cand = cand + " ";
            cand = cand + parts[z];
        }
        if (getStringWidth(cand) <= maxW) k = q;
    }
    a = "";
    b = "";
    for (z = 0; z < parts.length; z++) {
        if (z < k) {
            if (z > 0) a = a + " ";
            a = a + parts[z];
        } else {
            if (z > k) b = b + " ";
            b = b + parts[z];
        }
    }
    wl1 = a;
    wl2 = b;
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
