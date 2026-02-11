/**
 * Plate Halo Detector v4
 *
 * Uses RANSAC-affine grid fitting:
 * 1. Find colony blobs at strict & gentle thresholds
 * 2. Estimate grid spacing from nearest-neighbour distances
 * 3. Use Hough-style vector voting to find the two grid basis vectors
 * 4. RANSAC: pick blob triplets, hypothesize grid labels from basis vectors,
 *    solve affine transform, count inliers — keep the best
 * 5. Refine with least-squares on all inliers
 * 6. Generate all 96 well positions from the refined affine
 * 7. Halo detection by comparing gentle vs strict blob sizes
 */

const HaloDetector = (() => {

    // ============================================================
    //  Image processing
    // ============================================================

    function toGrayscale(imageData) {
        const { data, width, height } = imageData;
        const gray = new Float32Array(width * height);
        for (let i = 0; i < width * height; i++) {
            gray[i] = 0.299 * data[i * 4] + 0.587 * data[i * 4 + 1] + 0.114 * data[i * 4 + 2];
        }
        return gray;
    }

    function gaussianBlur(gray, w, h, radius) {
        const sigma = Math.max(radius / 2, 1);
        const ks = radius * 2 + 1;
        const kernel = new Float32Array(ks);
        let sum = 0;
        for (let i = 0; i < ks; i++) {
            const x = i - radius;
            kernel[i] = Math.exp(-(x * x) / (2 * sigma * sigma));
            sum += kernel[i];
        }
        for (let i = 0; i < ks; i++) kernel[i] /= sum;

        const temp = new Float32Array(w * h);
        for (let y = 0; y < h; y++) {
            for (let x = 0; x < w; x++) {
                let val = 0;
                for (let k = -radius; k <= radius; k++) {
                    val += gray[y * w + Math.min(Math.max(x + k, 0), w - 1)] * kernel[k + radius];
                }
                temp[y * w + x] = val;
            }
        }
        const result = new Float32Array(w * h);
        for (let y = 0; y < h; y++) {
            for (let x = 0; x < w; x++) {
                let val = 0;
                for (let k = -radius; k <= radius; k++) {
                    val += temp[Math.min(Math.max(y + k, 0), h - 1) * w + x] * kernel[k + radius];
                }
                result[y * w + x] = val;
            }
        }
        return result;
    }

    function localMean(gray, w, h, blockRadius) {
        const integral = new Float64Array(w * h);
        for (let y = 0; y < h; y++) {
            let rowSum = 0;
            for (let x = 0; x < w; x++) {
                rowSum += gray[y * w + x];
                integral[y * w + x] = rowSum + (y > 0 ? integral[(y - 1) * w + x] : 0);
            }
        }
        const getI = (x, y) => (x < 0 || y < 0) ? 0 : integral[Math.min(y, h - 1) * w + Math.min(x, w - 1)];
        const result = new Float32Array(w * h);
        for (let y = 0; y < h; y++) {
            for (let x = 0; x < w; x++) {
                const x0 = Math.max(x - blockRadius - 1, -1);
                const y0 = Math.max(y - blockRadius - 1, -1);
                const x1 = Math.min(x + blockRadius, w - 1);
                const y1 = Math.min(y + blockRadius, h - 1);
                const area = (x1 - x0) * (y1 - y0);
                result[y * w + x] = (getI(x1, y1) - getI(x0, y1) - getI(x1, y0) + getI(x0, y0)) / area;
            }
        }
        return result;
    }

    // ============================================================
    //  Connected components
    // ============================================================

    function connectedComponents(binary, w, h) {
        const labels = new Int32Array(w * h);
        let nextLabel = 1;
        const parent = [0];
        function find(x) { while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; }
        function union(a, b) { a = find(a); b = find(b); if (a !== b) parent[Math.max(a, b)] = Math.min(a, b); }

        for (let y = 0; y < h; y++) {
            for (let x = 0; x < w; x++) {
                if (!binary[y * w + x]) continue;
                const up = y > 0 ? labels[(y - 1) * w + x] : 0;
                const left = x > 0 ? labels[y * w + x - 1] : 0;
                if (!up && !left) { labels[y * w + x] = nextLabel; parent.push(nextLabel); nextLabel++; }
                else if (up && !left) { labels[y * w + x] = up; }
                else if (!up && left) { labels[y * w + x] = left; }
                else { labels[y * w + x] = Math.min(up, left); union(up, left); }
            }
        }

        const components = new Map();
        for (let y = 0; y < h; y++) {
            for (let x = 0; x < w; x++) {
                const l = labels[y * w + x]; if (!l) continue;
                const root = find(l);
                if (!components.has(root)) components.set(root, { sumX: 0, sumY: 0, count: 0, minX: w, minY: h, maxX: 0, maxY: 0 });
                const c = components.get(root);
                c.sumX += x; c.sumY += y; c.count++;
                c.minX = Math.min(c.minX, x); c.minY = Math.min(c.minY, y);
                c.maxX = Math.max(c.maxX, x); c.maxY = Math.max(c.maxY, y);
            }
        }
        return components;
    }

    // ============================================================
    //  Blob detection
    // ============================================================

    function findBlobs(blurred, w, h) {
        const blockRadius = Math.max(Math.round(Math.min(w, h) * 0.06), 15);
        const mean = localMean(blurred, w, h, blockRadius);

        const binaryStrict = new Uint8Array(w * h);
        for (let i = 0; i < w * h; i++) binaryStrict[i] = (mean[i] - blurred[i]) > 12 ? 1 : 0;

        const binaryGentle = new Uint8Array(w * h);
        for (let i = 0; i < w * h; i++) binaryGentle[i] = (mean[i] - blurred[i]) > 3 ? 1 : 0;

        const minDim = Math.min(w, h);
        const minArea = Math.pow(minDim * 0.008, 2) * Math.PI;
        const maxArea = Math.pow(minDim * 0.07, 2) * Math.PI;

        function extractBlobs(comps, minA, maxA, maxAspect) {
            const blobs = [];
            for (const [, c] of comps) {
                if (c.count < minA || c.count > maxA) continue;
                const bw = c.maxX - c.minX + 1, bh = c.maxY - c.minY + 1;
                if (Math.max(bw, bh) / Math.min(bw, bh) > maxAspect) continue;
                blobs.push({ cx: c.sumX / c.count, cy: c.sumY / c.count, area: c.count, radius: Math.sqrt(c.count / Math.PI) });
            }
            return blobs;
        }

        return {
            strictBlobs: extractBlobs(connectedComponents(binaryStrict, w, h), minArea, maxArea, 2.5),
            gentleBlobs: extractBlobs(connectedComponents(binaryGentle, w, h), minArea, maxArea * 4, 3)
        };
    }

    // ============================================================
    //  Math utilities
    // ============================================================

    function median(arr) {
        if (!arr.length) return 0;
        const s = arr.slice().sort((a, b) => a - b);
        const m = s.length >> 1;
        return s.length & 1 ? s[m] : (s[m - 1] + s[m]) / 2;
    }

    function invert3x3(M) {
        const [[a, b, c], [d, e, f], [g, h, i]] = M;
        const det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
        if (Math.abs(det) < 1e-12) return null;
        const inv = 1 / det;
        return [
            [(e * i - f * h) * inv, (c * h - b * i) * inv, (b * f - c * e) * inv],
            [(f * g - d * i) * inv, (a * i - c * g) * inv, (c * d - a * f) * inv],
            [(d * h - e * g) * inv, (b * g - a * h) * inv, (a * e - b * d) * inv]
        ];
    }

    function solveAffineLeastSquares(matches, blobs) {
        // matches = [{r, c, blobIdx}, ...]
        // Solve affine: x = a*c + b*r + tx,  y = d*c + e*r + ty
        const ATA = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
        const ATbx = [0, 0, 0], ATby = [0, 0, 0];

        for (const m of matches) {
            const row = [m.c, m.r, 1];
            const bx = blobs[m.blobIdx].cx, by = blobs[m.blobIdx].cy;
            for (let i = 0; i < 3; i++) {
                ATbx[i] += row[i] * bx;
                ATby[i] += row[i] * by;
                for (let j = 0; j < 3; j++) ATA[i][j] += row[i] * row[j];
            }
        }

        const inv = invert3x3(ATA);
        if (!inv) return null;

        return {
            a:  inv[0][0] * ATbx[0] + inv[0][1] * ATbx[1] + inv[0][2] * ATbx[2],
            b:  inv[1][0] * ATbx[0] + inv[1][1] * ATbx[1] + inv[1][2] * ATbx[2],
            tx: inv[2][0] * ATbx[0] + inv[2][1] * ATbx[1] + inv[2][2] * ATbx[2],
            d:  inv[0][0] * ATby[0] + inv[0][1] * ATby[1] + inv[0][2] * ATby[2],
            e:  inv[1][0] * ATby[0] + inv[1][1] * ATby[1] + inv[1][2] * ATby[2],
            ty: inv[2][0] * ATby[0] + inv[2][1] * ATby[1] + inv[2][2] * ATby[2]
        };
    }

    function applyAffine(af, r, c) {
        return { x: af.a * c + af.b * r + af.tx, y: af.d * c + af.e * r + af.ty };
    }

    // ============================================================
    //  RANSAC Grid Fitting
    // ============================================================

    function fitGrid(blobs, numRows, numCols, imgW, imgH) {
        if (blobs.length < 10) return fallbackGrid(numRows, numCols, imgW, imgH);

        // --- Step 1: Estimate grid spacing ---
        const nnDists = blobs.map((b, i) => {
            let best = Infinity;
            for (let j = 0; j < blobs.length; j++) {
                if (i === j) continue;
                const d = Math.hypot(b.cx - blobs[j].cx, b.cy - blobs[j].cy);
                if (d < best) best = d;
            }
            return best;
        });
        const spacing = median(nnDists);

        // --- Step 2: Find grid basis vectors via Hough-style voting ---
        // Collect all pairwise vectors that are close to 1x grid spacing
        const candVecs = [];
        for (let i = 0; i < blobs.length; i++) {
            for (let j = i + 1; j < blobs.length; j++) {
                const dx = blobs[j].cx - blobs[i].cx;
                const dy = blobs[j].cy - blobs[i].cy;
                const d = Math.hypot(dx, dy);
                if (d > spacing * 0.6 && d < spacing * 1.4) {
                    candVecs.push({ dx, dy });
                }
            }
        }

        // Bin vectors into angular histogram, find the dominant angle
        const nBins = 360;
        const aHist = new Float32Array(nBins);
        for (const v of candVecs) {
            let a = Math.atan2(v.dy, v.dx);
            if (a < 0) a += Math.PI;           // fold to [0, PI)
            if (a >= Math.PI) a -= Math.PI;
            const bin = Math.floor(a / Math.PI * nBins) % nBins;
            aHist[bin]++;
            aHist[(bin + 1) % nBins] += 0.5;
            aHist[(bin - 1 + nBins) % nBins] += 0.5;
        }
        let bestBin = 0;
        for (let i = 1; i < nBins; i++) if (aHist[i] > aHist[bestBin]) bestBin = i;
        const dominantAngle = (bestBin + 0.5) / nBins * Math.PI;

        // Determine which is the "column" direction (more horizontal) and "row" direction
        let colAngle, rowAngle;
        if (dominantAngle < Math.PI / 4 || dominantAngle > 3 * Math.PI / 4) {
            colAngle = dominantAngle;
            rowAngle = dominantAngle + Math.PI / 2;
        } else {
            rowAngle = dominantAngle;
            colAngle = dominantAngle - Math.PI / 2;
        }
        // Normalize colAngle to point rightward
        if (Math.cos(colAngle) < 0) colAngle += Math.PI;
        // Normalize rowAngle to point downward
        if (Math.sin(rowAngle) < 0) rowAngle += Math.PI;

        // Compute mean basis vectors from candidates near each direction
        function meanVecNear(angle, vecs, spacingEst) {
            const cos0 = Math.cos(angle), sin0 = Math.sin(angle);
            let sx = 0, sy = 0, n = 0;
            for (const v of vecs) {
                // Project onto the target direction
                const proj = v.dx * cos0 + v.dy * sin0;
                if (proj < spacingEst * 0.5) continue; // wrong direction
                const perp = Math.abs(-v.dx * sin0 + v.dy * cos0);
                if (perp > spacingEst * 0.3) continue; // too far off axis
                sx += v.dx; sy += v.dy; n++;
            }
            return n > 0 ? { dx: sx / n, dy: sy / n } : { dx: spacingEst * cos0, dy: spacingEst * sin0 };
        }

        const colVec = meanVecNear(colAngle, candVecs, spacing);
        const rowVec = meanVecNear(rowAngle, candVecs, spacing);

        // --- Step 3: RANSAC ---
        const epsilon = spacing * 0.35;
        const maxIter = 300;
        let bestInlierCount = 0;
        let bestAffine = null;
        let bestMatches = [];

        // Helper: given a blob index and the two basis vectors, compute grid (r,c)
        // relative to another blob
        function estimateRC(fromBlob, toBlob) {
            const dx = blobs[toBlob].cx - blobs[fromBlob].cx;
            const dy = blobs[toBlob].cy - blobs[fromBlob].cy;
            // Solve: dx = dc * colVec.dx + dr * rowVec.dx
            //        dy = dc * colVec.dy + dr * rowVec.dy
            const det = colVec.dx * rowVec.dy - colVec.dy * rowVec.dx;
            if (Math.abs(det) < 1e-8) return null;
            const dc = (dx * rowVec.dy - dy * rowVec.dx) / det;
            const dr = (dy * colVec.dx - dx * colVec.dy) / det;
            return { dr: Math.round(dr), dc: Math.round(dc) };
        }

        function countInliers(affine) {
            const matched = [];
            const usedBlobs = new Set();
            for (let r = 0; r < numRows; r++) {
                for (let c = 0; c < numCols; c++) {
                    const { x, y } = applyAffine(affine, r, c);
                    let bestDist = epsilon, bestIdx = -1;
                    for (let i = 0; i < blobs.length; i++) {
                        if (usedBlobs.has(i)) continue;
                        const d = Math.hypot(blobs[i].cx - x, blobs[i].cy - y);
                        if (d < bestDist) { bestDist = d; bestIdx = i; }
                    }
                    if (bestIdx >= 0) {
                        matched.push({ r, c, blobIdx: bestIdx, dist: bestDist });
                        usedBlobs.add(bestIdx);
                    }
                }
            }
            return matched;
        }

        // Seeded RNG for reproducibility
        let rngState = 42;
        function rng() { rngState = (rngState * 1664525 + 1013904223) & 0x7fffffff; return rngState / 0x7fffffff; }
        function randInt(n) { return Math.floor(rng() * n); }

        for (let iter = 0; iter < maxIter; iter++) {
            // Pick anchor blob
            const ai = randInt(blobs.length);

            // Pick two other blobs that are neighbours
            const bi = randInt(blobs.length);
            if (bi === ai) continue;
            const ci = randInt(blobs.length);
            if (ci === ai || ci === bi) continue;

            // Estimate grid offsets relative to anchor
            const rcB = estimateRC(ai, bi);
            const rcC = estimateRC(ai, ci);
            if (!rcB || !rcC) continue;

            // Skip degenerate cases (all on same row or column)
            if (rcB.dr === 0 && rcC.dr === 0) continue;
            if (rcB.dc === 0 && rcC.dc === 0) continue;

            // Try different anchor positions: the anchor could be at any grid cell.
            // We try a few candidates that would place the grid reasonably.
            // Project anchor position onto grid axes to guess its (r,c)
            const projC = (blobs[ai].cx * colVec.dx + blobs[ai].cy * colVec.dy) / (colVec.dx ** 2 + colVec.dy ** 2);
            const projR = (blobs[ai].cx * rowVec.dx + blobs[ai].cy * rowVec.dy) / (rowVec.dx ** 2 + rowVec.dy ** 2);

            // Find what (r0,c0) for the anchor would place the top-left grid cell
            // closest to the image top-left
            const minProjC = projC - (numCols - 1);
            const minProjR = projR - (numRows - 1);

            // Try the anchor at a few plausible positions
            for (let anchorR = 0; anchorR < numRows; anchorR += Math.max(1, numRows - 1)) {
                for (let anchorC = 0; anchorC < numCols; anchorC += Math.max(1, numCols - 1)) {
                    const src = [
                        { r: anchorR, c: anchorC },
                        { r: anchorR + rcB.dr, c: anchorC + rcB.dc },
                        { r: anchorR + rcC.dr, c: anchorC + rcC.dc }
                    ];

                    // Check bounds
                    if (src.some(s => s.r < 0 || s.r >= numRows || s.c < 0 || s.c >= numCols)) continue;

                    const dst = [
                        { x: blobs[ai].cx, y: blobs[ai].cy },
                        { x: blobs[bi].cx, y: blobs[bi].cy },
                        { x: blobs[ci].cx, y: blobs[ci].cy }
                    ];

                    // Solve affine from these 3 correspondences
                    const affine = solveAffineLeastSquares(
                        [{ r: src[0].r, c: src[0].c, blobIdx: ai },
                         { r: src[1].r, c: src[1].c, blobIdx: bi },
                         { r: src[2].r, c: src[2].c, blobIdx: ci }],
                        blobs
                    );
                    if (!affine) continue;

                    // Sanity check: the affine should produce spacings close to our estimate
                    const p00 = applyAffine(affine, 0, 0);
                    const p01 = applyAffine(affine, 0, 1);
                    const p10 = applyAffine(affine, 1, 0);
                    const colSpacing = Math.hypot(p01.x - p00.x, p01.y - p00.y);
                    const rowSpacing = Math.hypot(p10.x - p00.x, p10.y - p00.y);
                    if (colSpacing < spacing * 0.5 || colSpacing > spacing * 2) continue;
                    if (rowSpacing < spacing * 0.5 || rowSpacing > spacing * 2) continue;

                    const matched = countInliers(affine);
                    if (matched.length > bestInlierCount) {
                        bestInlierCount = matched.length;
                        bestAffine = affine;
                        bestMatches = matched;
                    }
                }
            }

            // Early exit if we matched almost everything
            if (bestInlierCount >= blobs.length - 2 && bestInlierCount >= 80) break;
        }

        // --- Step 4: Refine with all inliers ---
        if (bestMatches.length >= 3) {
            const refined = solveAffineLeastSquares(bestMatches, blobs);
            if (refined) {
                bestAffine = refined;
                // Re-count with refined affine
                bestMatches = countInliers(bestAffine);
                bestInlierCount = bestMatches.length;
            }
        }

        // --- Step 5: Generate final grid ---
        if (!bestAffine) return fallbackGrid(numRows, numCols, imgW, imgH);

        const wells = [];
        for (let r = 0; r < numRows; r++) {
            for (let c = 0; c < numCols; c++) {
                const { x, y } = applyAffine(bestAffine, r, c);
                wells.push({ gridRow: r, gridCol: c, cx: x, cy: y });
            }
        }

        // Compute spacings for halo sizing
        const p00 = applyAffine(bestAffine, 0, 0);
        const p01 = applyAffine(bestAffine, 0, 1);
        const p10 = applyAffine(bestAffine, 1, 0);

        return {
            wells,
            gridAngle: Math.atan2(p01.y - p00.y, p01.x - p00.x),
            hSpacing: Math.hypot(p01.x - p00.x, p01.y - p00.y),
            vSpacing: Math.hypot(p10.x - p00.x, p10.y - p00.y),
            medianColonyRadius: median(blobs.map(b => b.radius)),
            inlierCount: bestInlierCount,
            totalBlobs: blobs.length,
            affine: bestAffine
        };
    }

    function fallbackGrid(numRows, numCols, imgW, imgH) {
        const margin = Math.min(imgW, imgH) * 0.08;
        const cellW = (imgW - 2 * margin) / numCols;
        const cellH = (imgH - 2 * margin) / numRows;
        const wells = [];
        for (let r = 0; r < numRows; r++) {
            for (let c = 0; c < numCols; c++) {
                wells.push({ gridRow: r, gridCol: c, cx: margin + (c + 0.5) * cellW, cy: margin + (r + 0.5) * cellH });
            }
        }
        return { wells, gridAngle: 0, hSpacing: cellW, vSpacing: cellH, medianColonyRadius: Math.min(cellW, cellH) * 0.25, inlierCount: 0, totalBlobs: 0, affine: null };
    }

    // ============================================================
    //  Halo detection — radial intensity profile approach
    // ============================================================

    /**
     * Sample average intensity in a ring around (cx, cy).
     */
    function ringMean(gray, w, h, cx, cy, innerR, outerR) {
        let sum = 0, count = 0;
        const x0 = Math.max(0, Math.floor(cx - outerR));
        const x1 = Math.min(w - 1, Math.ceil(cx + outerR));
        const y0 = Math.max(0, Math.floor(cy - outerR));
        const y1 = Math.min(h - 1, Math.ceil(cy + outerR));
        for (let y = y0; y <= y1; y++) {
            for (let x = x0; x <= x1; x++) {
                const d = Math.hypot(x - cx, y - cy);
                if (d >= innerR && d <= outerR) {
                    sum += gray[y * w + x];
                    count++;
                }
            }
        }
        return count > 0 ? sum / count : 0;
    }

    function analyzeHalos(blurred, w, h, gridResult, strictBlobs, gentleBlobs, options) {
        const { sensitivity = 10 } = options;
        const { wells, hSpacing, vSpacing, medianColonyRadius } = gridResult;

        const cellSize = Math.min(Math.abs(hSpacing), Math.abs(vSpacing));
        const matchRadius = cellSize * 0.45;
        const colonyR = medianColonyRadius;

        function findNearestBlob(blobList, cx, cy, maxDist) {
            let best = null, bestD = maxDist;
            for (const b of blobList) {
                const d = Math.hypot(b.cx - cx, b.cy - cy);
                if (d < bestD) { bestD = d; best = b; }
            }
            return best;
        }

        const medianStrictR = median(strictBlobs.map(b => b.radius));
        const rowLabels = 'ABCDEFGH'.split('');

        // First pass: compute radial profile score for every colony
        const wellData = [];

        for (const well of wells) {
            const { cx, cy, gridRow, gridCol } = well;
            if (cx < 0 || cx >= w || cy < 0 || cy >= h) continue;

            const strictMatch = findNearestBlob(strictBlobs, cx, cy, matchRadius);
            const hasColony = strictMatch !== null;

            let profileScore = 0;
            let peakBrightness = 0;
            let bgBrightness = 0;

            if (hasColony) {
                // Sample radial intensity profile in concentric rings
                // from colony edge out to cell boundary
                const nRings = 10;
                const startR = colonyR * 1.1;
                const endR = cellSize * 0.48;
                const ringValues = [];

                for (let i = 0; i < nRings; i++) {
                    const r0 = startR + (endR - startR) * (i / nRings);
                    const r1 = startR + (endR - startR) * ((i + 1) / nRings);
                    ringValues.push(ringMean(blurred, w, h, cx, cy, r0, r1));
                }

                // Colony center brightness
                const colonyBright = ringMean(blurred, w, h, cx, cy, 0, colonyR * 0.8);

                // Background = outermost 3 rings average
                bgBrightness = (ringValues[nRings - 1] + ringValues[nRings - 2] + ringValues[nRings - 3]) / 3;

                // Peak brightness in the halo zone (inner half of rings)
                const haloRings = ringValues.slice(0, Math.ceil(nRings * 0.6));
                peakBrightness = Math.max(...haloRings);

                // Halo score = how much brighter the peak halo ring is vs background
                // A colony with a halo has: dark center, BRIGHT halo ring, then dimmer background
                // A colony without a halo has: dark center, then gradually reaches background
                profileScore = peakBrightness - bgBrightness;
            }

            // Also compute the gentle blob size ratio as before
            const gentleMatch = findNearestBlob(gentleBlobs, cx, cy, matchRadius);
            const gentleR = gentleMatch ? gentleMatch.radius : 0;
            const sizeRatio = hasColony && medianStrictR > 0 ? gentleR / medianStrictR : 0;

            wellData.push({
                well, gridRow, gridCol, cx, cy,
                hasColony, strictMatch, gentleR,
                sizeRatio, profileScore, peakBrightness, bgBrightness
            });
        }

        // Compute threshold from profile scores of all colonies
        // The idea: colonies without halos have profileScore near 0 or slightly negative
        // Colonies with halos have a clear positive profileScore
        const colonyScores = wellData.filter(d => d.hasColony).map(d => d.profileScore);
        colonyScores.sort((a, b) => a - b);

        // Use sensitivity to set the threshold
        // At sensitivity=10 (default): threshold at roughly the 60th percentile
        // Higher sensitivity = lower threshold = more detections
        const percentile = Math.max(0.1, Math.min(0.9, 1.0 - sensitivity / 30));
        const autoThreshold = colonyScores[Math.floor(colonyScores.length * percentile)] || 2;
        // Clamp to a minimum of 1.5 to avoid noise
        const profileThreshold = Math.max(1.5, autoThreshold);

        const results = [];
        for (const wd of wellData) {
            const { well, gridRow, gridCol, cx, cy, hasColony, gentleR, sizeRatio, profileScore } = wd;

            const hasHalo = hasColony && profileScore > profileThreshold;

            results.push({
                row: gridRow, col: gridCol,
                cx: Math.round(cx), cy: Math.round(cy),
                wellName: `${rowLabels[gridRow] || '?'}${gridCol + 1}`,
                rowLabel: rowLabels[gridRow] || '?', colNum: gridCol + 1,
                hasColony, hasHalo,
                haloScore: hasColony ? Math.round(profileScore * 10) / 10 : 0,
                gentleRadius: Math.round(gentleR),
                strictRadius: wd.strictMatch ? Math.round(wd.strictMatch.radius) : 0,
                sizeRatio: Math.round(sizeRatio * 100) / 100,
                profileScore: Math.round(profileScore * 10) / 10
            });
        }

        return {
            wells: results,
            colonyRadius: Math.round(medianStrictR),
            haloRadius: Math.round(cellSize * 0.45),
            medianStrictR,
            profileThreshold: Math.round(profileThreshold * 10) / 10
        };
    }

    // ============================================================
    //  Main entry
    // ============================================================

    function detect(imageData, options = {}) {
        const { orientation = 'landscape' } = options;
        const w = imageData.width, h = imageData.height;

        const gray = toGrayscale(imageData);
        const blurred = gaussianBlur(gray, w, h, 3);

        const { strictBlobs, gentleBlobs } = findBlobs(blurred, w, h);

        const numRows = orientation === 'landscape' ? 8 : 12;
        const numCols = orientation === 'landscape' ? 12 : 8;
        const gridResult = fitGrid(strictBlobs, numRows, numCols, w, h);

        const haloResult = analyzeHalos(blurred, w, h, gridResult, strictBlobs, gentleBlobs, options);

        return {
            plateRegion: { x: 0, y: 0, w, h },
            ...haloResult,
            gridAngle: gridResult.gridAngle,
            hSpacing: gridResult.hSpacing,
            vSpacing: gridResult.vSpacing,
            inlierCount: gridResult.inlierCount,
            totalBlobs: gridResult.totalBlobs,
            strictBlobs, gentleBlobs
        };
    }

    return { detect };
})();
