(() => {
    const dropZone = document.getElementById('drop-zone');
    const fileInput = document.getElementById('file-input');
    const controlsSection = document.getElementById('controls-section');
    const resultsSection = document.getElementById('results-section');
    const resultCanvas = document.getElementById('result-canvas');
    const workCanvas = document.getElementById('work-canvas');
    const detectBtn = document.getElementById('detect-btn');
    const downloadBtn = document.getElementById('download-btn');
    const copyBtn = document.getElementById('copy-btn');
    const sensitivityInput = document.getElementById('sensitivity');
    const minHaloInput = document.getElementById('min-halo-size');
    const orientationInput = document.getElementById('plate-orientation');
    const sensitivityVal = document.getElementById('sensitivity-val');
    const minHaloSizeVal = document.getElementById('min-halo-size-val');
    const wellCount = document.getElementById('well-count');
    const plateMap = document.getElementById('plate-map');
    const resultsBody = document.querySelector('#results-table tbody');

    let currentImage = null;
    let lastResult = null;
    let showDebug = false;

    // --- File upload ---
    dropZone.addEventListener('click', () => fileInput.click());
    dropZone.addEventListener('dragover', e => {
        e.preventDefault();
        dropZone.classList.add('drag-over');
    });
    dropZone.addEventListener('dragleave', () => {
        dropZone.classList.remove('drag-over');
    });
    dropZone.addEventListener('drop', e => {
        e.preventDefault();
        dropZone.classList.remove('drag-over');
        if (e.dataTransfer.files.length > 0) {
            handleFile(e.dataTransfer.files[0]);
        }
    });
    fileInput.addEventListener('change', () => {
        if (fileInput.files.length > 0) {
            handleFile(fileInput.files[0]);
        }
    });

    function handleFile(file) {
        if (!file.type.startsWith('image/')) return;
        const reader = new FileReader();
        reader.onload = e => {
            const img = new Image();
            img.onload = () => {
                currentImage = img;
                controlsSection.classList.remove('hidden');
                resultsSection.classList.add('hidden');
                runDetection();
            };
            img.src = e.target.result;
        };
        reader.readAsDataURL(file);
    }

    // --- Controls ---
    sensitivityInput.addEventListener('input', () => {
        sensitivityVal.textContent = sensitivityInput.value;
    });
    minHaloInput.addEventListener('input', () => {
        minHaloSizeVal.textContent = minHaloInput.value;
    });
    detectBtn.addEventListener('click', runDetection);

    const debugToggle = document.getElementById('debug-toggle');
    if (debugToggle) {
        debugToggle.addEventListener('change', () => {
            showDebug = debugToggle.checked;
            if (lastResult) {
                drawAnnotated(currentImage, lastResult);
            }
        });
    }

    // --- Detection ---
    function runDetection() {
        if (!currentImage) return;

        detectBtn.textContent = 'Processing...';
        detectBtn.disabled = true;

        requestAnimationFrame(() => {
            setTimeout(() => {
                try {
                    doDetection();
                } finally {
                    detectBtn.textContent = 'Detect Halos';
                    detectBtn.disabled = false;
                }
            }, 50);
        });
    }

    function doDetection() {
        const img = currentImage;

        workCanvas.width = img.naturalWidth;
        workCanvas.height = img.naturalHeight;
        const wCtx = workCanvas.getContext('2d');
        wCtx.drawImage(img, 0, 0);
        const imageData = wCtx.getImageData(0, 0, img.naturalWidth, img.naturalHeight);

        const options = {
            sensitivity: parseInt(sensitivityInput.value),
            minHaloSize: parseInt(minHaloInput.value),
            orientation: orientationInput.value
        };
        const result = HaloDetector.detect(imageData, options);
        lastResult = result;

        drawAnnotated(img, result);
        buildResults(result, options.orientation);

        resultsSection.classList.remove('hidden');
    }

    // --- Annotated image ---
    function drawAnnotated(img, result) {
        const canvas = resultCanvas;
        canvas.width = img.naturalWidth;
        canvas.height = img.naturalHeight;
        const ctx = canvas.getContext('2d');

        ctx.drawImage(img, 0, 0);

        const { wells, colonyRadius, haloRadius, strictBlobs, gentleBlobs } = result;

        // Debug: show detected blobs
        if (showDebug) {
            // Gentle blobs (halo extent) in cyan
            if (gentleBlobs) {
                for (const blob of gentleBlobs) {
                    ctx.beginPath();
                    ctx.arc(blob.cx, blob.cy, blob.radius, 0, Math.PI * 2);
                    ctx.strokeStyle = 'rgba(0, 200, 255, 0.4)';
                    ctx.lineWidth = 1;
                    ctx.stroke();
                }
            }
            // Strict blobs (colony centers) in orange
            if (strictBlobs) {
                for (const blob of strictBlobs) {
                    ctx.beginPath();
                    ctx.arc(blob.cx, blob.cy, blob.radius, 0, Math.PI * 2);
                    ctx.strokeStyle = 'rgba(255, 165, 0, 0.7)';
                    ctx.lineWidth = 1.5;
                    ctx.stroke();

                    ctx.beginPath();
                    ctx.arc(blob.cx, blob.cy, 2, 0, Math.PI * 2);
                    ctx.fillStyle = 'orange';
                    ctx.fill();
                }
            }
        }

        // Draw well annotations
        for (const well of wells) {
            if (!well.hasColony && !showDebug) continue;

            const { cx, cy, hasHalo, hasColony, gentleRadius } = well;

            if (hasHalo) {
                // Draw green circle sized to the gentle blob (actual halo extent)
                const drawR = gentleRadius || haloRadius;
                ctx.beginPath();
                ctx.arc(cx, cy, drawR, 0, Math.PI * 2);
                ctx.strokeStyle = '#22c55e';
                ctx.lineWidth = 3;
                ctx.stroke();

                // Well label above
                ctx.font = `bold ${Math.max(12, colonyRadius * 0.8)}px sans-serif`;
                ctx.fillStyle = '#22c55e';
                ctx.textAlign = 'center';
                ctx.fillText(well.wellName, cx, cy - drawR - 4);
            } else if (hasColony) {
                // Gray circle for non-halo colonies
                ctx.beginPath();
                ctx.arc(cx, cy, colonyRadius * 1.2, 0, Math.PI * 2);
                ctx.strokeStyle = 'rgba(148, 163, 184, 0.5)';
                ctx.lineWidth = 1.5;
                ctx.stroke();
            }

            // Debug: show grid crosshair and well name for all positions
            if (showDebug) {
                const crossSize = 4;
                ctx.strokeStyle = 'rgba(67, 97, 238, 0.5)';
                ctx.lineWidth = 1;
                ctx.beginPath();
                ctx.moveTo(cx - crossSize, cy);
                ctx.lineTo(cx + crossSize, cy);
                ctx.moveTo(cx, cy - crossSize);
                ctx.lineTo(cx, cy + crossSize);
                ctx.stroke();

                ctx.font = '9px sans-serif';
                ctx.fillStyle = 'rgba(67, 97, 238, 0.7)';
                ctx.textAlign = 'center';
                ctx.fillText(well.wellName, cx, cy + haloRadius + 10);
            }
        }

        // Debug info overlay
        if (showDebug && result.gridAngle !== undefined) {
            ctx.font = '14px monospace';
            ctx.fillStyle = 'rgba(255, 255, 255, 0.9)';
            ctx.textAlign = 'left';
            const angleDeg = (result.gridAngle * 180 / Math.PI).toFixed(1);
            const info = [
                `Grid angle: ${angleDeg}deg  |  RANSAC inliers: ${result.inlierCount ?? '?'}/${result.totalBlobs ?? '?'}`,
                `Strict blobs: ${strictBlobs ? strictBlobs.length : 0}  Gentle blobs: ${gentleBlobs ? gentleBlobs.length : 0}`,
                `H-spacing: ${result.hSpacing?.toFixed(1)}px  V-spacing: ${result.vSpacing?.toFixed(1)}px`,
                `Size ratio threshold: ${result.profileThreshold ?? '?'}`
            ];
            // Background for readability
            ctx.fillStyle = 'rgba(0, 0, 0, 0.6)';
            ctx.fillRect(5, 5, 450, info.length * 18 + 8);
            ctx.fillStyle = 'rgba(255, 255, 255, 0.95)';
            info.forEach((line, i) => ctx.fillText(line, 10, 20 + i * 18));
        }
    }

    // --- Results table and plate map ---
    function buildResults(result, orientation) {
        const { wells } = result;
        const haloWells = wells.filter(w => w.hasHalo);

        wellCount.textContent = haloWells.length;

        haloWells.sort((a, b) => a.wellName.localeCompare(b.wellName, undefined, { numeric: true }));

        resultsBody.innerHTML = '';
        for (const w of haloWells) {
            const tr = document.createElement('tr');
            tr.innerHTML = `
                <td><strong>${w.wellName}</strong></td>
                <td>${w.rowLabel}</td>
                <td>${w.colNum}</td>
                <td>${w.profileScore ?? w.haloScore}</td>
            `;
            resultsBody.appendChild(tr);
        }

        const numRows = 8;
        const numCols = 12;
        const rowLabels = 'ABCDEFGH'.split('');

        plateMap.innerHTML = '';

        const thead = document.createElement('thead');
        const headRow = document.createElement('tr');
        headRow.innerHTML = '<th></th>';
        for (let c = 1; c <= numCols; c++) {
            headRow.innerHTML += `<th>${c}</th>`;
        }
        thead.appendChild(headRow);
        plateMap.appendChild(thead);

        const tbody = document.createElement('tbody');
        for (let r = 0; r < numRows; r++) {
            const tr = document.createElement('tr');
            tr.innerHTML = `<th>${rowLabels[r]}</th>`;

            for (let c = 0; c < numCols; c++) {
                const well = wells.find(w => w.row === r && w.col === c);
                const td = document.createElement('td');
                if (well) {
                    if (well.hasHalo) {
                        td.className = 'halo-well';
                        td.textContent = '\u2713';
                        td.title = `${well.wellName} - Size ratio: ${well.sizeRatio}`;
                    } else if (well.hasColony) {
                        td.className = 'normal-colony';
                        td.textContent = '\u25CB';
                        td.title = `${well.wellName} - Colony (ratio: ${well.sizeRatio})`;
                    } else {
                        td.className = 'no-colony';
                        td.textContent = '\u00B7';
                        td.title = `${well.wellName} - Empty`;
                    }
                }
                tr.appendChild(td);
            }
            tbody.appendChild(tr);
        }
        plateMap.appendChild(tbody);
    }

    // --- Download ---
    downloadBtn.addEventListener('click', () => {
        const link = document.createElement('a');
        link.download = 'annotated_plate.png';
        link.href = resultCanvas.toDataURL('image/png');
        link.click();
    });

    // --- Copy well list ---
    copyBtn.addEventListener('click', () => {
        const rows = resultsBody.querySelectorAll('tr');
        const lines = ['Well\tRow\tColumn\tSize Ratio'];
        rows.forEach(row => {
            const cells = row.querySelectorAll('td');
            lines.push(Array.from(cells).map(c => c.textContent.trim()).join('\t'));
        });
        navigator.clipboard.writeText(lines.join('\n')).then(() => {
            const orig = copyBtn.textContent;
            copyBtn.textContent = 'Copied!';
            setTimeout(() => { copyBtn.textContent = orig; }, 2000);
        });
    });
})();
