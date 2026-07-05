/* KmerLab front-end. Vanilla JS, no frameworks. Talks to the local Flask API. */
const KmerLab = (function () {
  "use strict";

  function $(id) { return document.getElementById(id); }

  const BASE_COLORS = {
    A: "var(--base-a)", C: "var(--base-c)",
    G: "var(--base-g)", T: "var(--base-t)", other: "var(--base-n)",
  };

  function showError(msg) {
    const panel = $("error-panel");
    panel.textContent = msg;
    panel.classList.remove("hidden");
    panel.scrollIntoView({ behavior: "smooth", block: "center" });
  }
  function clearError() { $("error-panel").classList.add("hidden"); }
  function setLoading(on) { $("loading").classList.toggle("hidden", !on); }
  function fmt(n) { return typeof n === "number" ? n.toLocaleString() : n; }

  // Wrap each base of a k-mer/sequence in a colour-coded span.
  function colorizeSeq(seq) {
    let out = "";
    for (const ch of seq.toUpperCase()) {
      const cls = "ACGT".includes(ch) ? "b-" + ch : "b-other";
      out += `<span class="${cls}">${ch}</span>`;
    }
    return out;
  }

  // ------------------------------------------------------------- Hero ribbon
  function initRibbon() {
    const el = $("ribbon");
    if (!el) return;
    const bases = "ACGT";
    let seq = "";
    for (let i = 0; i < 140; i++) seq += bases[Math.floor(Math.random() * 4)];
    // Duplicate so the -50% scroll animation loops seamlessly.
    el.innerHTML = colorizeSeq(seq) + colorizeSeq(seq);
  }

  // ------------------------------------------------------------- Dropzone UX
  function wireDropzone(zoneId, inputId) {
    const zone = $(zoneId);
    const input = $(inputId);
    if (!zone || !input) return;
    const title = zone.querySelector(".dz-title");
    const original = title ? title.innerHTML : "";

    function showName(name) {
      if (title) title.innerHTML = `<span class="dz-file">✓ ${name}</span>`;
    }
    input.addEventListener("change", () => {
      if (input.files && input.files.length) showName(input.files[0].name);
      else if (title) title.innerHTML = original;
    });
    ["dragenter", "dragover"].forEach((ev) =>
      zone.addEventListener(ev, (e) => { e.preventDefault(); zone.classList.add("drag"); })
    );
    ["dragleave", "drop"].forEach((ev) =>
      zone.addEventListener(ev, (e) => { e.preventDefault(); zone.classList.remove("drag"); })
    );
    zone.addEventListener("drop", (e) => {
      if (e.dataTransfer.files && e.dataTransfer.files.length) {
        input.files = e.dataTransfer.files;
        showName(input.files[0].name);
      }
    });
  }

  async function loadSampleList(selectEl) {
    try {
      const res = await fetch("/api/samples");
      const files = await res.json();
      files.forEach((f) => {
        const opt = document.createElement("option");
        opt.value = f; opt.textContent = f;
        selectEl.appendChild(opt);
      });
    } catch (e) { /* samples optional */ }
  }

  async function attachInput(formData, fileEl, sampleEl, field) {
    if (fileEl.files && fileEl.files.length > 0) {
      formData.append(field, fileEl.files[0]);
      return true;
    }
    if (sampleEl && sampleEl.value) {
      const res = await fetch("/samples/" + encodeURIComponent(sampleEl.value));
      if (!res.ok) throw new Error("Could not load sample: " + sampleEl.value);
      formData.append(field + "_text", await res.text());
      return true;
    }
    return false;
  }

  function card(label, value, sub) {
    return `<div class="card"><div class="card-value">${value}</div>` +
           `<div class="card-label">${label}</div>` +
           (sub ? `<p>${sub}</p>` : "") + `</div>`;
  }

  function renderComposition(comp) {
    const box = $("composition");
    if (!box || !comp) return;
    const order = ["A", "C", "G", "T", "other"];
    const frac = comp.fractions;
    const counts = comp.counts;
    const bar = order
      .filter((b) => frac[b] > 0)
      .map((b) => {
        const pct = (frac[b] * 100);
        const tiny = pct < 6 ? " tiny" : "";
        return `<div class="comp-seg${tiny}" style="flex: 0 0 ${pct}%; background:${BASE_COLORS[b]}" ` +
               `title="${b}: ${counts[b].toLocaleString()} (${pct.toFixed(1)}%)">${b}</div>`;
      }).join("");
    $("comp-bar").innerHTML = bar;
    $("comp-legend").innerHTML = order
      .map((b) => `<span><span class="comp-dot" style="background:${BASE_COLORS[b]}"></span>` +
                  `${b === "other" ? "N / other" : b} · ${(frac[b] * 100).toFixed(1)}%</span>`)
      .join("");
    box.classList.remove("hidden");
  }

  // ---------------------------------------------------------------- Analyzer
  function initAnalyzer() {
    const form = $("analyze-form");
    if (!form) return;
    const sampleSelect = $("sample-select");
    loadSampleList(sampleSelect);
    wireDropzone("dropzone", "file");

    async function buildFormData() {
      const fd = new FormData();
      const ok = await attachInput(fd, $("file"), sampleSelect, "file");
      if (!ok) throw new Error("Choose a file or a sample first.");
      fd.append("k", $("k").value);
      fd.append("top_n", $("top_n").value);
      fd.append("canonical", $("canonical").checked);
      fd.append("include_ambiguous", $("include_ambiguous").checked);
      fd.append("fcgr", $("fcgr").checked);
      return fd;
    }

    form.addEventListener("submit", async (e) => {
      e.preventDefault();
      clearError();
      $("results").classList.add("hidden");
      setLoading(true);
      try {
        const fd = await buildFormData();
        const res = await fetch("/api/analyze", { method: "POST", body: fd });
        const data = await res.json();
        if (!res.ok) throw new Error(data.error || "Analysis failed.");
        renderResults(data);
      } catch (err) {
        showError(err.message);
      } finally {
        setLoading(false);
      }
    });

    function renderResults(data) {
      const s = data.summary;
      renderComposition(s.base_composition);

      $("summary-cards").innerHTML = [
        card("Format", s.format.toUpperCase()),
        card("Sequences", fmt(s.total_sequences)),
        card("Total bases", fmt(s.total_bases)),
        card("Valid k-mers", fmt(s.total_kmers)),
        card("Unique k-mers", fmt(s.unique_kmers)),
        card("GC content", (s.gc_content * 100).toFixed(1) + "%"),
        card("Skipped k-mers", fmt(s.skipped_kmers)),
        card("Invalid bases", fmt(s.invalid_bases)),
      ].join("");

      const warn = $("warnings");
      if (s.warnings && s.warnings.length) {
        warn.innerHTML = "<strong>Note:</strong> " + s.warnings.join(" ");
        warn.classList.remove("hidden");
      } else warn.classList.add("hidden");

      $("chart-top").src = data.charts.top_bar;
      $("chart-hist").src = data.charts.histogram;
      if (data.charts.fcgr) {
        $("chart-fcgr").src = data.charts.fcgr;
        $("fcgr-box").classList.remove("hidden");
      } else {
        $("fcgr-box").classList.add("hidden");
      }

      const maxCount = s.top_kmers.length ? s.top_kmers[0].count : 1;
      const tbody = $("kmer-table").querySelector("tbody");
      tbody.innerHTML = s.top_kmers.map((row, i) => {
        const w = Math.max(4, (row.count / maxCount) * 90);
        return `<tr><td class="rank">${i + 1}</td>` +
          `<td class="kmer">${colorizeSeq(row.kmer)}</td>` +
          `<td class="count bar-cell">${fmt(row.count)}<span class="bar" style="width:${w}px"></span></td></tr>`;
      }).join("");

      $("results").classList.remove("hidden");
      $("results").scrollIntoView({ behavior: "smooth", block: "start" });
    }

    async function download(endpoint) {
      clearError();
      try {
        const res = await fetch(endpoint, { method: "POST", body: await buildFormData() });
        if (!res.ok) throw new Error((await res.json()).error || "Export failed.");
        triggerDownload(res);
      } catch (err) { showError(err.message); }
    }

    $("export-csv").addEventListener("click", () => download("/api/export/csv"));
    $("export-json").addEventListener("click", () => download("/api/export/json"));
  }

  // ---------------------------------------------------------------- Compare
  function initCompare() {
    const form = $("compare-form");
    if (!form) return;
    const sampleA = $("sample-a"), sampleB = $("sample-b");
    loadSampleList(sampleA);
    loadSampleList(sampleB);
    wireDropzone("dropzone-a", "file_a");
    wireDropzone("dropzone-b", "file_b");

    async function buildFormData() {
      const fd = new FormData();
      const okA = await attachInput(fd, $("file_a"), sampleA, "file_a");
      const okB = await attachInput(fd, $("file_b"), sampleB, "file_b");
      if (!okA || !okB) throw new Error("Provide both File A and File B.");
      fd.append("k", $("k-cmp").value);
      fd.append("canonical", $("canonical-cmp").checked);
      fd.append("include_ambiguous", $("ambiguous-cmp").checked);
      return fd;
    }

    form.addEventListener("submit", async (e) => {
      e.preventDefault();
      clearError();
      $("compare-results").classList.add("hidden");
      setLoading(true);
      try {
        const res = await fetch("/api/compare", { method: "POST", body: await buildFormData() });
        const data = await res.json();
        if (!res.ok) throw new Error(data.error || "Comparison failed.");
        renderCompare(data);
      } catch (err) {
        showError(err.message);
      } finally {
        setLoading(false);
      }
    });

    function renderCompare(data) {
      const c = data.comparison;
      $("compare-cards").innerHTML = [
        card("Jaccard", c.jaccard_similarity.toFixed(4), "shared sets"),
        card("Cosine", c.cosine_similarity.toFixed(4), "frequency vectors"),
        card("Shared k-mers", fmt(c.shared_kmers)),
        card("Unique to A", fmt(c.unique_to_a)),
        card("Unique to B", fmt(c.unique_to_b)),
      ].join("");

      const fa = data.file_a, fb = data.file_b;
      $("file-stats").innerHTML =
        card("File A · " + fa.format.toUpperCase(), fmt(fa.unique_kmers) + " k-mers",
             fmt(fa.sequences) + " sequences · " + fmt(fa.bases) + " bases") +
        card("File B · " + fb.format.toUpperCase(), fmt(fb.unique_kmers) + " k-mers",
             fmt(fb.sequences) + " sequences · " + fmt(fb.bases) + " bases");

      $("chart-compare").src = data.chart;
      $("compare-results").classList.remove("hidden");
      $("compare-results").scrollIntoView({ behavior: "smooth", block: "start" });
    }

    async function download(fmtType) {
      clearError();
      try {
        const fd = await buildFormData();
        fd.append("format", fmtType);
        const res = await fetch("/api/compare/export", { method: "POST", body: fd });
        if (!res.ok) throw new Error((await res.json()).error || "Export failed.");
        triggerDownload(res);
      } catch (err) { showError(err.message); }
    }

    $("export-cmp-json").addEventListener("click", () => download("json"));
    $("export-cmp-csv").addEventListener("click", () => download("csv"));
  }

  async function triggerDownload(res) {
    const blob = await res.blob();
    const disposition = res.headers.get("Content-Disposition") || "";
    const match = disposition.match(/filename=([^;]+)/);
    const filename = match ? match[1].trim() : "download";
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url; a.download = filename;
    document.body.appendChild(a); a.click(); a.remove();
    URL.revokeObjectURL(url);
  }

  return { initRibbon, initAnalyzer, initCompare };
})();
