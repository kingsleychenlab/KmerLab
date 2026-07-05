/* KmerLab front-end. Vanilla JS, no frameworks. Talks to the local Flask API. */
const KmerLab = (function () {
  "use strict";

  function $(id) { return document.getElementById(id); }

  function showError(msg) {
    const panel = $("error-panel");
    panel.textContent = msg;
    panel.classList.remove("hidden");
  }
  function clearError() { $("error-panel").classList.add("hidden"); }

  function setLoading(on) {
    $("loading").classList.toggle("hidden", !on);
  }

  function fmt(n) {
    return typeof n === "number" ? n.toLocaleString() : n;
  }

  // Populate a <select> with the bundled sample filenames.
  async function loadSampleList(selectEl) {
    try {
      const res = await fetch("/api/samples");
      const files = await res.json();
      files.forEach((f) => {
        const opt = document.createElement("option");
        opt.value = f;
        opt.textContent = f;
        selectEl.appendChild(opt);
      });
    } catch (e) {
      /* samples are optional; ignore */
    }
  }

  // Attach either the uploaded file or the selected sample's text to formData.
  async function attachInput(formData, fileEl, sampleEl, field) {
    if (fileEl.files && fileEl.files.length > 0) {
      formData.append(field, fileEl.files[0]);
      return true;
    }
    if (sampleEl && sampleEl.value) {
      const res = await fetch("/samples/" + encodeURIComponent(sampleEl.value));
      if (!res.ok) throw new Error("Could not load sample: " + sampleEl.value);
      const text = await res.text();
      formData.append(field + "_text", text);
      return true;
    }
    return false;
  }

  function card(label, value) {
    return `<div class="card"><div class="card-value">${value}</div><div class="card-label">${label}</div></div>`;
  }

  // ---------------------------------------------------------------- Analyzer
  function initAnalyzer() {
    const form = $("analyze-form");
    const sampleSelect = $("sample-select");
    loadSampleList(sampleSelect);

    // Remember the last-built FormData inputs for exports.
    async function buildFormData() {
      const fd = new FormData();
      const ok = await attachInput(fd, $("file"), sampleSelect, "file");
      if (!ok) throw new Error("Please choose a file or a sample first.");
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
      const cards = [
        card("Format", s.format.toUpperCase()),
        card("Sequences", fmt(s.total_sequences)),
        card("Total bases", fmt(s.total_bases)),
        card("Valid k-mers", fmt(s.total_kmers)),
        card("Unique k-mers", fmt(s.unique_kmers)),
        card("GC content", (s.gc_content * 100).toFixed(2) + "%"),
        card("Skipped k-mers", fmt(s.skipped_kmers)),
        card("Invalid bases", fmt(s.invalid_bases)),
      ].join("");
      $("summary-cards").innerHTML = cards;

      const warn = $("warnings");
      if (s.warnings && s.warnings.length) {
        warn.innerHTML = "<strong>Warnings:</strong> " + s.warnings.join(" ");
        warn.classList.remove("hidden");
      } else {
        warn.classList.add("hidden");
      }

      $("chart-top").src = data.charts.top_bar;
      $("chart-hist").src = data.charts.histogram;
      if (data.charts.fcgr) {
        $("chart-fcgr").src = data.charts.fcgr;
        $("fcgr-box").classList.remove("hidden");
      } else {
        $("fcgr-box").classList.add("hidden");
      }

      const tbody = $("kmer-table").querySelector("tbody");
      tbody.innerHTML = s.top_kmers
        .map((row, i) => `<tr><td>${i + 1}</td><td>${row.kmer}</td><td>${fmt(row.count)}</td></tr>`)
        .join("");

      $("results").classList.remove("hidden");
    }

    async function download(endpoint) {
      clearError();
      try {
        const fd = await buildFormData();
        const res = await fetch(endpoint, { method: "POST", body: fd });
        if (!res.ok) {
          const data = await res.json();
          throw new Error(data.error || "Export failed.");
        }
        triggerDownload(res);
      } catch (err) {
        showError(err.message);
      }
    }

    $("export-csv").addEventListener("click", () => download("/api/export/csv"));
    $("export-json").addEventListener("click", () => download("/api/export/json"));
  }

  // ---------------------------------------------------------------- Compare
  function initCompare() {
    const form = $("compare-form");
    const sampleA = $("sample-a");
    const sampleB = $("sample-b");
    loadSampleList(sampleA);
    loadSampleList(sampleB);

    async function buildFormData() {
      const fd = new FormData();
      const okA = await attachInput(fd, $("file_a"), sampleA, "file_a");
      const okB = await attachInput(fd, $("file_b"), sampleB, "file_b");
      if (!okA || !okB) throw new Error("Please provide both File A and File B.");
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
        const fd = await buildFormData();
        const res = await fetch("/api/compare", { method: "POST", body: fd });
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
        card("Jaccard similarity", c.jaccard_similarity.toFixed(4)),
        card("Cosine similarity", c.cosine_similarity.toFixed(4)),
        card("Shared k-mers", fmt(c.shared_kmers)),
        card("Unique to A", fmt(c.unique_to_a)),
        card("Unique to B", fmt(c.unique_to_b)),
      ].join("");

      const fa = data.file_a, fb = data.file_b;
      $("file-stats").innerHTML =
        `<div class="card"><div class="card-label">File A (${fa.format.toUpperCase()})</div>
           <p>${fmt(fa.sequences)} sequences · ${fmt(fa.bases)} bases · ${fmt(fa.unique_kmers)} unique k-mers</p></div>
         <div class="card"><div class="card-label">File B (${fb.format.toUpperCase()})</div>
           <p>${fmt(fb.sequences)} sequences · ${fmt(fb.bases)} bases · ${fmt(fb.unique_kmers)} unique k-mers</p></div>`;

      $("chart-compare").src = data.chart;
      $("compare-results").classList.remove("hidden");
    }

    async function download(fmtType) {
      clearError();
      try {
        const fd = await buildFormData();
        fd.append("format", fmtType);
        const res = await fetch("/api/compare/export", { method: "POST", body: fd });
        if (!res.ok) {
          const data = await res.json();
          throw new Error(data.error || "Export failed.");
        }
        triggerDownload(res);
      } catch (err) {
        showError(err.message);
      }
    }

    $("export-cmp-json").addEventListener("click", () => download("json"));
    $("export-cmp-csv").addEventListener("click", () => download("csv"));
  }

  // Trigger a browser download from a fetch Response with Content-Disposition.
  async function triggerDownload(res) {
    const blob = await res.blob();
    const disposition = res.headers.get("Content-Disposition") || "";
    const match = disposition.match(/filename=([^;]+)/);
    const filename = match ? match[1].trim() : "download";
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  }

  return { initAnalyzer, initCompare };
})();
