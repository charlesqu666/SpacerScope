# Copyright (C) 2026 charlesqu666
# SPDX-License-Identifier: AGPL-3.0
#pragma once

#include <string>

namespace web_assets {
inline const std::string& index_html() {
    static const std::string html = [] {
        std::string value;
        value += R"HTML(
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>SpacerScope Web</title>
  <style>
    :root {
      --bg: #f4f7fb;
      --panel: #ffffff;
      --line: #dbe5f3;
      --text: #1f2937;
      --muted: #64748b;
      --brand: #0f766e;
      --brand-2: #1d4ed8;
      --danger: #dc2626;
      --warn: #b45309;
      --ok: #059669;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      font-family: "Segoe UI", Arial, sans-serif;
      background: linear-gradient(180deg, #edf4ff 0%, var(--bg) 100%);
      color: var(--text);
    }
    .page {
      width: min(1360px, 94vw);
      margin: 24px auto;
      background: var(--panel);
      border-radius: 18px;
      box-shadow: 0 12px 36px rgba(15, 23, 42, 0.08);
      overflow: hidden;
    }
    header {
      padding: 22px 28px;
      background: linear-gradient(135deg, var(--brand-2), var(--brand));
      color: white;
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 18px;
      flex-wrap: wrap;
    }
    header h1 {
      margin: 0;
      font-size: 28px;
    }
    header p {
      margin: 6px 0 0;
      opacity: 0.92;
      font-size: 14px;
    }
    .badge {
      display: inline-block;
      padding: 5px 10px;
      border-radius: 999px;
      font-size: 12px;
      font-weight: 700;
      background: rgba(255,255,255,0.16);
      border: 1px solid rgba(255,255,255,0.22);
    }
    .banner {
      padding: 12px 28px;
      background: #fff7ed;
      border-bottom: 1px solid #fed7aa;
      color: var(--warn);
      font-size: 13px;
    }
    .grid {
      display: grid;
      grid-template-columns: 420px 1fr;
      gap: 0;
      min-height: 760px;
    }
    .sidebar {
      border-right: 1px solid var(--line);
      padding: 20px;
      background: #fbfdff;
    }
    .content {
      padding: 20px;
      display: grid;
      grid-template-rows: auto auto 1fr auto;
      gap: 16px;
    }
    .card {
      background: white;
      border: 1px solid var(--line);
      border-radius: 14px;
      padding: 16px;
    }
    .card h2 {
      margin: 0 0 12px;
      font-size: 17px;
    }
    .row {
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 12px;
      margin-bottom: 12px;
    }
    .row.single { grid-template-columns: 1fr; }
    label {
      display: block;
      font-size: 13px;
      color: var(--muted);
      margin-bottom: 6px;
      font-weight: 600;
    }
    input[type="text"],
    input[type="number"],
    input[type="file"],
    select,
    textarea {
      width: 100%;
      padding: 10px 12px;
      border: 1px solid #cbd5e1;
      border-radius: 10px;
      font: inherit;
      background: white;
      color: var(--text);
    }
    textarea {
      min-height: 86px;
      resize: vertical;
    }
    .toggle {
      display: inline-flex;
      background: #eef4ff;
      border-radius: 999px;
      padding: 4px;
      gap: 4px;
      margin-bottom: 14px;
    }
    .toggle button {
      border: none;
      background: transparent;
      padding: 8px 14px;
      border-radius: 999px;
      cursor: pointer;
      font-weight: 600;
      color: var(--muted);
    }
    .toggle button.active {
      background: white;
      color: var(--brand-2);
      box-shadow: 0 2px 8px rgba(15, 23, 42, 0.08);
    }
    .actions {
      display: flex;
      gap: 10px;
      flex-wrap: wrap;
      margin-top: 12px;
    }
    .btn {
      border: none;
      border-radius: 10px;
      padding: 10px 16px;
      font: inherit;
      font-weight: 700;
      cursor: pointer;
      color: white;
      background: var(--brand-2);
    }
    .btn.secondary { background: #334155; }
    .btn.danger { background: var(--danger); }
    .btn:disabled {
      opacity: 0.55;
      cursor: not-allowed;
    }
    .kv {
      display: grid;
      grid-template-columns: 180px 1fr;
      gap: 8px 14px;
      font-size: 14px;
    }
    .kv .k {
      color: var(--muted);
      font-weight: 600;
    }
    .status {
      display: inline-block;
      padding: 4px 10px;
      border-radius: 999px;
      font-size: 12px;
      font-weight: 700;
      text-transform: uppercase;
      letter-spacing: 0.04em;
      background: #e2e8f0;
      color: #334155;
    }
    .status.running { background: #dbeafe; color: #1d4ed8; }
    .status.completed { background: #dcfce7; color: var(--ok); }
    .status.failed { background: #fee2e2; color: var(--danger); }
    .status.canceled { background: #fef3c7; color: var(--warn); }
    .metric-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 12px;
    }
    .metric {
      border: 1px solid var(--line);
      border-radius: 12px;
      padding: 12px;
      background: #f8fbff;
    }
    .metric .label {
      font-size: 12px;
      color: var(--muted);
      text-transform: uppercase;
      letter-spacing: 0.05em;
    }
    .metric .value {
      margin-top: 8px;
      font-size: 20px;
      font-weight: 800;
    }
    .progress-wrap {
      border: 1px solid var(--line);
      border-radius: 12px;
      padding: 12px;
      background: #fbfdff;
    }
    progress {
      width: 100%;
      height: 18px;
    }
    pre {
      margin: 0;
      padding: 14px;
      min-height: 320px;
      max-height: 420px;
      overflow: auto;
      border-radius: 12px;
      background: #0f172a;
      color: #e2e8f0;
      font-family: Consolas, Monaco, monospace;
      font-size: 12px;
      line-height: 1.5;
      white-space: pre-wrap;
      word-break: break-word;
    }
    .artifact-list {
      display: flex;
      gap: 10px;
      flex-wrap: wrap;
    }
    .artifact-list a {
      text-decoration: none;
      color: var(--brand-2);
      font-weight: 700;
    }
    .hidden { display: none; }
    .small {
      font-size: 12px;
      color: var(--muted);
      margin-top: 6px;
    }
    @media (max-width: 1080px) {
      .grid { grid-template-columns: 1fr; }
      .sidebar { border-right: none; border-bottom: 1px solid var(--line); }
    }
  </style>
</head>
<body>
  <div class="page">
    <header>
      <div>
        <h1 id="appTitle">SpacerScope Web</h1>
      </div>
      <div style="display:flex;align-items:center;gap:10px;flex-wrap:wrap;">
        <label for="langSelect" style="margin:0;color:white;font-size:12px;font-weight:700;" id="langLabel">Language</label>
        <select id="langSelect" style="width:auto;min-width:96px;padding:8px 10px;border:none;border-radius:999px;background:rgba(255,255,255,0.18);color:white;font-weight:700;">
          <option value="en" style="color:#111827;">English</option>
          <option value="zh" style="color:#111827;">中文</option>
        </select>
        <span class="badge" id="serverBadge">loading...</span>
      </div>
    </header>
    <div class="banner" id="lanWarning">
      LAN mode in v1 has no authentication. Only bind to non-local addresses on a trusted network.
    </div>
    <div class="grid">
      <aside class="sidebar">
        <div class="card">
          <h2 id="runTaskTitle">Run Task</h2>
          <div class="toggle" id="inputToggle">
            <button type="button" data-mode="path" class="active" id="inputModePath">Server Paths</button>
            <button type="button" data-mode="upload" id="inputModeUpload">Upload Files</button>
          </div>
          <form id="runForm">
            <div class="row">
              <div>
                <label for="mode" id="modeLabel">Mode</label>
                <select id="mode" name="mode">
                  <option value="fasta">fasta</option>
                  <option value="fasta-cut">fasta-cut</option>
                  <option value="anno">anno</option>
                  <option value="anno-cut">anno-cut</option>
                </select>
              </div>
              <div>
                <label for="threads" id="threadsLabel">Threads</label>
                <input id="threads" name="threads" type="number" min="1" value="8">
              </div>
            </div>

            <div id="pathFields">
              <div class="row single">
                <div>
                  <label for="refPath" id="refPathLabel">Reference FASTA path</label>
                  <input id="refPath" name="ref_path" type="text" placeholder="/data/ref.fa">
                </div>
              </div>
              <div class="row single query-only">
                <div>
                  <label for="queryPath" id="queryPathLabel">Query FASTA path</label>
                  <input id="queryPath" name="query_path" type="text" placeholder="/data/query.fa">
                </div>
              </div>
              <div class="row single anno-only hidden">
                <div>
                  <label for="annotationPath" id="annotationPathLabel">Annotation GFF/GFF3 path</label>
                  <input id="annotationPath" name="annotation_path" type="text" placeholder="/data/genes.gff3">
                </div>
              </div>
            </div>

            <div id="uploadFields" class="hidden">
              <div class="row single">
                <div>
                  <label for="refUpload" id="refUploadLabel">Reference FASTA upload</label>
                  <input id="refUpload" name="ref_upload" type="file" accept=".fa,.fasta,.fna,.fa.gz,.fasta.gz,.fna.gz">
                </div>
              </div>
              <div class="row single query-only">
                <div>
                  <label for="queryUpload" id="queryUploadLabel">Query FASTA upload</label>
                  <input id="queryUpload" name="query_upload" type="file" accept=".fa,.fasta,.fna,.fa.gz,.fasta.gz,.fna.gz">
                </div>
              </div>
              <div class="row single anno-only hidden">
                <div>
                  <label for="annotationUpload" id="annotationUploadLabel">Annotation upload</label>
                  <input id="annotationUpload" name="annotation_upload" type="file" accept=".gff,.gff3,.gtf">
                </div>
              </div>
            </div>

            <div class="row anno-only hidden">
              <div>
                <label for="geneId" id="geneIdLabel">Gene ID</label>
                <input id="geneId" name="gene_id" type="text" placeholder="gene:Gene1">
              </div>
              <div>
                <label for="seqType" id="seqTypeLabel">Sequence type</label>
                <select id="seqType" name="seq_type">
                  <option value="gene">gene</option>
                  <option value="mrna">mrna</option>
                  <option value="cds">cds</option>
                </select>
              </div>
            </div>

            <div class="row">
              <div>
                <label for="pam" id="pamLabel">PAM</label>
                <input id="pam" name="pam" type="text" value="NGG">
              </div>
              <div>
                <label for="altPam" id="altPamLabel">Alternate PAM</label>
                <input id="altPam" name="alt_pam" type="text" value="NAG">
              </div>
            </div>
            <div class="row">
              <div>
                <label for="spacerLen" id="spacerLenLabel">Spacer length</label>
                <input id="spacerLen" name="spacer_len" type="number" min="1" value="20">
              </div>
              <div>
                <label for="direction" id="directionLabel">Direction</label>
                <select id="direction" name="direction">
                  <option value="up">up</option>
                  <option value="down">down</option>
                </select>
              </div>
            </div>
            <div class="row">
              <div>
                <label for="mismatch" id="mismatchLabel">Max mismatches</label>
                <input id="mismatch" name="mismatch" type="number" min="0" value="4">
              </div>
              <div>
)HTML";
        value += R"HTML(
                <label for="indel" id="indelLabel">Max indels</label>
                <input id="indel" name="indel" type="number" min="0" value="2">
              </div>
            </div>
            <div class="row">
              <div>
                <label for="dust" id="dustLabel">DUST threshold</label>
                <input id="dust" name="dust_threshold" type="number" min="0" value="12">
              </div>
              <div>
                <label>&nbsp;</label>
                <div>
                  <label style="display:flex;gap:8px;align-items:center;margin:0;" id="rawOutputText">
                    <input id="rawOutput" name="raw_output" type="checkbox">
                    Raw TSV only
                  </label>
                  <label style="display:flex;gap:8px;align-items:center;margin:8px 0 0;" id="generateHtmlText">
                    <input id="generateHtml" name="generate_html" type="checkbox" checked>
                    Generate HTML
                  </label>
                </div>
              </div>
            </div>
            <div class="actions">
              <button class="btn" id="runBtn" type="submit">Run</button>
              <button class="btn danger" id="cancelBtn" type="button">Cancel</button>
              <button class="btn secondary" id="refreshBtn" type="button">Refresh Status</button>
            </div>
            <div class="small" id="outputHint">
              Output files are stored in the current web task workspace and exposed through the artifact links.
            </div>
          </form>
        </div>
      </aside>

      <main class="content">
        <section class="card">
          <h2 id="taskStatusTitle">Task Status</h2>
          <div class="kv">
            <div class="k" id="stateKey">State</div><div><span id="statusBadge" class="status">idle</span></div>
            <div class="k" id="modeKey">Mode</div><div id="statusMode">-</div>
            <div class="k" id="workspaceKey">Workspace</div><div id="statusWorkspace">-</div>
            <div class="k" id="lastErrorKey">Last error</div><div id="statusError">-</div>
          </div>
        </section>

        <section class="card">
          <h2 id="memoryEstimateTitle">Memory Estimate</h2>
          <div class="metric-grid">
            <div class="metric"><div class="label" id="memChromsLabel">Chromosomes</div><div class="value" id="memChroms">-</div></div>
            <div class="metric"><div class="label" id="memTotalLabel">Total Genome</div><div class="value" id="memTotal">-</div></div>
            <div class="metric"><div class="label" id="memLongestLabel">Longest Chromosome</div><div class="value" id="memLongest">-</div></div>
            <div class="metric"><div class="label" id="memPeakLabel">Peak Estimate</div><div class="value" id="memPeak">-</div></div>
          </div>
          <div class="small" id="memModel">No estimate yet.</div>
        </section>

        <section class="card">
          <h2 id="progressTitle">Progress</h2>
          <div class="progress-wrap">
            <div id="progressLabel">Waiting for task.</div>
            <progress id="progressBar" max="100" value="0"></progress>
          </div>
        </section>

        <section class="card">
          <h2 id="logTitle">Log</h2>
          <pre id="logBox"></pre>
        </section>

        <section class="card">
          <h2 id="artifactsTitle">Artifacts</h2>
          <div class="artifact-list" id="artifactList">
            <span class="small">No artifacts yet.</span>
          </div>
        </section>
      </main>
    </div>
  </div>

  <script>
    const state = {
      inputMode: 'path',
      currentRunId: null,
      eventSource: null,
      lang: 'en',
    };

    const I18N = {
      en: {
        appTitle: 'SpacerScope Web',
        langLabel: 'Language',
        serverBadgeLoading: 'loading...',
        lanWarning: 'LAN mode in v1 has no authentication. Only bind to non-local addresses on a trusted network.',
        runTaskTitle: 'Run Task',
        inputModePath: 'Server Paths',
        inputModeUpload: 'Upload Files',
        modeLabel: 'Mode',
        threadsLabel: 'Threads',
        refPathLabel: 'Reference FASTA path',
        queryPathLabel: 'Query FASTA path',
        annotationPathLabel: 'Annotation GFF/GFF3 path',
        refUploadLabel: 'Reference FASTA upload',
        queryUploadLabel: 'Query FASTA upload',
        annotationUploadLabel: 'Annotation upload',
        geneIdLabel: 'Gene ID',
        seqTypeLabel: 'Sequence type',
        pamLabel: 'PAM',
        altPamLabel: 'Alternate PAM',
        spacerLenLabel: 'Spacer length',
        directionLabel: 'Direction',
        mismatchLabel: 'Max mismatches',
        indelLabel: 'Max indels',
        dustLabel: 'DUST threshold',
        rawOutputText: 'Raw TSV only',
        generateHtmlText: 'Generate HTML',
        outputHint: 'Output files are stored in the current web task workspace and exposed through the artifact links.',
        taskStatusTitle: 'Task Status',
        stateKey: 'State',
        modeKey: 'Mode',
        workspaceKey: 'Workspace',
        lastErrorKey: 'Last error',
        memoryEstimateTitle: 'Memory Estimate',
        memChromsLabel: 'Chromosomes',
        memTotalLabel: 'Total Genome',
        memLongestLabel: 'Longest Chromosome',
        memPeakLabel: 'Peak Estimate',
        memModelEmpty: 'No estimate yet.',
        progressTitle: 'Progress',
        progressWaiting: 'Waiting for task.',
        logTitle: 'Log',
        artifactsTitle: 'Artifacts',
        noArtifacts: 'No artifacts yet.',
        runBtn: 'Run',
        cancelBtn: 'Cancel',
        refreshBtn: 'Refresh Status',
        webTaskAccepted: '[web] task accepted',
        webCancelRequested: '[web] cancel requested',
        webEventReconnect: '[web] event stream disconnected; retrying...',
        runRequestFailed: 'Run request failed.',
        cancelFailed: 'Cancel failed.',
        completed: 'Completed',
        memoryPeak: (gib) => `[memory] peak ~ ${gib.toFixed(2)} GiB`,
        serverBadgeFormat: (payload) => `${payload.product} ${payload.version} @ ${payload.bind}:${payload.port}`,
      },
      zh: {
        appTitle: 'SpacerScope 网页版',
        langLabel: '语言',
        serverBadgeLoading: '加载中...',
        lanWarning: 'v1 局域网模式不提供认证。只有在可信网络中才应绑定到非本机地址。',
        runTaskTitle: '运行任务',
        inputModePath: '服务器路径',
        inputModeUpload: '上传文件',
        modeLabel: '模式',
        threadsLabel: '线程数',
        refPathLabel: '参考基因组 FASTA 路径',
        queryPathLabel: '查询 FASTA 路径',
        annotationPathLabel: '注释 GFF/GFF3 路径',
        refUploadLabel: '上传参考基因组 FASTA',
        queryUploadLabel: '上传查询 FASTA',
        annotationUploadLabel: '上传注释文件',
        geneIdLabel: '基因 ID',
        seqTypeLabel: '序列类型',
        pamLabel: 'PAM',
        altPamLabel: '备用 PAM',
        spacerLenLabel: 'spacer 长度',
        directionLabel: '方向',
        mismatchLabel: '最大错配数',
        indelLabel: '最大 indel 数',
        dustLabel: 'DUST 阈值',
        rawOutputText: '仅输出原始 TSV',
        generateHtmlText: '生成 HTML 报告',
        outputHint: '输出文件会保存在当前 Web 任务工作目录中，并通过下方产物链接提供下载。',
        taskStatusTitle: '任务状态',
        stateKey: '状态',
        modeKey: '模式',
        workspaceKey: '工作目录',
        lastErrorKey: '最近错误',
        memoryEstimateTitle: '内存估算',
        memChromsLabel: '染色体数',
        memTotalLabel: '基因组总长度',
        memLongestLabel: '最长染色体',
        memPeakLabel: '峰值估算',
        memModelEmpty: '暂未生成估算。',
        progressTitle: '进度',
        progressWaiting: '等待任务开始。',
        logTitle: '日志',
        artifactsTitle: '结果产物',
        noArtifacts: '暂无产物。',
        runBtn: '运行',
        cancelBtn: '取消',
        refreshBtn: '刷新状态',
        webTaskAccepted: '[web] 已接受任务',
        webCancelRequested: '[web] 已请求取消任务',
        webEventReconnect: '[web] 事件流已断开，正在重连...',
        runRequestFailed: '任务提交失败。',
        cancelFailed: '取消失败。',
        completed: '已完成',
        memoryPeak: (gib) => `[memory] 峰值约 ${gib.toFixed(2)} GiB`,
        serverBadgeFormat: (payload) => `${payload.product} ${payload.version} @ ${payload.bind}:${payload.port}`,
      }
    };

    const logBox = document.getElementById('logBox');
    const statusBadge = document.getElementById('statusBadge');
    const statusMode = document.getElementById('statusMode');
    const statusWorkspace = document.getElementById('statusWorkspace');
    const statusError = document.getElementById('statusError');
    const artifactList = document.getElementById('artifactList');
    const progressBar = document.getElementById('progressBar');
    const progressLabel = document.getElementById('progressLabel');
    const memChroms = document.getElementById('memChroms');
    const memTotal = document.getElementById('memTotal');
    const memLongest = document.getElementById('memLongest');
    const memPeak = document.getElementById('memPeak');
    const memModel = document.getElementById('memModel');

    function t(key) {
      return (I18N[state.lang] || I18N.en)[key];
    }

    function applyLanguage() {
      document.documentElement.lang = state.lang === 'zh' ? 'zh-CN' : 'en';
      const langSelect = document.getElementById('langSelect');
      langSelect.value = state.lang;
      if (langSelect.options.length >= 2) {
        langSelect.options[0].text = 'English';
        langSelect.options[1].text = '中文';
      }
      document.getElementById('appTitle').textContent = t('appTitle');
      document.getElementById('langLabel').textContent = t('langLabel');
      document.getElementById('lanWarning').textContent = t('lanWarning');
      document.getElementById('runTaskTitle').textContent = t('runTaskTitle');
      document.getElementById('inputModePath').textContent = t('inputModePath');
      document.getElementById('inputModeUpload').textContent = t('inputModeUpload');
      document.getElementById('modeLabel').textContent = t('modeLabel');
      document.getElementById('threadsLabel').textContent = t('threadsLabel');
      document.getElementById('refPathLabel').textContent = t('refPathLabel');
      document.getElementById('queryPathLabel').textContent = t('queryPathLabel');
      document.getElementById('annotationPathLabel').textContent = t('annotationPathLabel');
      document.getElementById('refUploadLabel').textContent = t('refUploadLabel');
      document.getElementById('queryUploadLabel').textContent = t('queryUploadLabel');
      document.getElementById('annotationUploadLabel').textContent = t('annotationUploadLabel');
      document.getElementById('geneIdLabel').textContent = t('geneIdLabel');
      document.getElementById('seqTypeLabel').textContent = t('seqTypeLabel');
      document.getElementById('pamLabel').textContent = t('pamLabel');
      document.getElementById('altPamLabel').textContent = t('altPamLabel');
      document.getElementById('spacerLenLabel').textContent = t('spacerLenLabel');
      document.getElementById('directionLabel').textContent = t('directionLabel');
      document.getElementById('mismatchLabel').textContent = t('mismatchLabel');
      document.getElementById('indelLabel').textContent = t('indelLabel');
      document.getElementById('dustLabel').textContent = t('dustLabel');
      document.getElementById('outputHint').textContent = t('outputHint');
      document.getElementById('taskStatusTitle').textContent = t('taskStatusTitle');
      document.getElementById('stateKey').textContent = t('stateKey');
      document.getElementById('modeKey').textContent = t('modeKey');
      document.getElementById('workspaceKey').textContent = t('workspaceKey');
      document.getElementById('lastErrorKey').textContent = t('lastErrorKey');
      document.getElementById('memoryEstimateTitle').textContent = t('memoryEstimateTitle');
      document.getElementById('memChromsLabel').textContent = t('memChromsLabel');
      document.getElementById('memTotalLabel').textContent = t('memTotalLabel');
)HTML";
        value += R"HTML(
      document.getElementById('memLongestLabel').textContent = t('memLongestLabel');
      document.getElementById('memPeakLabel').textContent = t('memPeakLabel');
      document.getElementById('progressTitle').textContent = t('progressTitle');
      document.getElementById('logTitle').textContent = t('logTitle');
      document.getElementById('artifactsTitle').textContent = t('artifactsTitle');
      document.getElementById('runBtn').textContent = t('runBtn');
      document.getElementById('cancelBtn').textContent = t('cancelBtn');
      document.getElementById('refreshBtn').textContent = t('refreshBtn');
      document.getElementById('serverBadge').textContent =
        document.getElementById('serverBadge').dataset.ready === '1'
          ? document.getElementById('serverBadge').textContent
          : t('serverBadgeLoading');
      if (memModel.dataset.ready !== '1') memModel.textContent = t('memModelEmpty');
      if (progressLabel.dataset.ready !== '1') progressLabel.textContent = t('progressWaiting');
      if (artifactList.dataset.ready !== '1') artifactList.innerHTML = `<span class="small">${t('noArtifacts')}</span>`;
      const rawText = document.getElementById('rawOutputText');
      rawText.lastChild.textContent = ' ' + t('rawOutputText');
      const htmlText = document.getElementById('generateHtmlText');
      htmlText.lastChild.textContent = ' ' + t('generateHtmlText');
      try { localStorage.setItem('spacerscope-web-lang', state.lang); } catch (_) {}
    }

    function appendLog(line) {
      logBox.textContent += line + "\n";
      logBox.scrollTop = logBox.scrollHeight;
    }

    function setStatus(payload) {
      const stateName = payload.state || 'idle';
      statusBadge.textContent = stateName;
      statusBadge.className = 'status ' + stateName;
      statusMode.textContent = payload.mode || '-';
      statusWorkspace.textContent = payload.workspace || '-';
      statusError.textContent = payload.last_error || '-';

      const artifacts = payload.artifacts || {};
      const links = [];
      ['tsv', 'html', 'stdout', 'stderr'].forEach((kind) => {
        if (artifacts[kind]) {
          links.push(`<a href="/api/artifacts/${kind}" target="_blank">${kind}</a>`);
        }
      });
      artifactList.dataset.ready = links.length ? '1' : '0';
      artifactList.innerHTML = links.length ? links.join('') : `<span class="small">${t('noArtifacts')}</span>`;
    }

    function setMemoryEstimate(payload) {
      memChroms.textContent = payload.chromosomes ?? '-';
      memTotal.textContent = payload.total_genome_mbp != null ? payload.total_genome_mbp.toFixed(2) + ' Mb' : '-';
      memLongest.textContent = payload.longest_chromosome_mbp != null ? payload.longest_chromosome_mbp.toFixed(2) + ' Mb' : '-';
      memPeak.textContent = payload.peak_gib != null ? payload.peak_gib.toFixed(2) + ' GiB' : '-';
      memModel.dataset.ready = '1';
      memModel.textContent = payload.model || t('memModelEmpty');
    }

    function setProgress(label, percent) {
      progressLabel.dataset.ready = '1';
      progressLabel.textContent = label;
      progressBar.value = percent;
    }

    function updateModeFields() {
      const mode = document.getElementById('mode').value;
      const isAnno = mode === 'anno' || mode === 'anno-cut';
      document.querySelectorAll('.anno-only').forEach((el) => el.classList.toggle('hidden', !isAnno));
      document.querySelectorAll('.query-only').forEach((el) => el.classList.toggle('hidden', isAnno));
    }

    function updateInputMode() {
      const isUpload = state.inputMode === 'upload';
      document.getElementById('pathFields').classList.toggle('hidden', isUpload);
      document.getElementById('uploadFields').classList.toggle('hidden', !isUpload);
      document.querySelectorAll('#pathFields input').forEach((el) => { el.disabled = isUpload; });
      document.querySelectorAll('#uploadFields input').forEach((el) => { el.disabled = !isUpload; });
      document.querySelectorAll('#inputToggle button').forEach((btn) => {
        btn.classList.toggle('active', btn.dataset.mode === state.inputMode);
      });
    }

    function bindEvents() {
      document.getElementById('mode').addEventListener('change', updateModeFields);
      document.querySelectorAll('#inputToggle button').forEach((btn) => {
        btn.addEventListener('click', () => {
          state.inputMode = btn.dataset.mode;
          updateInputMode();
        });
      });
      document.getElementById('rawOutput').addEventListener('change', (evt) => {
        if (evt.target.checked) {
          document.getElementById('generateHtml').checked = false;
        }
      });
      document.getElementById('generateHtml').addEventListener('change', (evt) => {
        if (evt.target.checked) {
          document.getElementById('rawOutput').checked = false;
        }
      });
      document.getElementById('langSelect').addEventListener('change', (evt) => {
        state.lang = evt.target.value === 'zh' ? 'zh' : 'en';
        applyLanguage();
      });
      document.getElementById('refreshBtn').addEventListener('click', refreshStatus);
      document.getElementById('cancelBtn').addEventListener('click', cancelTask);
      document.getElementById('runForm').addEventListener('submit', runTask);
    }

    async function refreshInfo() {
      const response = await fetch('/api/info');
      const payload = await response.json();
      document.getElementById('serverBadge').dataset.ready = '1';
      document.getElementById('serverBadge').textContent = t('serverBadgeFormat')(payload);
      document.getElementById('threads').value = payload.defaults.threads;
    }

    async function refreshStatus() {
      const response = await fetch('/api/status');
      const payload = await response.json();
      setStatus(payload);
      if (payload.memory_estimate) setMemoryEstimate(payload.memory_estimate);
    }

    function ensureEventSource() {
      if (state.eventSource) return;
      const evt = new EventSource('/api/events');
      evt.onmessage = (message) => {
        const payload = JSON.parse(message.data);
        handleEvent(payload);
      };
      evt.onerror = () => {
        appendLog(t('webEventReconnect'));
      };
      state.eventSource = evt;
    }

    function handleEvent(payload) {
      if (payload.event === 'status') {
        setStatus(payload);
        return;
      }
      if (payload.event === 'memory_estimate') {
        setMemoryEstimate(payload);
        appendLog(t('memoryPeak')(payload.peak_gib));
        return;
      }
      if (payload.event === 'step_start') {
        appendLog(`[step] ${payload.step}: ${payload.title}`);
        setProgress(`${payload.step}: ${payload.title}`, progressBar.value);
        return;
      }
      if (payload.event === 'step_done') {
        appendLog(`[done] ${payload.step} (${payload.seconds.toFixed(3)}s)`);
        return;
      }
      if (payload.event === 'progress') {
        setProgress(`${payload.step}: ${payload.done}/${payload.total} ${payload.unit}`, payload.percent || 0);
        appendLog(`[progress] ${payload.step} ${payload.done}/${payload.total} ${payload.unit} (${payload.percent}%)`);
        return;
      }
      if (payload.event === 'stderr') {
        appendLog(`[stderr] ${payload.message}`);
        return;
      }
      if (payload.event === 'error') {
        appendLog(`[error] ${payload.message}`);
        return;
      }
      if (payload.event === 'complete') {
        appendLog(`[complete] total ${payload.total_seconds.toFixed(3)}s`);
        setProgress(t('completed'), 100);
        refreshStatus();
        return;
      }
      if (payload.event === 'log') {
        appendLog(payload.message);
        return;
      }
    }

    async function runTask(evt) {
      evt.preventDefault();
      logBox.textContent = '';
      setProgress(state.lang === 'zh' ? '正在提交任务...' : 'Submitting task...', 0);
      const form = evt.target;
      const data = new FormData(form);
      data.set('input_mode', state.inputMode);
      data.set('raw_output', document.getElementById('rawOutput').checked ? '1' : '0');
      data.set('generate_html', document.getElementById('generateHtml').checked ? '1' : '0');

      const response = await fetch('/api/run', {
        method: 'POST',
        body: data
      });
      const payload = await response.json();
      if (!response.ok) {
        appendLog('[error] ' + (payload.error || t('runRequestFailed')));
        return;
      }
      state.currentRunId = payload.run_id;
      appendLog(t('webTaskAccepted'));
      await refreshStatus();
    }

    async function cancelTask() {
      const response = await fetch('/api/cancel', { method: 'POST' });
      const payload = await response.json();
      if (!response.ok) {
        appendLog('[error] ' + (payload.error || t('cancelFailed')));
        return;
      }
      appendLog(t('webCancelRequested'));
      await refreshStatus();
    }

    try {
      const savedLang = localStorage.getItem('spacerscope-web-lang');
      if (savedLang === 'zh' || savedLang === 'en') state.lang = savedLang;
    } catch (_) {}
    bindEvents();
    applyLanguage();
    updateModeFields();
    updateInputMode();
    refreshInfo().catch((err) => appendLog('[error] ' + err.message));
    refreshStatus().catch((err) => appendLog('[error] ' + err.message));
    ensureEventSource();
  </script>
</body>
</html>
)HTML";
        return value;
    }();
    return html;
}

}  // namespace web_assets

