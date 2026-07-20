import { useState, useEffect } from 'react';
import { motion } from 'framer-motion';
import NetworkGraph from '../components/NetworkGraph';
import {
  fetchInitialNetwork,
  fetchRewiring,
  fetchEagatLeaderboard,
  fetchEagatGene,
} from '../lib/api';

const HubRewiring = () => {
  const [geneInput, setGeneInput] = useState('');
  const [selectedGene, setSelectedGene] = useState(null);
  const [tissue, setTissue] = useState('normal');
  const [networkData, setNetworkData] = useState(null);
  const [initialNetwork, setInitialNetwork] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [topGenes, setTopGenes] = useState([]);
  const [networkSource, setNetworkSource] = useState('');
  const [hubMetrics, setHubMetrics] = useState(null);

  useEffect(() => {
    fetchInitialNetworkData();
    fetchTopGenes();
  }, []);

  const fetchTopGenes = async () => {
    try {
      const data = await fetchEagatLeaderboard(15);
      if (data?.genes?.length) setTopGenes(data.genes);
    } catch (err) {
      console.warn('EA-GAT leaderboard not available:', err.message);
    }
  };

  const fetchInitialNetworkData = async () => {
    try {
      const data = await fetchInitialNetwork(150);
      setInitialNetwork(data);
      setNetworkData(data);
      setNetworkSource(data.source || '');
    } catch (err) {
      console.error('Error loading initial network:', err);
      setError(err.message);
    }
  };

  const loadHubMetrics = async (gene, fromNetwork = null) => {
    // Prefer metrics embedded on the hub node from the rewiring payload
    const hubNode = fromNetwork?.nodes?.find(
      (n) => n.type === 'hub' || n.id === gene
    );
    if (hubNode?.embeddingDrift != null || hubNode?.egatProb != null) {
      setHubMetrics({
        gene,
        rank: hubNode.rank,
        egatProb: hubNode.egatProb,
        embeddingDrift: hubNode.embeddingDrift,
        sksDisruption: hubNode.sksDisruption,
      });
      return;
    }
    try {
      const info = await fetchEagatGene(gene);
      setHubMetrics(info);
    } catch {
      setHubMetrics(null);
    }
  };

  const handleAnalyze = async () => {
    if (!geneInput.trim()) return;

    setLoading(true);
    setError(null);

    try {
      const data = await fetchRewiring(geneInput.trim());
      setNetworkData(data);
      setSelectedGene(geneInput.trim().toUpperCase());
      await loadHubMetrics(geneInput.trim().toUpperCase(), data);
    } catch (err) {
      setError(err.message);
      setNetworkData(initialNetwork);
    } finally {
      setLoading(false);
    }
  };

  const handleNodeClick = (geneId) => {
    if (geneId) {
      setGeneInput(geneId);
      setSelectedGene(geneId);
      fetchRewiringForGene(geneId);
    }
  };

  const fetchRewiringForGene = async (gene) => {
    setLoading(true);
    setError(null);

    try {
      const data = await fetchRewiring(gene);
      setNetworkData(data);
      setSelectedGene(gene);
      await loadHubMetrics(gene, data);
    } catch (err) {
      setError(err.message);
      setNetworkData(initialNetwork);
    } finally {
      setLoading(false);
    }
  };

  const calculateStats = () => {
    // Prefer server-side stats (weight-delta τ classification) when present
    if (networkData?.stats) {
      return {
        total: networkData.stats.total ?? 0,
        gain: networkData.stats.gained ?? 0,
        loss: networkData.stats.lost ?? 0,
        maintained: networkData.stats.maintained ?? 0,
        tau: networkData.tau ?? 0.3,
      };
    }
    if (!networkData?.links) return { total: 0, gain: 0, loss: 0, maintained: 0, tau: 0.3 };

    return networkData.links.reduce(
      (acc, link) => {
        acc.total++;
        if (link.status === 'gained') acc.gain++;
        else if (link.status === 'lost') acc.loss++;
        else if (link.status === 'maintained') acc.maintained++;
        return acc;
      },
      { total: 0, gain: 0, loss: 0, maintained: 0, tau: networkData.tau ?? 0.3 }
    );
  };

  const stats = calculateStats();
  const mode =
    selectedGene && networkData?.nodes?.some((n) => n.type === 'hub')
      ? 'ego'
      : 'full';

  return (
    <div className="h-full flex flex-col p-6 gap-6">
      {/* Controls Panel */}
      <div className="flex-none">
        <div className="glass-card p-6">
          <div className="flex flex-wrap items-center gap-2 mb-4">
            <span className="text-xs uppercase tracking-wider text-[#A6AEB8]">
              Edge style driven by
            </span>
            <span className="px-2 py-0.5 rounded text-[10px] font-medium uppercase tracking-wider bg-[#00F0FF]/10 border border-[#00F0FF]/30 text-[#00F0FF]">
              EA-GAT Attention (α_ij)
            </span>
            {networkSource && (
              <span className="px-2 py-0.5 rounded text-[10px] font-medium uppercase tracking-wider bg-white/5 border border-white/10 text-[#A6AEB8]">
                skeleton: {networkSource}
              </span>
            )}
          </div>

          <div className="flex items-end gap-4">
            <div className="flex-1">
              <label className="block text-sm text-[#A6AEB8] mb-2">Gene Symbol</label>
              <input
                type="text"
                value={geneInput}
                onChange={(e) => setGeneInput(e.target.value.toUpperCase())}
                onKeyPress={(e) => e.key === 'Enter' && handleAnalyze()}
                placeholder="e.g., PARK2, TP53, EGFR"
                className="w-full px-4 py-3 rounded-xl bg-white/[0.03] border border-white/[0.08]
                           text-[#F2F4F8] placeholder-[#A6AEB8]/50
                           focus:outline-none focus:border-[#00F0FF]/50 focus:ring-1 focus:ring-[#00F0FF]/30
                           transition-all duration-200"
              />
            </div>

            <button
              onClick={handleAnalyze}
              disabled={loading || !geneInput.trim()}
              className="px-6 py-3 rounded-xl bg-[#00F0FF]/10 border border-[#00F0FF]/40
                         text-[#00F0FF] font-medium hover:bg-[#00F0FF]/20 hover:border-[#00F0FF]/60
                         disabled:opacity-50 disabled:cursor-not-allowed
                         transition-all duration-200 flex items-center gap-2"
            >
              {loading ? (
                <>
                  <svg className="w-5 h-5 animate-spin" fill="none" viewBox="0 0 24 24">
                    <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
                    <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z" />
                  </svg>
                  Analyzing...
                </>
              ) : (
                <>
                  <svg className="w-5 h-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" />
                  </svg>
                  Analyze Hub
                </>
              )}
            </button>

            <div className="toggle-pill">
              <div
                className="toggle-slider"
                style={{
                  left: tissue === 'normal' ? '4px' : 'calc(50% + 2px)',
                  background: tissue === 'normal'
                    ? 'linear-gradient(135deg, rgba(0, 240, 255, 0.3), rgba(0, 240, 255, 0.1))'
                    : 'linear-gradient(135deg, rgba(255, 45, 141, 0.3), rgba(255, 45, 141, 0.1))',
                  border: `1px solid ${tissue === 'normal' ? 'rgba(0, 240, 255, 0.4)' : 'rgba(255, 45, 141, 0.4)'}`,
                }}
              />
              <button
                onClick={() => setTissue('normal')}
                className={`toggle-option py-2 ${tissue === 'normal' ? 'text-[#F2F4F8]' : 'text-[#A6AEB8]'}`}
              >
                <span className="flex items-center justify-center gap-2">
                  <span
                    className="w-2 h-2 rounded-full"
                    style={{ backgroundColor: '#00F0FF', boxShadow: '0 0 8px #00F0FF' }}
                  />
                  Normal
                </span>
              </button>
              <button
                onClick={() => setTissue('tumor')}
                className={`toggle-option py-2 ${tissue === 'tumor' ? 'text-[#F2F4F8]' : 'text-[#A6AEB8]'}`}
              >
                <span className="flex items-center justify-center gap-2">
                  <span
                    className="w-2 h-2 rounded-full"
                    style={{ backgroundColor: '#FF2D8D', boxShadow: '0 0 8px #FF2D8D' }}
                  />
                  Tumor
                </span>
              </button>
            </div>
          </div>

          {/* Quick-pick EA-GAT top genes */}
          {topGenes.length > 0 && (
            <div className="mt-4 flex flex-wrap items-center gap-2">
              <span className="text-xs text-[#A6AEB8] mr-1">Top Δh:</span>
              {topGenes.slice(0, 10).map((g) => (
                <button
                  key={g.gene}
                  onClick={() => {
                    setGeneInput(g.gene);
                    setSelectedGene(g.gene);
                    fetchRewiringForGene(g.gene);
                  }}
                  className={`px-2.5 py-1 rounded-lg text-xs font-mono-data border transition-all duration-200
                    ${selectedGene === g.gene
                      ? 'bg-[#00F0FF]/15 border-[#00F0FF]/40 text-[#00F0FF]'
                      : 'bg-white/[0.03] border-white/[0.08] text-[#A6AEB8] hover:text-[#F2F4F8] hover:border-white/20'
                    }`}
                  title={`Δh=${(g.embeddingDrift ?? 0).toFixed(3)}  p=${(g.egatProb ?? 0).toFixed(2)}`}
                >
                  #{g.rank} {g.gene}
                </button>
              ))}
            </div>
          )}

          {error && (
            <div className="mt-4 p-3 rounded-lg bg-[#FF2D8D]/10 border border-[#FF2D8D]/30 text-[#FF2D8D] text-sm">
              {error}
            </div>
          )}
        </div>
      </div>

      {/* Main Content */}
      <div className="flex-1 flex gap-6 min-h-0">
        <div className="flex-1 glass-card network-container">
          <NetworkGraph
            tissue={tissue}
            selectedGene={selectedGene}
            onNodeClick={handleNodeClick}
            width="100%"
            height="100%"
            networkData={networkData}
            mode={mode}
          />
        </div>

        {/* Stats / EA-GAT Panel */}
        {selectedGene && mode === 'ego' && (
          <motion.div
            className="w-80 glass-card p-6 flex flex-col gap-4 overflow-y-auto"
            initial={{ opacity: 0, x: 40 }}
            animate={{ opacity: 1, x: 0 }}
            transition={{ duration: 0.3 }}
          >
            <div>
              <h3 className="font-display text-2xl font-bold text-[#FFD166] mb-1">
                {selectedGene}
              </h3>
              <p className="text-xs text-[#A6AEB8] uppercase tracking-wider">
                Hub Rewiring + EA-GAT
              </p>
            </div>

            {/* EA-GAT metrics */}
            {hubMetrics && (
              <div className="space-y-2">
                <h4 className="text-xs font-medium text-[#F2F4F8] uppercase tracking-wider">
                  Deep Learning Scores
                </h4>
                <div className="grid grid-cols-1 gap-2">
                  <div className="p-3 rounded-xl bg-white/[0.03] border border-white/[0.06]">
                    <div className="text-[10px] text-[#A6AEB8] uppercase tracking-wider mb-0.5">
                      Embedding Drift (Δh)
                    </div>
                    <div className="font-mono-data text-xl font-bold text-[#00F0FF]">
                      {hubMetrics.embeddingDrift != null
                        ? Number(hubMetrics.embeddingDrift).toFixed(4)
                        : '—'}
                    </div>
                  </div>
                  <div className="p-3 rounded-xl bg-white/[0.03] border border-white/[0.06]">
                    <div className="text-[10px] text-[#A6AEB8] uppercase tracking-wider mb-0.5">
                      EA-GAT Probability
                    </div>
                    <div className="font-mono-data text-xl font-bold text-[#00F0FF]">
                      {hubMetrics.egatProb != null
                        ? Number(hubMetrics.egatProb).toFixed(3)
                        : '—'}
                    </div>
                  </div>
                  <div className="p-3 rounded-xl bg-white/[0.03] border border-white/[0.06]">
                    <div className="text-[10px] text-[#A6AEB8] uppercase tracking-wider mb-0.5">
                      SKS Disruption
                    </div>
                    <div className="font-mono-data text-lg font-bold text-[#A6AEB8]">
                      {hubMetrics.sksDisruption != null && hubMetrics.sksDisruption > 0
                        ? Number(hubMetrics.sksDisruption).toFixed(3)
                        : 'awaiting SKS'}
                    </div>
                  </div>
                  {hubMetrics.rank != null && (
                    <div className="text-xs text-[#A6AEB8]">
                      Global rank by Δh:{' '}
                      <span className="text-[#F2F4F8] font-mono-data">#{hubMetrics.rank}</span>
                    </div>
                  )}
                </div>
              </div>
            )}

            <div className="space-y-3">
              <div className="p-4 rounded-xl bg-white/[0.03] border border-white/[0.06]">
                <div className="text-xs text-[#A6AEB8] uppercase tracking-wider mb-1">
                  Total Connections
                </div>
                <div className="font-mono-data text-3xl font-bold text-[#00F0FF]">
                  {stats.total}
                </div>
              </div>

              <div className="grid grid-cols-3 gap-2">
                <div className="p-2 rounded-lg bg-white/[0.03] border border-white/[0.06] text-center">
                  <div className="text-[10px] text-[#A6AEB8]">Gained</div>
                  <div className="font-mono-data text-lg font-bold text-[#FF2D8D]">{stats.gain}</div>
                  <div className="text-[9px] text-[#A6AEB8]/70 mt-0.5">Δw &gt; τ</div>
                </div>
                <div className="p-2 rounded-lg bg-white/[0.03] border border-white/[0.06] text-center">
                  <div className="text-[10px] text-[#A6AEB8]">Lost</div>
                  <div className="font-mono-data text-lg font-bold text-[#00F0FF]">{stats.loss}</div>
                  <div className="text-[9px] text-[#A6AEB8]/70 mt-0.5">Δw &lt; −τ</div>
                </div>
                <div className="p-2 rounded-lg bg-white/[0.03] border border-white/[0.06] text-center">
                  <div className="text-[10px] text-[#A6AEB8]">Stable</div>
                  <div className="font-mono-data text-lg font-bold text-[#9B59D6]">{stats.maintained}</div>
                  <div className="text-[9px] text-[#A6AEB8]/70 mt-0.5">|Δw| ≤ τ</div>
                </div>
              </div>
              <div className="text-[10px] text-[#A6AEB8]/80 font-mono-data">
                τ = {stats.tau ?? 0.3} · Δw = w_tumor − w_normal (Pearson distance)
              </div>
            </div>

            <div className="pt-4 border-t border-white/10">
              <h4 className="text-sm font-medium text-[#F2F4F8] uppercase tracking-wider mb-2">
                Neighbors
              </h4>
              <div className="space-y-2 max-h-40 overflow-y-auto">
                {networkData?.nodes
                  ?.filter((n) => n.type !== 'hub')
                  ?.slice(0, 12)
                  .map((node) => (
                    <div
                      key={node.id}
                      className="flex items-center justify-between p-2 rounded-lg bg-white/[0.02] hover:bg-white/[0.05] cursor-pointer transition-colors"
                      onClick={() => {
                        setGeneInput(node.id);
                        fetchRewiringForGene(node.id);
                      }}
                    >
                      <div className="flex items-center gap-2">
                        <div
                          className="w-2 h-2 rounded-full"
                          style={{
                            backgroundColor:
                              node.status === 'lost'
                                ? '#00F0FF'
                                : node.status === 'gained'
                                  ? '#FF2D8D'
                                  : '#9B59D6',
                          }}
                        />
                        <span className="text-sm text-[#F2F4F8] font-mono-data">{node.id}</span>
                      </div>
                      {node.embeddingDrift != null && (
                        <span className="text-[10px] font-mono-data text-[#A6AEB8]">
                          Δh {Number(node.embeddingDrift).toFixed(2)}
                        </span>
                      )}
                    </div>
                  ))}
              </div>
            </div>
          </motion.div>
        )}

        {/* When no hub selected: show EA-GAT top list */}
        {!selectedGene && topGenes.length > 0 && (
          <motion.div
            className="w-80 glass-card p-6 flex flex-col gap-3 overflow-y-auto"
            initial={{ opacity: 0, x: 40 }}
            animate={{ opacity: 1, x: 0 }}
          >
            <div>
              <h3 className="font-display text-lg font-bold text-[#F2F4F8]">
                EA-GAT Leaderboard
              </h3>
              <p className="text-xs text-[#A6AEB8] mt-1">
                Ranked by Embedding Drift (Δh). Click to analyze hub.
              </p>
            </div>
            <div className="space-y-1.5">
              {topGenes.map((g) => (
                <button
                  key={g.gene}
                  onClick={() => {
                    setGeneInput(g.gene);
                    setSelectedGene(g.gene);
                    fetchRewiringForGene(g.gene);
                  }}
                  className="w-full flex items-center justify-between p-2.5 rounded-lg
                             bg-white/[0.02] border border-white/[0.06]
                             hover:border-[#00F0FF]/30 hover:bg-[#00F0FF]/5 transition-all text-left"
                >
                  <div className="flex items-center gap-2">
                    <span className="text-xs font-mono-data text-[#A6AEB8] w-6">
                      #{g.rank}
                    </span>
                    <span className="font-mono-data text-sm font-bold text-[#F2F4F8]">
                      {g.gene}
                    </span>
                  </div>
                  <span className="font-mono-data text-xs text-[#00F0FF]">
                    {(g.embeddingDrift ?? 0).toFixed(3)}
                  </span>
                </button>
              ))}
            </div>
          </motion.div>
        )}
      </div>
    </div>
  );
};

export default HubRewiring;
