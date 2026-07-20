import { useState, useEffect } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import NetworkGraph from '../components/NetworkGraph';
import { topRewiredGenes } from '../data/mockData';
import {
  fetchInitialNetwork,
  fetchEagatGene,
} from '../lib/api';

const NetworkExplorer = () => {
  const [tissue, setTissue] = useState('normal');
  const [selectedGene, setSelectedGene] = useState(null);
  const [isSidebarOpen, setIsSidebarOpen] = useState(false);
  const [networkData, setNetworkData] = useState(null);
  const [geneDetails, setGeneDetails] = useState(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);

  useEffect(() => {
    let cancelled = false;
    (async () => {
      try {
        setLoading(true);
        const data = await fetchInitialNetwork(150);
        if (!cancelled) {
          setNetworkData(data);
          setError(null);
        }
      } catch (err) {
        if (!cancelled) {
          setError(err.message);
          // Fallback: empty graph; sidebar can still use mock ranks
          setNetworkData(null);
        }
      } finally {
        if (!cancelled) setLoading(false);
      }
    })();
    return () => { cancelled = true; };
  }, []);

  const handleNodeClick = async (geneId) => {
    setSelectedGene(geneId);
    setIsSidebarOpen(!!geneId);
    if (!geneId) {
      setGeneDetails(null);
      return;
    }
    try {
      const info = await fetchEagatGene(geneId);
      setGeneDetails({
        gene: info.gene || geneId,
        rank: info.rank,
        egatProb: info.egatProb,
        embeddingDrift: info.embeddingDrift,
        sksDisruption: info.sksDisruption,
        description: `EA-GAT embedding drift rank #${info.rank}`,
      });
    } catch {
      const mock = topRewiredGenes.find((g) => g.gene === geneId);
      setGeneDetails(
        mock || {
          gene: geneId,
          rank: '—',
          egatProb: null,
          embeddingDrift: null,
          sksDisruption: null,
          description: 'Gene present in network (no EA-GAT score).',
        }
      );
    }
  };

  const handleReset = () => {
    setSelectedGene(null);
    setIsSidebarOpen(false);
    setGeneDetails(null);
  };

  return (
    <section id="explorer" className="relative w-full min-h-screen bg-[#0A0A0A] py-12 px-4 sm:px-8">
      <div className="absolute inset-0 gradient-mesh pointer-events-none opacity-30" />

      <div className="relative z-10 max-w-[1800px] mx-auto">
        <motion.div
          className="flex flex-col sm:flex-row items-start sm:items-center justify-between gap-4 mb-8"
          initial={{ opacity: 0, y: -20 }}
          whileInView={{ opacity: 1, y: 0 }}
          viewport={{ once: true }}
          transition={{ duration: 0.5 }}
        >
          <div>
            <h2 className="font-display text-[clamp(24px,3vw,36px)] font-bold text-[#F2F4F8]">
              Network Explorer
            </h2>
            <p className="mt-1 text-sm text-[#A6AEB8]">
              Top genes by EA-GAT Embedding Drift · edges scaled by attention α_ij
            </p>
            {error && (
              <p className="mt-1 text-xs text-[#FF2D8D]">API: {error}</p>
            )}
          </div>

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
        </motion.div>

        <div className="flex flex-col lg:flex-row gap-6">
          <motion.div
            className="flex-1 glass-card network-container"
            style={{ height: '70vh', minHeight: '500px' }}
            initial={{ opacity: 0, x: -40 }}
            whileInView={{ opacity: 1, x: 0 }}
            viewport={{ once: true }}
            transition={{ duration: 0.6, delay: 0.1 }}
          >
            {loading ? (
              <div className="w-full h-full flex items-center justify-center text-[#A6AEB8] text-sm">
                Loading EA-GAT skeleton network…
              </div>
            ) : (
              <NetworkGraph
                tissue={tissue}
                selectedGene={selectedGene}
                onNodeClick={handleNodeClick}
                width="100%"
                height="100%"
                networkData={networkData}
                mode="full"
              />
            )}
          </motion.div>

          <AnimatePresence mode="wait">
            {isSidebarOpen && geneDetails && (
              <motion.div
                className="w-full lg:w-80 glass-card p-6"
                initial={{ opacity: 0, x: 40 }}
                animate={{ opacity: 1, x: 0 }}
                exit={{ opacity: 0, x: 40 }}
                transition={{ duration: 0.3, ease: [0.25, 0.46, 0.45, 0.94] }}
              >
                <div className="flex items-center justify-between mb-6">
                  <div>
                    <h3 className="font-display text-2xl font-bold text-[#FFD166]">
                      {geneDetails.gene}
                    </h3>
                    <span className="text-xs text-[#A6AEB8] uppercase tracking-wider">
                      Rank #{geneDetails.rank ?? '—'}
                    </span>
                  </div>
                  <button
                    onClick={handleReset}
                    className="w-8 h-8 rounded-lg border border-white/10 flex items-center justify-center
                               text-[#A6AEB8] hover:text-[#F2F4F8] hover:border-white/20 transition-colors"
                  >
                    <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
                    </svg>
                  </button>
                </div>

                <div className="mb-6">
                  <div className="flex items-center justify-between mb-2">
                    <span className="text-sm text-[#A6AEB8]">EA-GAT Probability</span>
                    <span className="font-mono-data text-lg font-bold text-[#00F0FF]">
                      {geneDetails.egatProb != null
                        ? Number(geneDetails.egatProb).toFixed(2)
                        : '—'}
                    </span>
                  </div>
                  {geneDetails.egatProb != null && (
                    <div className="progress-bar">
                      <div
                        className="progress-bar bg-gradient-to-r from-[#00F0FF] to-[#00F0FF]/50"
                        style={{ width: `${Number(geneDetails.egatProb) * 100}%` }}
                      />
                    </div>
                  )}
                </div>

                <div className="space-y-4 mb-6">
                  <h4 className="text-sm font-medium text-[#F2F4F8] uppercase tracking-wider">
                    Deep Learning Scores
                  </h4>

                  <div className="p-4 rounded-xl bg-white/[0.03] border border-white/[0.06]">
                    <div className="flex items-center justify-between mb-2">
                      <span className="text-xs text-[#A6AEB8]">Embedding Drift (Δh)</span>
                      <span className="font-mono-data text-sm font-bold text-[#00F0FF]">
                        {geneDetails.embeddingDrift != null
                          ? Number(geneDetails.embeddingDrift).toFixed(4)
                          : '—'}
                      </span>
                    </div>
                  </div>

                  <div className="p-4 rounded-xl bg-white/[0.03] border border-white/[0.06]">
                    <div className="flex items-center justify-between mb-2">
                      <span className="text-xs text-[#A6AEB8]">SKS Disruption</span>
                      <span className="font-mono-data text-sm font-bold text-[#A6AEB8]">
                        {geneDetails.sksDisruption != null && geneDetails.sksDisruption > 0
                          ? Number(geneDetails.sksDisruption).toFixed(3)
                          : 'awaiting SKS'}
                      </span>
                    </div>
                  </div>
                </div>

                <div className="mb-6">
                  <h4 className="text-sm font-medium text-[#F2F4F8] uppercase tracking-wider mb-2">
                    Description
                  </h4>
                  <p className="text-sm text-[#A6AEB8] leading-relaxed">
                    {geneDetails.description}
                  </p>
                </div>
              </motion.div>
            )}
          </AnimatePresence>

          {!isSidebarOpen && (
            <motion.div
              className="hidden lg:flex w-80 glass-card p-6 items-center justify-center"
              initial={{ opacity: 0 }}
              animate={{ opacity: 1 }}
              transition={{ duration: 0.3 }}
            >
              <div className="text-center">
                <div className="w-16 h-16 mx-auto mb-4 rounded-full bg-white/[0.03] border border-white/[0.08]
                                flex items-center justify-center">
                  <svg className="w-8 h-8 text-[#A6AEB8]/50" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                          d="M13 10V3L4 14h7v7l9-11h-7z" />
                  </svg>
                </div>
                <p className="text-sm text-[#A6AEB8]">
                  Click on a node to view EA-GAT Embedding Drift metrics
                </p>
              </div>
            </motion.div>
          )}
        </div>
      </div>
    </section>
  );
};

export default NetworkExplorer;
