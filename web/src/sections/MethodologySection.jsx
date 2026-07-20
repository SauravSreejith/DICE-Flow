import { useState, useEffect, useRef } from 'react';
import { motion } from 'framer-motion';
import { methodologySteps } from '../data/mockData';

const MethodologySection = () => {
  const [activeStep, setActiveStep] = useState(0);
  const intervalRef = useRef(null);

  // Auto-advance active step
  useEffect(() => {
    intervalRef.current = setInterval(() => {
      setActiveStep(prev => (prev + 1) % methodologySteps.length);
    }, 4000);

    return () => clearInterval(intervalRef.current);
  }, []);

  const getStepIcon = (iconName) => {
    const icons = {
      network: (
        <svg className="w-6 h-6" fill="none" viewBox="0 0 24 24" stroke="currentColor">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} 
                d="M19 11H5m14 0a2 2 0 012 2v6a2 2 0 01-2 2H5a2 2 0 01-2-2v-6a2 2 0 012-2m14 0V9a2 2 0 00-2-2M5 11V9a2 2 0 012-2m0 0V5a2 2 0 012-2h6a2 2 0 012 2v2M7 7h10" />
        </svg>
      ),
      calculator: (
        <svg className="w-6 h-6" fill="none" viewBox="0 0 24 24" stroke="currentColor">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} 
                d="M9 7h6m0 10v-3m-3 3h.01M9 17h.01M9 14h.01M12 14h.01M15 11h.01M12 11h.01M9 11h.01M7 21h10a2 2 0 002-2V5a2 2 0 00-2-2H7a2 2 0 00-2 2v14a2 2 0 002 2z" />
        </svg>
      ),
      'git-merge': (
        <svg className="w-6 h-6" fill="none" viewBox="0 0 24 24" stroke="currentColor">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} 
                d="M8.684 13.342C8.886 12.938 9 12.482 9 12c0-.482-.114-.938-.316-1.342m0 2.684a3 3 0 110-2.684m0 2.684l6.632 3.316m-6.632-6l6.632-3.316m0 0a3 3 0 105.367-2.684 3 3 0 00-5.367 2.684zm0 9.316a3 3 0 105.368 2.684 3 3 0 00-5.368-2.684z" />
        </svg>
      ),
      download: (
        <svg className="w-6 h-6" fill="none" viewBox="0 0 24 24" stroke="currentColor">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} 
                d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-4l-4 4m0 0l-4-4m4 4V4" />
        </svg>
      )
    };
    return icons[iconName] || icons.network;
  };

  return (
    <section id="method" className="relative w-full min-h-screen bg-[#0A0A0A] py-20 px-4 sm:px-8 lg:px-16">
      {/* Background gradient */}
      <div className="absolute inset-0 gradient-mesh pointer-events-none opacity-30" />

      <div className="relative z-10 max-w-7xl mx-auto">
        {/* Header */}
        <motion.div
          className="text-center mb-16"
          initial={{ opacity: 0, y: -20 }}
          whileInView={{ opacity: 1, y: 0 }}
          viewport={{ once: true }}
          transition={{ duration: 0.5 }}
        >
          <h2 className="font-display text-[clamp(28px,3.5vw,44px)] font-bold text-[#F2F4F8]">
            From data to candidates
          </h2>
          <p className="mt-3 text-[#A6AEB8] max-w-xl mx-auto">
            Our pipeline transforms raw expression data into actionable gene rankings
          </p>
        </motion.div>

        {/* Pipeline Steps */}
        <div className="relative">
          {/* Connecting Line */}
          <div className="hidden lg:block absolute top-[72px] left-[12%] right-[12%] h-[1px] bg-white/10">
            <motion.div
              className="h-full bg-gradient-to-r from-[#00F0FF] via-[#00F0FF] to-transparent"
              initial={{ width: '0%' }}
              whileInView={{ width: '100%' }}
              viewport={{ once: true }}
              transition={{ duration: 2, delay: 0.5, ease: [0.25, 0.46, 0.45, 0.94] }}
            />
          </div>

          {/* Steps Grid */}
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6">
            {methodologySteps.map((step, index) => (
              <motion.div
                key={step.id}
                className={`relative glass-card p-6 cursor-pointer transition-all duration-300 ${
                  activeStep === index 
                    ? 'border-[#00F0FF]/40 bg-[#00F0FF]/5' 
                    : 'border-white/[0.08] hover:border-white/20'
                }`}
                initial={{ opacity: 0, y: 60 }}
                whileInView={{ opacity: 1, y: 0 }}
                viewport={{ once: true }}
                transition={{ 
                  duration: 0.6, 
                  delay: 0.1 + index * 0.1,
                  ease: [0.25, 0.46, 0.45, 0.94]
                }}
                onClick={() => {
                  setActiveStep(index);
                  // Reset interval
                  clearInterval(intervalRef.current);
                  intervalRef.current = setInterval(() => {
                    setActiveStep(prev => (prev + 1) % methodologySteps.length);
                  }, 4000);
                }}
              >
                {/* Step Number */}
                <div 
                  className={`absolute -top-3 -left-3 w-8 h-8 rounded-full flex items-center justify-center
                              text-sm font-bold font-mono-data transition-colors duration-300 ${
                    activeStep === index
                      ? 'bg-[#00F0FF] text-[#0A0A0A]'
                      : 'bg-white/10 text-[#A6AEB8]'
                  }`}
                >
                  {step.id}
                </div>

                {/* Icon */}
                <div 
                  className={`w-12 h-12 rounded-xl flex items-center justify-center mb-4 transition-colors duration-300 ${
                    activeStep === index
                      ? 'bg-[#00F0FF]/20 text-[#00F0FF]'
                      : 'bg-white/5 text-[#A6AEB8]'
                  }`}
                >
                  {getStepIcon(step.icon)}
                </div>

                {/* Title */}
                <h3 
                  className={`font-display text-lg font-bold mb-2 transition-colors duration-300 ${
                    activeStep === index ? 'text-[#F2F4F8]' : 'text-[#A6AEB8]'
                  }`}
                >
                  {step.title}
                </h3>

                {/* Description */}
                <p className="text-sm text-[#A6AEB8]/80 leading-relaxed">
                  {step.description}
                </p>

                {/* Active Indicator */}
                {activeStep === index && (
                  <motion.div
                    className="absolute bottom-0 left-0 right-0 h-[2px] bg-[#00F0FF]"
                    layoutId="activeStepIndicator"
                    transition={{ duration: 0.3 }}
                  />
                )}
              </motion.div>
            ))}
          </div>
        </div>

        {/* Technical Details */}
        <motion.div
          className="mt-16 glass-card p-8"
          initial={{ opacity: 0, y: 40 }}
          whileInView={{ opacity: 1, y: 0 }}
          viewport={{ once: true }}
          transition={{ duration: 0.6, delay: 0.5 }}
        >
          <div className="grid grid-cols-1 md:grid-cols-3 gap-8">
            <div>
              <h4 className="text-sm font-medium text-[#F2F4F8] uppercase tracking-wider mb-4">
                Deep Learning Architecture
              </h4>
              <ul className="space-y-2">
                {['Edge-Augmented GAT (EA-GAT)', 'Anisotropic Attention Routing', 'Contrastive Embedding Drift', 'Heterophily-Aware Learning'].map((item) => (
                  <li key={item} className="flex items-center gap-2 text-sm text-[#A6AEB8]">
                    <svg className="w-4 h-4 text-[#00F0FF]" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 10V3L4 14h7v7l9-11h-7z" />
                    </svg>
                    {item}
                  </li>
                ))}
              </ul>
            </div>

            <div>
              <h4 className="text-sm font-medium text-[#F2F4F8] uppercase tracking-wider mb-4">
                Statistical Methods
              </h4>
              <ul className="space-y-2">
                {['IQR Outlier Detection', 'Z-score Normalization', 'Ensemble Aggregation', 'Bootstrap Confidence'].map((item) => (
                  <li key={item} className="flex items-center gap-2 text-sm text-[#A6AEB8]">
                    <svg className="w-4 h-4 text-[#00F0FF]" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 13l4 4L19 7" />
                    </svg>
                    {item}
                  </li>
                ))}
              </ul>
            </div>

            <div>
              <h4 className="text-sm font-medium text-[#F2F4F8] uppercase tracking-wider mb-4">
                Output Formats
              </h4>
              <ul className="space-y-2">
                {['TSV Gene Lists', 'Cytoscape Networks', 'JSON Reports', 'PDF Visualizations'].map((item) => (
                  <li key={item} className="flex items-center gap-2 text-sm text-[#A6AEB8]">
                    <svg className="w-4 h-4 text-[#00F0FF]" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 13l4 4L19 7" />
                    </svg>
                    {item}
                  </li>
                ))}
              </ul>
            </div>
          </div>
        </motion.div>

        {/* Code Snippet Preview */}
        <motion.div
          className="mt-8 glass-card overflow-hidden"
          initial={{ opacity: 0, y: 30 }}
          whileInView={{ opacity: 1, y: 0 }}
          viewport={{ once: true }}
          transition={{ duration: 0.6, delay: 0.6 }}
        >
          <div className="flex items-center gap-2 px-4 py-3 bg-white/[0.03] border-b border-white/[0.06]">
            <div className="w-3 h-3 rounded-full bg-[#FF2D8D]/60" />
            <div className="w-3 h-3 rounded-full bg-[#FFD166]/60" />
            <div className="w-3 h-3 rounded-full bg-[#00F0FF]/60" />
            <span className="ml-4 text-xs text-[#A6AEB8] font-mono-data">pipeline.py</span>
          </div>
          <div className="p-4 overflow-x-auto">
            <pre className="text-sm font-mono-data">
              <code className="text-[#A6AEB8]">
                <span className="text-[#FF2D8D]">from</span> dice_duan <span className="text-[#FF2D8D]">import</span> NetworkAnalyzer{'\n'}
                <span className="text-[#FF2D8D]">from</span> dice_duan.centrality <span className="text-[#FF2D8D]">import</span> compute_differential_centrality{'\n'}
                {'\n'}
                <span className="text-[#A6AEB8]"># Load networks</span>{'\n'}
                normal_net = NetworkAnalyzer.load(<span className="text-[#00F0FF]">'normal.tsv'</span>){'\n'}
                tumor_net = NetworkAnalyzer.load(<span className="text-[#00F0FF]">'tumor.tsv'</span>){'\n'}
                {'\n'}
                <span className="text-[#A6AEB8]"># Compute ensemble scores</span>{'\n'}
                results = compute_differential_centrality(normal_net, tumor_net){'\n'}
                results.export(<span className="text-[#00F0FF]">'rewired_genes.tsv'</span>){'\n'}
              </code>
            </pre>
          </div>
        </motion.div>
      </div>
    </section>
  );
};

export default MethodologySection;
