import { useState } from 'react';
import { motion } from 'framer-motion';
import HubRewiring from './HubRewiring';
import PathTracer from './PathTracer';
import VirtualKnockout from './VirtualKnockout';
import LeaderboardSection from './LeaderboardSection';

const MainApp = () => {
  const [activeTab, setActiveTab] = useState('hub');

  const tabs = [
    { id: 'hub', label: 'Hub Rewiring', icon: 'network' },
    { id: 'ranks', label: 'EA-GAT Ranks', icon: 'bolt' },
    { id: 'path', label: 'Path Tracer', icon: 'route' },
    { id: 'knockout', label: 'Virtual Knockout', icon: 'chart' }
  ];

  const getIcon = (iconName) => {
    const icons = {
      network: (
        <svg className="w-5 h-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                d="M19 11H5m14 0a2 2 0 012 2v6a2 2 0 01-2 2H5a2 2 0 01-2-2v-6a2 2 0 012-2m14 0V9a2 2 0 00-2-2M5 11V9a2 2 0 012-2m0 0V5a2 2 0 012-2h6a2 2 0 012 2v2M7 7h10" />
        </svg>
      ),
      route: (
        <svg className="w-5 h-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                d="M9 20l-5.447-2.724A1 1 0 013 16.382V5.618a1 1 0 011.447-.894L9 7m0 13l6-3m-6 3V7m6 10l4.553 2.276A1 1 0 0021 18.382V7.618a1 1 0 00-.553-.894L15 4m0 13V4m0 0L9 7" />
        </svg>
      ),
      chart: (
        <svg className="w-5 h-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                d="M9 19v-6a2 2 0 00-2-2H5a2 2 0 00-2 2v6a2 2 0 002 2h2a2 2 0 002-2zm0 0V9a2 2 0 012-2h2a2 2 0 012 2v10m-6 0a2 2 0 002 2h2a2 2 0 002-2m0 0V5a2 2 0 012-2h2a2 2 0 012 2v14a2 2 0 01-2 2h-2a2 2 0 01-2-2z" />
        </svg>
      ),
      bolt: (
        <svg className="w-5 h-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                d="M13 10V3L4 14h7v7l9-11h-7z" />
        </svg>
      )
    };
    return icons[iconName];
  };

  return (
    <div className="relative h-screen flex flex-col z-10">
      {/* Header */}
      <div className="flex-none px-6 py-4 glass-card-sm border-b border-white/10">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-3">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-[#00F0FF]/20 to-[#FF2D8D]/20
                            border border-[#00F0FF]/30 flex items-center justify-center">
              <svg className="w-6 h-6 text-[#00F0FF]" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                      d="M13 10V3L4 14h7v7l9-11h-7z" />
              </svg>
            </div>
            <div>
              <h1 className="font-display text-xl font-bold text-[#F2F4F8]">DiCE-Duan Analyzer</h1>
              <p className="text-xs text-[#A6AEB8]">EA-GAT · Embedding Drift · Network Rewiring</p>
            </div>
          </div>

          {/* Tab Navigation */}
          <div className="flex items-center gap-2 glass-card-sm px-2 py-2">
            {tabs.map((tab) => (
              <button
                key={tab.id}
                onClick={() => setActiveTab(tab.id)}
                className={`
                  px-4 py-2 rounded-lg text-sm font-medium transition-all duration-200
                  flex items-center gap-2
                  ${activeTab === tab.id
                    ? 'bg-[#00F0FF]/10 text-[#00F0FF] border border-[#00F0FF]/30'
                    : 'text-[#A6AEB8] hover:text-[#F2F4F8] hover:bg-white/5'
                  }
                `}
              >
                {getIcon(tab.icon)}
                {tab.label}
              </button>
            ))}
          </div>
        </div>
      </div>

      {/* Content Area */}
      <div className="flex-1 overflow-hidden">
        {activeTab === 'hub' && <HubRewiring />}
        {activeTab === 'ranks' && (
          <div className="h-full overflow-y-auto">
            <LeaderboardSection
              onGeneHighlight={(gene) => {
                if (gene) {
                  // Jump to hub rewiring with this gene pre-selected via session flag
                  sessionStorage.setItem('dice_focus_gene', gene);
                }
              }}
            />
          </div>
        )}
        {activeTab === 'path' && <PathTracer />}
        {activeTab === 'knockout' && <VirtualKnockout />}
      </div>
    </div>
  );
};

export default MainApp;
