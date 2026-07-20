import { useEffect, useRef, useState, useCallback } from 'react';
import cytoscape from 'cytoscape';

const NetworkGraph = ({ 
  tissue = 'normal', 
  selectedGene, 
  onNodeClick, 
  highlightedGene,
  width = '100%',
  height = '100%',
  networkData = null,
  mode = 'full' // 'full' for initial network, 'ego' for hub rewiring, 'path' for path tracer
}) => {
  const cyRef = useRef(null);
  const containerRef = useRef(null);
  const [isReady, setIsReady] = useState(false);

  // Get status-based colors for hub rewiring mode
  const getStatusColor = (status) => {
    if (status === 'lost') return '#00F0FF'; // Blue for normal-only
    if (status === 'gained') return '#FF2D8D'; // Red for tumor-only
    if (status === 'maintained') return '#9B59D6'; // Purple for both
    return tissue === 'normal' ? '#00F0FF' : '#FF2D8D';
  };

  // Initialize Cytoscape
  useEffect(() => {
    if (!containerRef.current || cyRef.current || !networkData) return;

    const nodeColor = tissue === 'normal' ? '#00F0FF' : '#FF2D8D';
    const edgeColor = tissue === 'normal' ? 'rgba(0, 240, 255, 0.25)' : 'rgba(255, 45, 141, 0.25)';

    // Prepare elements from API data
    const elements = [
      ...networkData.nodes.map(node => {
        // Calculate size based on rank if available, otherwise use default
        let nodeSize = 16;
        if (node.rank) {
          // Formula: radius = 20 - (rank / 5)
          // Rank 1 = 20, Rank 50 = 10
          nodeSize = Math.max(10, 20 - (node.rank / 5));
        } else if (node.type === 'hub') {
          nodeSize = 24;
        }
        
        return {
          data: { 
            id: node.id, 
            label: node.label || node.id,
            type: node.type || 'normal',
            status: node.status,
            rank: node.rank,
            size: nodeSize,
            log2fc: node.log2fc,
            is_essential: node.is_essential
          }
        };
      }),
      ...networkData.links.map((edge, idx) => ({
        data: { 
          id: `e${idx}`, 
          source: edge.source, 
          target: edge.target, 
          weight: edge.weight || 0.5,
          // GAT attention coefficient (α_ij) from EA-GAT backend when available
          attention: edge.attention,
          directed: edge.directed === true || edge.directed === 1 || edge.directed === 'true',
          status: edge.status
        }
      }))
    ];

    cyRef.current = cytoscape({
      container: containerRef.current,
      elements,
      style: [
        {
          selector: 'node',
          style: {
            'background-color': (ele) => {
              if (mode === 'ego' && ele.data('status')) {
                return getStatusColor(ele.data('status'));
              }
              const l2fc = ele.data('log2fc');
              if (l2fc !== undefined && l2fc !== 0) {
                return l2fc > 0 ? '#EF4444' : '#3B82F6'; // TailWind Red-500 or Blue-500
              }
              return nodeColor;
            },
            'width': (ele) => ele.data('size') || 16,
            'height': (ele) => ele.data('size') || 16,
            'label': (ele) => ele.data('label'),
            'color': '#F2F4F8',
            'font-size': '10px',
            'font-family': 'JetBrains Mono, monospace',
            'text-valign': 'bottom',
            'text-halign': 'center',
            'text-margin-y': 6,
            'text-background-color': '#0A0A0A',
            'text-background-opacity': 0.8,
            'text-background-padding': '2px 4px',
            'text-background-shape': 'roundrectangle',
            'border-width': (ele) => {
              // Priority: is_essential -> selected -> none
              if (ele.selected()) return 3;
              if (ele.data('is_essential')) return 2;
              return 0;
            },
            'border-color': (ele) => {
              if (ele.selected()) return '#FFD166';
              if (ele.data('is_essential')) return '#10B981'; // Tailwind Emerald-500
              return '#000000';
            },
            'transition-property': 'background-color, border-width, border-color, width, height',
            'transition-duration': '0.3s'
          }
        },
        {
          selector: 'node[type="hub"]',
          style: {
            'width': 28,
            'height': 28,
            'font-size': '12px',
            'font-weight': 'bold',
            'border-width': 2,
            'border-color': '#FFD166',
            'border-opacity': 0.8
          }
        },
        {
          selector: 'node:selected',
          style: {
            'border-width': 3,
            'border-color': '#FFD166',
            'border-opacity': 1
          }
        },
        {
          selector: 'node.highlighted',
          style: {
            'border-width': 3,
            'border-color': '#FFD166',
            'border-opacity': 1,
            'background-color': '#FFD166'
          }
        },
        {
          selector: 'node.neighbor',
          style: {
            'border-width': 2,
            'border-color': nodeColor,
            'border-opacity': 0.8
          }
        },
        {
          selector: 'node.dimmed',
          style: {
            'opacity': 0.2
          }
        },
        {
          selector: 'edge',
          style: {
            // Prefer GAT attention score (α_ij) over raw physical weight when present
            'width': (ele) => {
              const attentionScore = ele.data('attention') || ele.data('weight') || 0.5;
              return Math.max(1, attentionScore * 5); // Thicker lines = High GAT Attention
            },
            'line-color': (ele) => {
              if (mode === 'ego' && ele.data('status')) {
                const status = ele.data('status');
                if (status === 'lost') return 'rgba(0, 240, 255, 0.5)';
                if (status === 'gained') return 'rgba(255, 45, 141, 0.5)';
                if (status === 'maintained') return 'rgba(155, 89, 214, 0.5)';
              }
              // Make high-attention edges glow brighter
              const alpha = ele.data('attention') ? 0.2 + (ele.data('attention') * 0.8) : 0.4;
              return tissue === 'normal' ? `rgba(0, 240, 255, ${alpha})` : `rgba(255, 45, 141, ${alpha})`;
            },
            'target-arrow-color': edgeColor,
            // Show arrows for directed edges (parsed from KEGG / merged network)
            'target-arrow-shape': (ele) => ele.data('directed') ? 'triangle' : 'none',
            'curve-style': 'bezier',
            'opacity': 0.6,
            'transition-property': 'line-color, width, opacity',
            'transition-duration': '0.3s'
          }
        },
        {
          selector: 'edge.highlighted',
          style: {
            'line-color': '#FFD166',
            'width': (ele) => {
              const attentionScore = ele.data('attention') || ele.data('weight') || 0.5;
              return Math.max(2, attentionScore * 5);
            },
            'opacity': 1
          }
        },
        {
          selector: 'edge.dimmed',
          style: {
            'opacity': 0.1
          }
        }
      ],
      layout: {
        name: mode === 'path' ? 'breadthfirst' : 'cose',
        padding: 30,
        nodeRepulsion: mode === 'ego' ? 8000 : 6000,
        nodeOverlap: 10,
        idealEdgeLength: mode === 'ego' ? 100 : 70,
        edgeElasticity: 200, // Increased for more bouncy movement
        nestingFactor: 5,
        gravity: mode === 'ego' ? 30 : 20, // Reduced gravity for more free movement
        numIter: mode === 'ego' ? 600 : 800,
        initialTemp: 200,
        coolingFactor: 0.95,
        minTemp: 1.0,
        fit: true,
        animate: true,
        animationDuration: 500,
        directed: mode === 'path',
        spacingFactor: mode === 'path' ? 1.5 : 1
      },
      minZoom: 0.3,
      maxZoom: 3,
      wheelSensitivity: 0.3
    });

    // Add event listeners
    cyRef.current.on('tap', 'node', (evt) => {
      const node = evt.target;
      onNodeClick?.(node.data('id'));
    });

    cyRef.current.on('mouseover', 'node', (evt) => {
      const node = evt.target;
      const currentSize = node.data('size') || 16;
      node.animate({
        style: { 
          'width': currentSize * 1.4,
          'height': currentSize * 1.4
        }
      }, { duration: 200 });
    });

    cyRef.current.on('mouseout', 'node', (evt) => {
      const node = evt.target;
      if (!node.selected()) {
        node.animate({
          style: { 
            'width': node.data('size') || 16,
            'height': node.data('size') || 16
          }
        }, { duration: 200 });
      }
    });

    setIsReady(true);

    return () => {
      if (cyRef.current) {
        cyRef.current.destroy();
        cyRef.current = null;
      }
    };
  }, [tissue, networkData, mode]); // Re-initialize when these change

  // Update colors when tissue changes
  useEffect(() => {
    if (!cyRef.current || !isReady) return;

    const nodeColor = tissue === 'normal' ? '#00F0FF' : '#FF2D8D';
    const edgeColor = tissue === 'normal' ? 'rgba(0, 240, 255, 0.25)' : 'rgba(255, 45, 141, 0.25)';

    cyRef.current.style()
      .selector('node')
      .style({ 'background-color': nodeColor })
      .selector('node.neighbor')
      .style({ 'border-color': nodeColor })
      .selector('edge')
      .style({ 
        'line-color': edgeColor,
        'target-arrow-color': edgeColor
      })
      .update();
  }, [tissue, isReady]);

  // Handle selected gene - auto-select when gene is selected
  useEffect(() => {
    if (!cyRef.current || !isReady) return;

    // Clear previous selection
    cyRef.current.nodes().removeClass(['highlighted', 'neighbor', 'dimmed']);
    cyRef.current.edges().removeClass(['highlighted', 'dimmed']);
    cyRef.current.nodes().unselect();

    if (selectedGene) {
      const selectedNode = cyRef.current.getElementById(selectedGene);
      
      if (selectedNode.length > 0) {
        // Select and highlight the node
        selectedNode.select();
        selectedNode.addClass('highlighted');

        if (mode === 'full') {
          // Get neighbors for full network mode
          const neighbors = selectedNode.neighborhood();
          neighbors.nodes().addClass('neighbor');
          neighbors.edges().addClass('highlighted');

          // Dim other elements
          cyRef.current.nodes().not(selectedNode).not(neighbors.nodes()).addClass('dimmed');
          cyRef.current.edges().not(neighbors.edges()).addClass('dimmed');

          // Center on selected node with more padding (less zoom)
          cyRef.current.animate({
            fit: {
              eles: selectedNode.union(neighbors),
              padding: 150 // Increased from 80 for less aggressive zoom
            }
          }, { duration: 500 });
        } else {
          // For ego network, just center on the selected node with generous padding
          cyRef.current.animate({
            fit: {
              eles: selectedNode,
              padding: 200 // Increased from 100 for less aggressive zoom
            }
          }, { duration: 500 });
        }
      }
    }
  }, [selectedGene, isReady, mode]);

  // Handle highlighted gene from leaderboard
  useEffect(() => {
    if (!cyRef.current || !isReady || !highlightedGene) return;

    const node = cyRef.current.getElementById(highlightedGene);
    if (node.length > 0) {
      const currentSize = node.data('size') || 16;
      node.animate({
        style: { 
          'width': currentSize * 1.6,
          'height': currentSize * 1.6
        }
      }, { duration: 200 });

      // Flash effect
      const originalColor = tissue === 'normal' ? '#00F0FF' : '#FF2D8D';
      node.animate({
        style: { 'background-color': '#FFD166' }
      }, { duration: 150 })
      .animate({
        style: { 'background-color': originalColor }
      }, { duration: 150 })
      .animate({
        style: { 'background-color': '#FFD166' }
      }, { duration: 150 })
      .animate({
        style: { 'background-color': originalColor }
      }, { duration: 150 });

      return () => {
        if (node.length > 0) {
          node.stop();
          node.animate({
            style: { 
              'width': node.data('size') || 16,
              'height': node.data('size') || 16,
              'background-color': originalColor
            }
          }, { duration: 200 });
        }
      };
    }
  }, [highlightedGene, tissue, isReady]);

  const handleReset = useCallback(() => {
    if (!cyRef.current) return;
    
    cyRef.current.nodes().removeClass(['highlighted', 'neighbor', 'dimmed']);
    cyRef.current.edges().removeClass(['highlighted', 'dimmed']);
    cyRef.current.nodes().unselect();
    
    cyRef.current.animate({
      fit: { padding: 30 }
    }, { duration: 500 });
    
    onNodeClick?.(null);
  }, [onNodeClick]);

  const handleZoomIn = useCallback(() => {
    if (!cyRef.current) return;
    cyRef.current.zoom(cyRef.current.zoom() * 1.2);
  }, []);

  const handleZoomOut = useCallback(() => {
    if (!cyRef.current) return;
    cyRef.current.zoom(cyRef.current.zoom() * 0.8);
  }, []);

  if (!networkData) {
    return (
      <div className="w-full h-full flex items-center justify-center">
        <div className="text-center">
          <svg className="w-12 h-12 mx-auto mb-4 animate-spin text-[#00F0FF]" fill="none" viewBox="0 0 24 24">
            <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
            <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z" />
          </svg>
          <p className="text-sm text-[#A6AEB8]">Loading network...</p>
        </div>
      </div>
    );
  }

  return (
    <div className="relative w-full h-full">
      <div 
        ref={containerRef} 
        className="w-full h-full"
        style={{ width, height }}
      />
      
      {/* Zoom Controls */}
      <div className="absolute bottom-4 right-4 flex flex-col gap-2">
        <button
          onClick={handleZoomIn}
          className="w-10 h-10 rounded-lg glass-card-sm flex items-center justify-center
                     text-[#A6AEB8] hover:text-[#00F0FF] hover:border-[#00F0FF]/30
                     transition-all duration-200"
        >
          <svg className="w-5 h-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M12 4v16m8-8H4" />
          </svg>
        </button>
        <button
          onClick={handleZoomOut}
          className="w-10 h-10 rounded-lg glass-card-sm flex items-center justify-center
                     text-[#A6AEB8] hover:text-[#00F0FF] hover:border-[#00F0FF]/30
                     transition-all duration-200"
        >
          <svg className="w-5 h-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M20 12H4" />
          </svg>
        </button>
        <button
          onClick={handleReset}
          className="w-10 h-10 rounded-lg glass-card-sm flex items-center justify-center
                     text-[#A6AEB8] hover:text-[#00F0FF] hover:border-[#00F0FF]/30
                     transition-all duration-200"
          title="Reset view"
        >
          <svg className="w-5 h-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
          </svg>
        </button>
      </div>
      
      {/* Legend - different for ego mode */}
      <div className="absolute bottom-4 left-4 glass-card-sm px-4 py-3">
        {mode === 'ego' ? (
          <div className="flex items-center gap-4 text-xs flex-wrap">
            <div className="flex items-center gap-2">
              <div className="w-3 h-3 rounded-full bg-[#FF2D8D]" />
              <span className="text-[#A6AEB8]">Gained (Δw &gt; τ)</span>
            </div>
            <div className="flex items-center gap-2">
              <div className="w-3 h-3 rounded-full bg-[#00F0FF]" />
              <span className="text-[#A6AEB8]">Lost (Δw &lt; −τ)</span>
            </div>
            <div className="flex items-center gap-2">
              <div className="w-3 h-3 rounded-full bg-[#9B59D6]" />
              <span className="text-[#A6AEB8]">Stable (|Δw| ≤ τ)</span>
            </div>
          </div>
        ) : (
          <div className="flex items-center gap-4 text-xs flex-wrap">
            <div className="flex items-center gap-2">
              <div className="w-3 h-3 rounded-full bg-[#EF4444]" />
              <span className="text-[#A6AEB8]">Upregulated Hub</span>
            </div>
            <div className="flex items-center gap-2">
              <div className="w-3 h-3 rounded-full bg-[#3B82F6]" />
              <span className="text-[#A6AEB8]">Downregulated Hub</span>
            </div>
            <div className="flex items-center gap-2">
              <div className="w-3 h-3 rounded-full border-2 border-[#10B981] bg-transparent" />
              <span className="text-[#A6AEB8]">Essential</span>
            </div>
            <div className="flex items-center gap-2">
              <div 
                className="w-3 h-3 rounded-full"
                style={{ backgroundColor: tissue === 'normal' ? '#00F0FF' : '#FF2D8D' }}
              />
              <span className="text-[#A6AEB8]">Network Gene</span>
            </div>
          </div>
        )}
      </div>
    </div>
  );
};

export default NetworkGraph;