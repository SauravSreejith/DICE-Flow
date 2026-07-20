/**
 * DiCE-Flow API client
 * Points at the FastAPI backend (default localhost:8000).
 */

const API_BASE = import.meta.env.VITE_API_BASE || 'http://localhost:8000';

async function getJson(path) {
  const res = await fetch(`${API_BASE}${path}`);
  if (!res.ok) {
    const detail = await res.text().catch(() => res.statusText);
    throw new Error(detail || `HTTP ${res.status}`);
  }
  return res.json();
}

/** Health / EA-GAT load status */
export function fetchHealth() {
  return getJson('/api/health');
}

/**
 * EA-GAT leaderboard ranked by Embedding Drift (Δh).
 * @param {number} limit
 * @returns {Promise<{genes: Array, source: string, meta: object}>}
 */
export function fetchEagatLeaderboard(limit = 50) {
  return getJson(`/api/eagat/leaderboard?limit=${limit}`);
}

/** Per-gene EA-GAT metrics */
export function fetchEagatGene(gene) {
  return getJson(`/api/eagat/gene/${encodeURIComponent(gene)}`);
}

/** Top-gene skeleton network with attention + directed edges */
export function fetchInitialNetwork(limit = 150) {
  return getJson(`/api/network/initial?limit=${limit}`);
}

/** Ego rewiring network for a hub gene */
export function fetchRewiring(gene) {
  return getJson(`/api/rewiring?gene=${encodeURIComponent(gene)}`);
}

export { API_BASE };
