import React, { useEffect, useMemo, useRef, useState } from "react";

/**
 * Render a v2.0 PredictionPayload whose `task.kind` is
 * `binary_per_position` or `scalar_per_position` — one score per genome
 * position drawn as a horizontal track, with protein-region annotations
 * stacked above.
 *
 * The component is data-driven from the payload alone: it does not assume
 * SARS-CoV-2, length 29903, or A/T/G/C labels. The frontend just renders
 * whatever the backend declared.
 */

const TRACK_HEIGHT = 110;
const ANNOTATION_HEIGHT = 24;
const MARGIN_LEFT = 60;
const MARGIN_RIGHT = 24;
const MARGIN_TOP = 36;
const MARGIN_BOTTOM = 40;

const PROTEIN_COLORS = [
  "#fecaca", "#fde68a", "#bbf7d0", "#bae6fd", "#c7d2fe",
  "#fbcfe8", "#fed7aa", "#a7f3d0", "#e9d5ff", "#fef08a",
];

function colorForProtein(name) {
  let h = 0;
  for (let i = 0; i < name.length; i += 1) {
    h = (h * 31 + name.charCodeAt(i)) | 0;
  }
  return PROTEIN_COLORS[Math.abs(h) % PROTEIN_COLORS.length];
}

function scoreColor(t) {
  // light gray (240,240,240) -> deep red (185, 28, 28)
  const clamped = Math.max(0, Math.min(1, t));
  const r = Math.round(240 + (185 - 240) * clamped);
  const g = Math.round(240 + (28 - 240) * clamped);
  const b = Math.round(240 + (28 - 240) * clamped);
  return `rgb(${r}, ${g}, ${b})`;
}

/**
 * Compress values into `targetBins` buckets by max-pooling. Max-pool (vs.
 * mean) preserves rare-but-strong signals at low zoom — important for binary
 * tracks where only a handful of positions carry signal.
 */
function maxPoolBins(values, targetBins) {
  if (values.length <= targetBins) {
    return values.map((value, i) => ({ start: i, end: i + 1, value }));
  }
  const binSize = values.length / targetBins;
  const bins = new Array(targetBins);
  for (let i = 0; i < targetBins; i += 1) {
    const lo = Math.floor(i * binSize);
    const hi = Math.min(values.length, Math.floor((i + 1) * binSize));
    let m = -Infinity;
    for (let j = lo; j < hi; j += 1) {
      if (values[j] > m) m = values[j];
    }
    bins[i] = { start: lo, end: hi, value: m === -Infinity ? 0 : m };
  }
  return bins;
}

function intersectInWindow(absStart, absEnd, windowStart, windowEnd) {
  const lo = Math.max(absStart, windowStart);
  const hi = Math.min(absEnd, windowEnd);
  return hi > lo ? [lo, hi] : null;
}

export default function PerPositionTrack({ payload }) {
  const containerRef = useRef(null);
  const [width, setWidth] = useState(1200);
  const [hover, setHover] = useState(null);

  useEffect(() => {
    if (!containerRef.current) return undefined;
    const observer = new ResizeObserver((entries) => {
      for (const entry of entries) {
        setWidth(Math.max(480, Math.floor(entry.contentRect.width)));
      }
    });
    observer.observe(containerRef.current);
    return () => observer.disconnect();
  }, []);

  const view = useMemo(() => {
    const totalLength = payload?.domain?.total_length ?? 0;
    const region = payload?.domain?.region ?? null;
    const windowStart = region ? region.start : 0;
    const windowEnd = region ? region.end : totalLength;
    return {
      values: payload?.predictions?.values ?? [],
      taskKind: payload?.task?.kind ?? "binary_per_position",
      valueKind: payload?.predictions?.value_kind ?? "score",
      valueRange: payload?.predictions?.value_range ?? null,
      annotations: payload?.annotations ?? {},
      totalLength,
      windowStart,
      windowEnd,
      hasRegion: Boolean(region),
    };
  }, [payload]);

  const innerWidth = Math.max(1, width - MARGIN_LEFT - MARGIN_RIGHT);
  const svgHeight = MARGIN_TOP + ANNOTATION_HEIGHT + TRACK_HEIGHT + MARGIN_BOTTOM;
  const span = Math.max(1, view.windowEnd - view.windowStart);
  const xOf = (absIdx) => MARGIN_LEFT + ((absIdx - view.windowStart) / span) * innerWidth;

  const targetBins = Math.min(view.values.length || 1, Math.max(64, Math.floor(innerWidth)));
  const bins = useMemo(() => maxPoolBins(view.values, targetBins), [view.values, targetBins]);

  let yMin = 0;
  let yMax = 1;
  if (Array.isArray(view.valueRange) && view.valueRange.length === 2) {
    [yMin, yMax] = view.valueRange;
  } else if (view.values.length > 0) {
    yMax = Math.max(...view.values, 0.0001);
  }
  const yScale = yMax - yMin || 1;

  const proteinRegions = view.annotations.protein_regions || {};
  const roi = view.annotations.region_of_interest || null;
  const proteinAnnotation = view.annotations.protein || null;

  const labelTrack = view.taskKind === "binary_per_position"
    ? "Mutation probability"
    : "Per-position score";
  const titleSuffix = view.hasRegion
    ? ` — Region [${view.windowStart}, ${view.windowEnd})`
    : ` — Whole sequence (length ${view.totalLength})`;

  return (
    <div className="w-full">
      <div className="flex items-baseline justify-between px-2 pb-2 text-xs text-gray-600">
        <span>
          <span className="font-semibold text-gray-800">{labelTrack}</span>
          {titleSuffix}
        </span>
        <span>
          value_kind: <span className="font-mono">{view.valueKind}</span>
          {proteinAnnotation ? <span className="ml-3">protein: <span className="font-mono">{proteinAnnotation}</span></span> : null}
        </span>
      </div>

      <div ref={containerRef} className="w-full">
        <svg width={width} height={svgHeight} role="img" aria-label="Per-position prediction track">
          {/* Y axis */}
          <text x={MARGIN_LEFT - 8} y={MARGIN_TOP + ANNOTATION_HEIGHT + 10} textAnchor="end" fontSize="10" fill="#555">{yMax.toFixed(2)}</text>
          <text x={MARGIN_LEFT - 8} y={MARGIN_TOP + ANNOTATION_HEIGHT + TRACK_HEIGHT} textAnchor="end" fontSize="10" fill="#555">{yMin.toFixed(2)}</text>
          <line x1={MARGIN_LEFT} x2={MARGIN_LEFT} y1={MARGIN_TOP + ANNOTATION_HEIGHT} y2={MARGIN_TOP + ANNOTATION_HEIGHT + TRACK_HEIGHT} stroke="#bbb" />

          {/* Protein region overlay strip */}
          {Object.entries(proteinRegions).map(([name, range]) => {
            // proteinRegions values are [start_1based, end_1based]. Convert to
            // 0-based half-open absolute coords, then clip to the visible window.
            if (!Array.isArray(range) || range.length < 2) return null;
            const absStart = Math.max(0, range[0] - 1);
            const absEnd = range[1];
            const clipped = intersectInWindow(absStart, absEnd, view.windowStart, view.windowEnd);
            if (!clipped) return null;
            const x1 = xOf(clipped[0]);
            const x2 = xOf(clipped[1]);
            const w = Math.max(1, x2 - x1);
            return (
              <g key={name}>
                <rect
                  x={x1}
                  y={MARGIN_TOP}
                  width={w}
                  height={ANNOTATION_HEIGHT}
                  fill={colorForProtein(name)}
                  opacity={0.65}
                  stroke="#999"
                  strokeWidth={0.5}
                />
                {w > 32 && (
                  <text x={x1 + w / 2} y={MARGIN_TOP + ANNOTATION_HEIGHT / 2 + 4} textAnchor="middle" fontSize="10" fill="#222">
                    {name}
                  </text>
                )}
              </g>
            );
          })}

          {/* Region-of-interest dashed outline */}
          {roi && (() => {
            const clipped = intersectInWindow(roi.start, roi.end, view.windowStart, view.windowEnd);
            if (!clipped) return null;
            const x1 = xOf(clipped[0]);
            const x2 = xOf(clipped[1]);
            return (
              <rect
                x={x1}
                y={MARGIN_TOP}
                width={Math.max(0, x2 - x1)}
                height={ANNOTATION_HEIGHT + TRACK_HEIGHT}
                fill="none"
                stroke="#0ea5e9"
                strokeWidth={1.5}
                strokeDasharray="4 3"
              />
            );
          })()}

          {/* Score bars */}
          {bins.length > 0 && (() => {
            const trackTop = MARGIN_TOP + ANNOTATION_HEIGHT;
            const barWidth = Math.max(0.5, innerWidth / bins.length);
            return bins.map((b, i) => {
              const t = (b.value - yMin) / yScale;
              const h = Math.max(0, Math.min(TRACK_HEIGHT, t * TRACK_HEIGHT));
              const x = MARGIN_LEFT + i * barWidth;
              return (
                <rect
                  key={i}
                  x={x}
                  y={trackTop + TRACK_HEIGHT - h}
                  width={barWidth + 0.5}
                  height={h}
                  fill={scoreColor(t)}
                  onMouseEnter={() => setHover({ bin: b, x })}
                  onMouseLeave={() => setHover(null)}
                />
              );
            });
          })()}

          {/* X axis */}
          <line
            x1={MARGIN_LEFT}
            x2={MARGIN_LEFT + innerWidth}
            y1={MARGIN_TOP + ANNOTATION_HEIGHT + TRACK_HEIGHT}
            y2={MARGIN_TOP + ANNOTATION_HEIGHT + TRACK_HEIGHT}
            stroke="#888"
          />
          <text x={MARGIN_LEFT} y={svgHeight - MARGIN_BOTTOM + 14} fontSize="10" fill="#555">
            {view.windowStart}
          </text>
          <text x={MARGIN_LEFT + innerWidth} y={svgHeight - MARGIN_BOTTOM + 14} textAnchor="end" fontSize="10" fill="#555">
            {view.windowEnd}
          </text>

          {/* Tooltip */}
          {hover && (() => {
            const absStart = hover.bin.start + view.windowStart;
            const absEnd = hover.bin.end + view.windowStart;
            const tx = Math.min(width - 170, Math.max(MARGIN_LEFT, hover.x + 6));
            return (
              <g pointerEvents="none">
                <rect x={tx} y={MARGIN_TOP - 20} width={160} height={28} fill="white" stroke="#888" rx={3} />
                <text x={tx + 6} y={MARGIN_TOP - 2} fontSize="10" fill="#222">
                  pos {absStart}{absStart === absEnd - 1 ? "" : `..${absEnd - 1}`}: max {hover.bin.value.toFixed(4)}
                </text>
              </g>
            );
          })()}
        </svg>
      </div>
    </div>
  );
}
