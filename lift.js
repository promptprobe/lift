/**
 * lift — A/B test sample size and peeking-risk planner.
 *
 * Two-sample binomial (unpooled) equal-allocation sample size:
 *
 *   n_per_variant = ((z_alpha + z_beta)^2 * (p1*(1-p1) + p2*(1-p2))) / delta^2
 *
 *   delta = |p2 - p1|
 *
 * z_alpha is the standard-normal quantile for the *tail* used by the test:
 *   two-sided:  z_{alpha/2}     (95% → alpha 0.05 → z_0.025 ≈ 1.960)
 *   one-sided:  z_alpha         (95% → alpha 0.05 → z_0.05  ≈ 1.645)
 *
 * z_beta = z_{1-power} wait, no: z_beta = Φ^{-1}(power) = Φ^{-1}(1-beta).
 *
 * Multi-variant (k = 2..6): Bonferroni on the alpha used for z_alpha only.
 *   comparisons = k - 1  (control vs each treatment)
 *   alpha_adj   = alpha / (k - 1)
 *
 * This file is DOM-free. Node: require('./lift.js'). Browser: window.Lift.
 */

'use strict';

// ---------------------------------------------------------------------------
// Exact standard-normal quantiles (spec table). Used as-is when the
// requested tail/power hits these points so tests and UI match textbooks.
// ---------------------------------------------------------------------------
var Z_90_TWO_SIDED = 1.6448536269514722; // alpha 0.10 two-sided = z_0.05
var Z_95_TWO_SIDED = 1.959963984540054;  // alpha 0.05 two-sided = z_0.025
var Z_99_TWO_SIDED = 2.5758293035489004; // alpha 0.01 two-sided = z_0.005
var Z_POWER_80 = 0.8416212335729143;     // beta 0.20
var Z_POWER_90 = 1.2815515655446004;     // beta 0.10

/**
 * Pocock-like n inflation for ~5 interim looks, two-sided alpha 0.05.
 *
 * Pocock (1977) constant-boundary designs use a larger critical z than
 * the fixed-sample z_0.025 ≈ 1.96. A classroom shortcut is 1.2× on z,
 * and because n ∝ z^2 that is 1.2^2 = 1.44× on sample size.
 *
 * ~2.0× (a crude always-valid / mixture bound) is too harsh here.
 * This is a fudge factor, NOT a full group-sequential or mSPRT design.
 */
var POCOCK_N_INFLATION = 1.44;

var DEFAULTS = {
  p1Pct: 4,
  mdePct: 15,
  mdeType: 'relative',
  direction: 'up',
  dailyVisitors: 800,
  variants: 2,
  significance: 95,
  power: 80,
  sided: 'two',
  peeking: false
};

var PRESETS = [
  {
    id: 'saas',
    label: 'SaaS',
    p1Pct: 4,
    mdePct: 15,
    mdeType: 'relative',
    dailyVisitors: 800,
    variants: 2
  },
  {
    id: 'ecom',
    label: '쇼핑',
    p1Pct: 2.2,
    mdePct: 12,
    mdeType: 'relative',
    dailyVisitors: 2500,
    variants: 2
  },
  {
    id: 'lead',
    label: '콘텐츠',
    p1Pct: 12,
    mdePct: 10,
    mdeType: 'relative',
    dailyVisitors: 400,
    variants: 2
  }
];

// ---------------------------------------------------------------------------
// Inverse normal (probit)
// ---------------------------------------------------------------------------

/**
 * Φ^{-1}(p) — Acklam's rational approximation (relative error < 1.15e-9),
 * with exact overrides for the spec table so 95%/80% hits 8155 on the
 * sanity check rather than a 1-off from a slightly different z.
 *
 * @param {number} p  in (0, 1)
 * @returns {number}
 */
function invNorm(p) {
  if (!(p > 0 && p < 1)) {
    throw new Error('invNorm: p must be in (0, 1)');
  }

  // Exact table hits (and their lower tails).
  var exact = exactProbit(p);
  if (exact != null) return exact;

  // Acklam 2002 coefficients.
  var a = [
    -3.969683028665376e+01,
    2.209460984245205e+02,
    -2.759285104469687e+02,
    1.383577509590705e+02,
    -3.066479806614716e+01,
    2.506628277459239e+00
  ];
  var b = [
    -5.447609879822406e+01,
    1.615858368580409e+02,
    -1.556989798598866e+02,
    6.680131188771972e+01,
    -1.328068155288572e+01
  ];
  var c = [
    -7.784894002430293e-03,
    -3.223964580411365e-01,
    -2.400758277161838e+00,
    -2.549732539343734e+00,
    4.374664141464858e+00,
    2.938163982698783e+00
  ];
  var d = [
    7.784695709041462e-03,
    3.224671290700398e-01,
    2.445134137142996e+00,
    3.754408661907416e+00
  ];

  var pLow = 0.02425;
  var pHigh = 1 - pLow;
  var q, r;

  if (p < pLow) {
    q = Math.sqrt(-2 * Math.log(p));
    return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
      ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  }
  if (p <= pHigh) {
    q = p - 0.5;
    r = q * q;
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
      (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
  }
  q = Math.sqrt(-2 * Math.log(1 - p));
  return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
    ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
}

function exactProbit(p) {
  var known = [
    [0.8, Z_POWER_80],
    [0.9, Z_POWER_90],
    [0.95, Z_90_TWO_SIDED],
    [0.975, Z_95_TWO_SIDED],
    [0.995, Z_99_TWO_SIDED]
  ];
  for (var i = 0; i < known.length; i++) {
    var pk = known[i][0];
    var z = known[i][1];
    if (Math.abs(p - pk) < 1e-12) return z;
    if (Math.abs(p - (1 - pk)) < 1e-12) return -z;
  }
  return null;
}

/**
 * z used for the alpha (Type I) term.
 * two-sided: z_{alpha/2}; one-sided: z_alpha.
 * Same table: 0.05 one-sided → 1.64485 (not 1.96).
 *
 * @param {number} alpha  in (0, 1), e.g. 0.05
 * @param {string} [sided='two']  'one' | 'two' | 'one-sided' | 'two-sided'
 * @returns {number}
 */
function zFromAlpha(alpha, sided) {
  var a = Number(alpha);
  if (!(a > 0 && a < 1)) {
    throw new Error('zFromAlpha: alpha must be in (0, 1)');
  }
  var side = normalizeSided(sided);
  var tail = side === 'one' ? a : a / 2;
  return invNorm(1 - tail);
}

function zFromPower(power) {
  var p = asUnit(power);
  if (!(p > 0 && p < 1)) {
    throw new Error('zFromPower: power must be in (0, 1)');
  }
  return invNorm(p);
}

function normalizeSided(sided) {
  var s = String(sided == null ? 'two' : sided).toLowerCase();
  if (s === 'one' || s === 'one-sided' || s === '1') return 'one';
  return 'two';
}

/** 80 → 0.80, 0.8 → 0.8. Values in (1, 100] treated as percents. */
function asUnit(x) {
  var n = Number(x);
  if (!Number.isFinite(n)) return NaN;
  return n > 1 ? n / 100 : n;
}

function significanceToAlpha(sig) {
  var s = Number(sig);
  if (!Number.isFinite(s)) return NaN;
  if (s > 1) return (100 - s) / 100; // 95 → 0.05
  return 1 - s;                       // 0.95 → 0.05
}

function fail(message) {
  return { ok: false, error: message };
}

function isProportion(x) {
  return Number.isFinite(x) && x > 0 && x < 1;
}

// ---------------------------------------------------------------------------
// Rates
// ---------------------------------------------------------------------------

/**
 * Resolve control / treatment rates from either raw p1/p2 or MDE inputs.
 *
 * Relative (default): p2 = p1 * (1 ± mde). mde=0.20 on p1=0.05 → 0.06.
 * Absolute:          p2 = p1 ± mde.       mde=0.01 on p1=0.05 → 0.06.
 * Default direction is improvement (p2 > p1).
 *
 * @returns {{ok:true, p1:number, p2:number}|{ok:false, error:string}}
 */
function resolveRates(input) {
  input = input || {};
  var p1 = input.p1;
  if (p1 == null && input.p1Pct != null) p1 = Number(input.p1Pct) / 100;
  p1 = Number(p1);

  if (!Number.isFinite(p1)) {
    return fail('기존 전환율을 넣으세요.');
  }
  if (p1 <= 0 || p1 >= 1) {
    return fail('기존 전환율은 0%와 100% 사이.');
  }

  var p2 = input.p2;
  if (p2 == null && input.p2Pct != null) p2 = Number(input.p2Pct) / 100;

  if (p2 == null) {
    var type = String(input.mdeType || 'relative').toLowerCase();
    if (type === 'abs') type = 'absolute';
    if (type === 'rel') type = 'relative';
    var dir = String(input.direction || 'up').toLowerCase();
    var sign = (dir === 'down' || dir === 'drop' || dir === '-' || dir === 'decrease') ? -1 : 1;

    var mde = input.mde;
    if (mde == null && input.mdePct != null) {
      mde = Number(input.mdePct) / 100;
    }
    mde = Number(mde);
    if (!Number.isFinite(mde) || mde <= 0) {
      return fail('최소 감지 효과는 0보다 큼.');
    }

    if (type === 'absolute') {
      p2 = p1 + sign * mde;
    } else {
      p2 = p1 * (1 + sign * mde);
    }
  }

  p2 = Number(p2);
  if (!Number.isFinite(p2)) {
    return fail('처리 전환율을 계산하지 못함.');
  }
  if (p2 <= 0 || p2 >= 1) {
    return fail('처리 전환율이 0%와 100% 사이여야 함. 효과를 줄이거나 방향을 바꾸세요.');
  }
  if (Math.abs(p2 - p1) === 0) {
    return fail('효과가 0이라 두 안이 같음.');
  }
  return { ok: true, p1: p1, p2: p2 };
}

// ---------------------------------------------------------------------------
// Sample size
// ---------------------------------------------------------------------------

/**
 * Unpooled two-sample binomial n per variant (ceiled).
 *
 * @param {object} opts
 * @param {number} opts.p1
 * @param {number} opts.p2
 * @param {number} [opts.alpha=0.05]
 * @param {number} [opts.significance]  90 | 95 | 99  (alternative to alpha)
 * @param {number} [opts.power=0.8]     0.8 or 80
 * @param {string} [opts.sided='two']
 * @param {number} [opts.variants=2]
 * @returns {object}
 */
function sampleSize(opts) {
  opts = opts || {};
  var p1 = Number(opts.p1);
  var p2 = Number(opts.p2);

  if (!isProportion(p1)) {
    return fail('기존 전환율은 0%와 100% 사이.');
  }
  if (!isProportion(p2)) {
    return fail('처리 전환율은 0%와 100% 사이.');
  }

  var delta = Math.abs(p2 - p1);
  if (!(delta > 0) || !Number.isFinite(delta)) {
    return fail('효과가 0이라 두 안이 같음.');
  }

  var variants = opts.variants == null ? 2 : Number(opts.variants);
  if (!Number.isInteger(variants) || variants < 2 || variants > 6) {
    return fail('안 개수는 2에서 6.');
  }

  var alpha;
  if (opts.alpha != null) {
    alpha = Number(opts.alpha);
    if (alpha > 1) alpha = alpha / 100;
  } else if (opts.significance != null) {
    alpha = significanceToAlpha(opts.significance);
  } else {
    alpha = 0.05;
  }
  if (!(alpha > 0 && alpha < 1)) {
    return fail('유의수준이 범위를 벗어남.');
  }

  var power = opts.power == null ? 0.8 : asUnit(opts.power);
  if (!(power > 0 && power < 1)) {
    return fail('검정력은 0%와 100% 사이.');
  }

  var sided = normalizeSided(opts.sided);
  var comparisons = variants - 1;
  var alphaAdj = alpha / comparisons;

  var zAlpha;
  var zBeta;
  try {
    zAlpha = zFromAlpha(alphaAdj, sided);
    zBeta = zFromPower(power);
  } catch (err) {
    return fail(err.message || 'z를 계산하지 못함.');
  }

  var varSum = p1 * (1 - p1) + p2 * (1 - p2);
  var nRaw = (Math.pow(zAlpha + zBeta, 2) * varSum) / (delta * delta);
  if (!Number.isFinite(nRaw) || nRaw < 0) {
    return fail('이 비율로는 n을 구하지 못함.');
  }
  var nPerVariant = Math.ceil(nRaw);

  return {
    ok: true,
    p1: p1,
    p2: p2,
    delta: delta,
    relLift: (p2 - p1) / p1,
    variants: variants,
    comparisons: comparisons,
    alpha: alpha,
    alphaAdj: alphaAdj,
    power: power,
    sided: sided,
    zAlpha: zAlpha,
    zBeta: zBeta,
    nRaw: nRaw,
    nPerVariant: nPerVariant,
    totalN: nPerVariant * variants
  };
}

// ---------------------------------------------------------------------------
// Full plan
// ---------------------------------------------------------------------------

/**
 * End-to-end planner used by the UI and tests.
 *
 * Accepts UI-style fields (p1Pct, mdePct, …) or math fields (p1, p2, alpha).
 */
function planTest(input) {
  input = input || {};
  var rates = resolveRates(input);
  if (!rates.ok) return rates;

  var daily = input.dailyVisitors == null ? 0 : Number(input.dailyVisitors);
  if (!Number.isFinite(daily) || daily < 0) {
    return fail('하루 방문자는 0 이상.');
  }

  var sized = sampleSize({
    p1: rates.p1,
    p2: rates.p2,
    alpha: input.alpha,
    significance: input.significance,
    power: input.power,
    sided: input.sided,
    variants: input.variants
  });
  if (!sized.ok) return sized;

  var days = null;
  if (daily > 0) {
    days = Math.ceil(sized.totalN / daily);
  }

  var peekingOn = Boolean(input.peeking);
  // Naive independent-look sketch: L = max(days, 1). If traffic is 0 we
  // cannot count calendar days, so L = 1 (a single end-of-test look).
  var looks = days == null ? 1 : Math.max(days, 1);
  var naiveFP = 1 - Math.pow(1 - sized.alpha, looks);
  var nPerConservative = Math.ceil(sized.nPerVariant * POCOCK_N_INFLATION);
  var totalConservative = nPerConservative * sized.variants;
  var daysConservative = daily > 0 ? Math.ceil(totalConservative / daily) : null;

  return {
    ok: true,
    p1: sized.p1,
    p2: sized.p2,
    delta: sized.delta,
    relLift: sized.relLift,
    nPerVariant: sized.nPerVariant,
    totalN: sized.totalN,
    nRaw: sized.nRaw,
    days: days,
    dailyVisitors: daily,
    expectedControl: sized.nPerVariant * sized.p1,
    expectedTreatment: sized.nPerVariant * sized.p2,
    variants: sized.variants,
    comparisons: sized.comparisons,
    alpha: sized.alpha,
    alphaAdj: sized.alphaAdj,
    power: sized.power,
    sided: sized.sided,
    zAlpha: sized.zAlpha,
    zBeta: sized.zBeta,
    peeking: {
      enabled: peekingOn,
      looks: looks,
      naiveFP: naiveFP,
      inflation: POCOCK_N_INFLATION,
      nPerConservative: nPerConservative,
      totalConservative: totalConservative,
      daysConservative: daysConservative
    }
  };
}

/**
 * n-per-variant vs relative MDE at the current baseline (for the sparkline).
 * Points whose treatment rate would leave (0, 1) are skipped.
 */
function mdeCurve(input, fromRel, toRel, steps) {
  fromRel = fromRel == null ? 0.05 : fromRel;
  toRel = toRel == null ? 0.4 : toRel;
  steps = steps == null ? 15 : steps;
  var points = [];
  var base = Object.assign({}, input);
  for (var i = 0; i <= steps; i++) {
    var mde = fromRel + (toRel - fromRel) * (i / steps);
    var trial = Object.assign({}, base, {
      mde: mde,
      mdePct: mde * 100,
      mdeType: 'relative',
      p2: undefined,
      p2Pct: undefined
    });
    var plan = planTest(trial);
    if (plan.ok) {
      points.push({ mde: mde, n: plan.nPerVariant, p2: plan.p2 });
    }
  }
  return points;
}

function formatPlan(plan, extras) {
  extras = extras || {};
  if (!plan || !plan.ok) {
    return plan && plan.error ? plan.error : '결과를 만들지 못함.';
  }
  var sigPct = Math.round((1 - plan.alpha) * 100);
  var pwrPct = Math.round(plan.power * 100);
  var side = plan.sided === 'one' ? '단측' : '양측';
  var relPct = (plan.relLift * 100);
  var dirWord = plan.relLift >= 0 ? '상승' : '하락';
  var daysBit = plan.days == null ? '' : ' · ' + plan.days + '일';
  var peekBit = (plan.peeking && plan.peeking.enabled)
    ? ' 중간에 보면 거짓양이 늘어납니다.'
    : '';
  return '안 ' + plan.variants +
    ' · 기존 ' + pct(plan.p1, 2) +
    ' · 상대 ' + dirWord + ' ' +
    abs(relPct).toFixed(Math.abs(relPct) >= 10 ? 0 : 1) + '%' +
    ' (' + pct(plan.p2, 2) + ')' +
    ' · ' + side + ' ' + sigPct + '%' +
    ' · 검정력 ' + pwrPct + '%' +
    ' · 안당 ' + plan.nPerVariant.toLocaleString('ko-KR') + '명' +
    ' · 전체 ' + plan.totalN.toLocaleString('ko-KR') + '명' +
    daysBit +
    ' · 기존 전환 ' + Math.round(plan.expectedControl).toLocaleString('ko-KR') +
    peekBit;
}

function sigLabel(alpha) {
  return (alpha * 100).toFixed(alpha < 0.01 ? 2 : 0).replace(/\.00$/, '') + '%';
}

function pct(p, digits) {
  digits = digits == null ? 2 : digits;
  return (p * 100).toFixed(digits) + '%';
}

function abs(x) { return Math.abs(x); }

function formatAlpha(a) {
  if (a >= 0.01) return a.toFixed(3).replace(/0+$/, '').replace(/\.$/, '');
  return a.toPrecision(2);
}

// ---------------------------------------------------------------------------
// Shareable hash
// ---------------------------------------------------------------------------

var HASH_KEYS = ['p1', 'mde', 'mt', 'dir', 'dv', 'k', 'sig', 'pwr', 'side', 'peek'];

function normalizeState(raw) {
  raw = raw || {};
  var state = {
    p1Pct: numOr(raw.p1Pct, DEFAULTS.p1Pct),
    mdePct: numOr(raw.mdePct, DEFAULTS.mdePct),
    mdeType: raw.mdeType === 'absolute' || raw.mdeType === 'abs' ? 'absolute' : 'relative',
    direction: raw.direction === 'down' ? 'down' : 'up',
    dailyVisitors: numOr(raw.dailyVisitors, DEFAULTS.dailyVisitors),
    variants: clampInt(raw.variants, 2, 6, DEFAULTS.variants),
    significance: pick(Number(raw.significance), [90, 95, 99], DEFAULTS.significance),
    power: pick(Number(raw.power), [80, 90], DEFAULTS.power),
    sided: normalizeSided(raw.sided),
    peeking: Boolean(raw.peeking)
  };
  return state;
}

function stateToPlanInput(state) {
  state = normalizeState(state);
  return {
    p1Pct: state.p1Pct,
    mdePct: state.mdePct,
    mdeType: state.mdeType,
    direction: state.direction,
    dailyVisitors: state.dailyVisitors,
    variants: state.variants,
    significance: state.significance,
    power: state.power,
    sided: state.sided,
    peeking: state.peeking
  };
}

/**
 * Compact querystring (no leading #). Stable key order.
 */
function serializeHash(state) {
  state = normalizeState(state);
  var parts = [
    'p1=' + trimNum(state.p1Pct),
    'mde=' + trimNum(state.mdePct),
    'mt=' + (state.mdeType === 'absolute' ? 'abs' : 'rel'),
    'dir=' + state.direction,
    'dv=' + trimNum(state.dailyVisitors),
    'k=' + state.variants,
    'sig=' + state.significance,
    'pwr=' + state.power,
    'side=' + state.sided,
    'peek=' + (state.peeking ? '1' : '0')
  ];
  return parts.join('&');
}

function parseHash(str) {
  if (str == null) return normalizeState(DEFAULTS);
  var s = String(str);
  if (s.charAt(0) === '#') s = s.slice(1);
  if (s.charAt(0) === '?') s = s.slice(1);
  if (!s) return normalizeState(DEFAULTS);

  var map = {};
  var bits = s.split('&');
  for (var i = 0; i < bits.length; i++) {
    if (!bits[i]) continue;
    var eq = bits[i].indexOf('=');
    var k = eq === -1 ? decode(bits[i]) : decode(bits[i].slice(0, eq));
    var v = eq === -1 ? '' : decode(bits[i].slice(eq + 1));
    map[k] = v;
  }

  return normalizeState({
    p1Pct: map.p1 != null ? Number(map.p1) : DEFAULTS.p1Pct,
    mdePct: map.mde != null ? Number(map.mde) : DEFAULTS.mdePct,
    mdeType: map.mt === 'abs' || map.mt === 'absolute' ? 'absolute' : 'relative',
    direction: map.dir === 'down' ? 'down' : 'up',
    dailyVisitors: map.dv != null ? Number(map.dv) : DEFAULTS.dailyVisitors,
    variants: map.k != null ? Number(map.k) : DEFAULTS.variants,
    significance: map.sig != null ? Number(map.sig) : DEFAULTS.significance,
    power: map.pwr != null ? Number(map.pwr) : DEFAULTS.power,
    sided: map.side,
    peeking: map.peek === '1' || map.peek === 'true'
  });
}

function decode(x) {
  try { return decodeURIComponent(x.replace(/\+/g, ' ')); }
  catch (e) { return x; }
}

function numOr(v, fallback) {
  var n = Number(v);
  return Number.isFinite(n) ? n : fallback;
}

function clampInt(v, lo, hi, fallback) {
  var n = Number(v);
  if (!Number.isFinite(n)) return fallback;
  n = Math.round(n);
  if (n < lo || n > hi) return fallback;
  return n;
}

function pick(v, allowed, fallback) {
  for (var i = 0; i < allowed.length; i++) {
    if (v === allowed[i]) return v;
  }
  return fallback;
}

function trimNum(n) {
  if (!Number.isFinite(n)) return '0';
  if (Math.abs(n - Math.round(n)) < 1e-9) return String(Math.round(n));
  return String(parseFloat(n.toPrecision(8)));
}

function naiveFalsePositive(alpha, looks) {
  var a = Number(alpha);
  var L = Math.max(1, Number(looks) || 1);
  if (!(a > 0 && a < 1)) return NaN;
  return 1 - Math.pow(1 - a, L);
}

// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------

var Lift = {
  sampleSize: sampleSize,
  planTest: planTest,
  resolveRates: resolveRates,
  zFromAlpha: zFromAlpha,
  zFromPower: zFromPower,
  invNorm: invNorm,
  mdeCurve: mdeCurve,
  formatPlan: formatPlan,
  serializeHash: serializeHash,
  parseHash: parseHash,
  normalizeState: normalizeState,
  stateToPlanInput: stateToPlanInput,
  naiveFalsePositive: naiveFalsePositive,
  POCOCK_N_INFLATION: POCOCK_N_INFLATION,
  DEFAULTS: DEFAULTS,
  PRESETS: PRESETS,
  HASH_KEYS: HASH_KEYS,
  Z_95_TWO_SIDED: Z_95_TWO_SIDED,
  Z_90_TWO_SIDED: Z_90_TWO_SIDED,
  Z_99_TWO_SIDED: Z_99_TWO_SIDED,
  Z_POWER_80: Z_POWER_80,
  Z_POWER_90: Z_POWER_90
};

if (typeof module !== 'undefined' && module.exports) {
  module.exports = Lift;
}
if (typeof window !== 'undefined') {
  window.Lift = Lift;
}
