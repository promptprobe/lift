'use strict';

var L = require('./lift.js');

var passed = 0;
var failed = 0;
var failures = [];

function assert(cond, msg) {
  if (cond) {
    passed += 1;
    console.log('ok   ' + msg);
  } else {
    failed += 1;
    failures.push(msg);
    console.log('FAIL ' + msg);
  }
}

function close(a, b, tol, msg) {
  var ok = Number.isFinite(a) && Number.isFinite(b) && Math.abs(a - b) <= tol;
  assert(ok, msg + ' (got ' + a + ', expected ~' + b + ', tol ' + tol + ')');
}

// 1. Sanity: p1=0.05, p2=0.06, two-sided 95%, 80% power, 2 variants
var sanity = L.sampleSize({
  p1: 0.05,
  p2: 0.06,
  alpha: 0.05,
  power: 0.8,
  sided: 'two',
  variants: 2
});
assert(sanity.ok, 'sanity: sampleSize ok');
var nSanity = sanity.nPerVariant;
var tolSanity = Math.max(50, 0.02 * 8156);
close(nSanity, 8156, tolSanity, 'sanity n per variant ~8156');

// 2. Formula expansion: ((1.959964+0.841621)^2 * (0.05*0.95 + 0.06*0.94)) / 0.01^2
var zA = 1.959963984540054;
var zB = 0.8416212335729143;
var expectedRaw = (Math.pow(zA + zB, 2) * (0.05 * 0.95 + 0.06 * 0.94)) / Math.pow(0.01, 2);
close(sanity.nRaw, expectedRaw, 1e-6, 'sanity nRaw matches unpooled formula');
assert(sanity.nPerVariant === Math.ceil(expectedRaw), 'sanity n is ceil of formula');

// 3. z table
close(L.zFromAlpha(0.05, 'two'), 1.959963984540054, 1e-12, 'z two-sided 0.05');
close(L.zFromAlpha(0.10, 'two'), 1.6448536269514722, 1e-12, 'z two-sided 0.10');
close(L.zFromAlpha(0.01, 'two'), 2.5758293035489004, 1e-12, 'z two-sided 0.01');
close(L.zFromAlpha(0.05, 'one'), 1.6448536269514722, 1e-12, 'z one-sided 0.05');
close(L.zFromPower(0.8), 0.8416212335729143, 1e-12, 'z power 80%');
close(L.zFromPower(90), 1.2815515655446004, 1e-12, 'z power 90%');

// 4. Bonferroni: 3 variants, same rates, alpha 0.05 → alpha_adj=0.025, z~2.2414
var bonf = L.sampleSize({
  p1: 0.05,
  p2: 0.06,
  alpha: 0.05,
  power: 0.8,
  sided: 'two',
  variants: 3
});
assert(bonf.ok, 'bonferroni: sampleSize ok');
close(bonf.alphaAdj, 0.025, 1e-12, 'bonferroni alpha_adj = 0.05/2');
close(bonf.zAlpha, 2.2414, 0.01, 'bonferroni z_alpha ~2.2414');
assert(bonf.nPerVariant > sanity.nPerVariant, 'bonferroni n increases vs 2-variant');

// 5. relative 20% of 5% → p2 0.06
var rel = L.resolveRates({ p1: 0.05, mde: 0.20, mdeType: 'relative', direction: 'up' });
assert(rel.ok && Math.abs(rel.p2 - 0.06) < 1e-12, 'relative 20% of 5% → p2 0.06');

// 6. absolute 1pp of 5% → p2 0.06
var absR = L.resolveRates({ p1: 0.05, mde: 0.01, mdeType: 'absolute', direction: 'up' });
assert(absR.ok && Math.abs(absR.p2 - 0.06) < 1e-12, 'absolute 1pp of 5% → p2 0.06');

var relPlan = L.planTest({ p1Pct: 5, mdePct: 20, mdeType: 'relative', significance: 95, power: 80, variants: 2 });
var absPlan = L.planTest({ p1Pct: 5, mdePct: 1, mdeType: 'absolute', significance: 95, power: 80, variants: 2 });
assert(relPlan.ok && absPlan.ok && relPlan.nPerVariant === absPlan.nPerVariant,
  'relative 20% and absolute 1pp yield the same n');

// 7. invalid p1=0 / 1 → error object, not NaN
var z0 = L.planTest({ p1: 0, p2: 0.06 });
var z1 = L.planTest({ p1: 1, p2: 0.06 });
var s0 = L.sampleSize({ p1: 0, p2: 0.06 });
var s1 = L.sampleSize({ p1: 1, p2: 0.06 });
assert(z0.ok === false && typeof z0.error === 'string' && !Number.isNaN(z0.nPerVariant),
  'p1=0 returns error object');
assert(z1.ok === false && typeof z1.error === 'string', 'p1=1 returns error object');
assert(s0.ok === false && typeof s0.error === 'string', 'sampleSize p1=0 returns error');
assert(s1.ok === false && typeof s1.error === 'string', 'sampleSize p1=1 returns error');
assert(!('nPerVariant' in z0) || !Number.isFinite(z0.nPerVariant), 'error result has no finite n');

// 8. p2 out of range / delta 0
var dropTooFar = L.resolveRates({ p1: 0.05, mde: 1.5, mdeType: 'relative', direction: 'down' });
assert(dropTooFar.ok === false, 'relative drop that drives p2 <= 0 errors');
var same = L.sampleSize({ p1: 0.05, p2: 0.05, alpha: 0.05, power: 0.8 });
assert(same.ok === false, 'delta 0 errors');
var p2bad = L.sampleSize({ p1: 0.05, p2: 1.2 });
assert(p2bad.ok === false, 'p2 > 1 errors');

// 9. variants bounds
assert(L.sampleSize({ p1: 0.05, p2: 0.06, variants: 1 }).ok === false, 'variants=1 errors');
assert(L.sampleSize({ p1: 0.05, p2: 0.06, variants: 7 }).ok === false, 'variants=7 errors');
assert(L.sampleSize({ p1: 0.05, p2: 0.06, variants: 6 }).ok === true, 'variants=6 ok');

// 10. one-sided needs fewer visitors than two-sided
var one = L.sampleSize({ p1: 0.05, p2: 0.06, alpha: 0.05, power: 0.8, sided: 'one', variants: 2 });
assert(one.ok && one.nPerVariant < sanity.nPerVariant, 'one-sided n < two-sided n');

// 11. higher power / tighter alpha increases n
var p90 = L.sampleSize({ p1: 0.05, p2: 0.06, alpha: 0.05, power: 0.9, sided: 'two', variants: 2 });
var a99 = L.sampleSize({ p1: 0.05, p2: 0.06, alpha: 0.01, power: 0.8, sided: 'two', variants: 2 });
assert(p90.ok && p90.nPerVariant > sanity.nPerVariant, '90% power n > 80% power n');
assert(a99.ok && a99.nPerVariant > sanity.nPerVariant, '99% sig n > 95% sig n');

// 12. days, expected conversions, peeking inflation
var planned = L.planTest({
  p1: 0.05,
  p2: 0.06,
  significance: 95,
  power: 80,
  variants: 2,
  dailyVisitors: 400,
  peeking: true
});
assert(planned.ok, 'planTest ok');
assert(planned.days === Math.ceil(planned.totalN / 400), 'days = ceil(totalN / daily)');
close(planned.expectedControl, planned.nPerVariant * 0.05, 1e-9, 'expected control conversions');
assert(planned.peeking.enabled === true, 'peeking flag on');
assert(planned.peeking.looks === Math.max(planned.days, 1), 'looks = max(days, 1)');
close(planned.peeking.naiveFP, 1 - Math.pow(0.95, planned.peeking.looks), 1e-9, 'naive FP = 1-(1-a)^L');
assert(planned.peeking.inflation === 1.44, 'Pocock inflation is 1.44');
assert(planned.peeking.nPerConservative === Math.ceil(planned.nPerVariant * 1.44),
  'conservative n = ceil(n * 1.44)');

// 13. hash serialize / parse roundtrip
var state = {
  p1Pct: 4,
  mdePct: 15,
  mdeType: 'relative',
  direction: 'up',
  dailyVisitors: 800,
  variants: 2,
  significance: 95,
  power: 80,
  sided: 'two',
  peeking: true
};
var hash = L.serializeHash(state);
var back = L.parseHash(hash);
assert(typeof hash === 'string' && hash.indexOf('p1=4') !== -1, 'serializeHash writes p1');
assert(back.p1Pct === 4 && back.mdePct === 15 && back.peeking === true && back.variants === 2,
  'parseHash roundtrip');
var hashed = L.parseHash('#' + hash);
assert(hashed.dailyVisitors === 800 && hashed.sided === 'two', 'parseHash strips leading #');
var fromEmpty = L.parseHash('');
assert(fromEmpty.p1Pct === 4 && fromEmpty.mdePct === 15, 'empty hash restores SaaS defaults');

var absHash = L.serializeHash({
  p1Pct: 5, mdePct: 1, mdeType: 'absolute', direction: 'down',
  dailyVisitors: 100, variants: 4, significance: 90, power: 90, sided: 'one', peeking: false
});
var absBack = L.parseHash(absHash);
assert(absBack.mdeType === 'absolute' && absBack.direction === 'down' && absBack.sided === 'one' &&
  absBack.significance === 90 && absBack.variants === 4, 'hash roundtrip absolute/one-sided/4-var');

// 14. drop direction
var drop = L.resolveRates({ p1: 0.10, mde: 0.20, mdeType: 'relative', direction: 'down' });
assert(drop.ok && Math.abs(drop.p2 - 0.08) < 1e-12, 'relative 20% drop of 10% → p2 0.08');

// 15. daily visitors 0 → days null, looks = 1
var noTraffic = L.planTest({ p1: 0.05, p2: 0.06, dailyVisitors: 0, peeking: true });
assert(noTraffic.ok && noTraffic.days == null && noTraffic.peeking.looks === 1,
  'zero traffic: days null, looks=1');
assert(L.planTest({ p1: 0.05, p2: 0.06, dailyVisitors: -3 }).ok === false, 'negative traffic errors');

// 16. mdeCurve produces increasing-MDE / decreasing-n points
var curve = L.mdeCurve({ p1Pct: 4, significance: 95, power: 80, variants: 2, dailyVisitors: 800 }, 0.05, 0.4, 7);
assert(curve.length >= 5, 'mdeCurve has points');
assert(curve[0].n > curve[curve.length - 1].n, 'larger relative MDE needs fewer visitors');

// 17. formatPlan is a non-empty English paragraph
var summary = L.formatPlan(planned);
assert(typeof summary === 'string' && summary.indexOf('per arm') !== -1 && summary.indexOf('peek') !== -1,
  'formatPlan mentions n and peeking');

console.log('');
console.log(passed + ' passed, ' + failed + ' failed');
if (failed) {
  console.log('failures:');
  failures.forEach(function (f) { console.log('  - ' + f); });
  process.exit(1);
}
console.log('sanity n per variant = ' + nSanity);
process.exit(0);
