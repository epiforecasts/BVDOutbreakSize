// Tests of the dashboard map's date logic, docs/src/public/zone_map/timeline.js.
// Run with `node --test docs/test/zone_map_timeline.test.js` (no packages needed).
"use strict";
const test = require("node:test");
const assert = require("node:assert/strict");
const T = require("../src/public/zone_map/timeline.js");

const day = T.dayOf;
const CUT = day("2026-10-20");
// Weekly counts ending on the cut-off and every seventh day before it.
function weekly(counts) {
  return counts.map((c, i) => ({ date: T.isoOf(CUT - 7 * (counts.length - 1 - i)), count: c }));
}
function series(map) {
  const ix = {};
  for (const k of Object.keys(map)) { ix[k] = T.indexSeries(map[k]); }
  return (s) => ix[s] || null;
}
const RT = [
  { date: "2026-10-01", median: 1.2, lower_90: 0.9, upper_90: 1.6 },
  { date: "2026-10-02", median: 1.1, lower_90: 0.8, upper_90: 1.5 },
  { date: "2026-10-20", median: 0.97, lower_90: 0.58, upper_90: 1.46 }
];

test("pointAt finds a day exactly or the last point on or before it", () => {
  const ix = T.indexSeries(RT);
  assert.equal(T.pointAt(ix, day("2026-10-02"), true).median, 1.1);
  assert.equal(T.pointAt(ix, day("2026-10-05"), true), null);
  assert.equal(T.pointAt(ix, day("2026-10-05")).median, 1.1);
  assert.equal(T.pointAt(ix, day("2026-09-30")), null);
  assert.equal(T.pointAt(ix, day("2026-11-01")).median, 0.97);
  assert.equal(T.pointAt(ix, day("2026-11-01"), true), null);
  assert.equal(T.pointAt(T.indexSeries([]), CUT), null);
});

test("weekEnd snaps to the cut-off's weekday on or before the date", () => {
  assert.equal(T.weekEnd(CUT, CUT), CUT);
  assert.equal(T.weekEnd(CUT - 1, CUT), CUT - 7);
  assert.equal(T.weekEnd(CUT - 7, CUT), CUT - 7);
  assert.equal(T.weekEnd(CUT - 8, CUT), CUT - 14);
});

test("weekly sums give the past 7, 14 and 28 days and the total to date", () => {
  // Beni's weekly cases, 1 September to 20 October, and its 553 to date.
  const w = weekly([28, 39, 56, 51, 56, 61, 73, 71]);
  assert.equal(T.weeksSum(w, CUT, 1, CUT), 71);
  assert.equal(T.weeksSum(w, CUT, 2, CUT), 144);
  assert.equal(T.weeksSum(w, CUT, 4, CUT), 261);
  assert.equal(T.weeksSum(w, CUT - 7, 2, CUT), 134);
  // Four weeks back from 15 September reach before 1 September.
  assert.equal(T.weeksSum(w, day("2026-09-15"), 4, CUT), null);
  assert.equal(T.totalAt(553, w, CUT, CUT), 553);
  assert.equal(T.totalAt(553, w, CUT - 7, CUT), 482);
  assert.equal(T.totalAt(553, w, day("2026-08-25"), CUT), 553 - 435);
  assert.equal(T.totalAt(553, w, day("2026-08-18"), CUT), null);
  // A series that stops short of the cut-off cannot be read back from it.
  assert.equal(T.totalAt(553, w.slice(0, -1), CUT - 7, CUT), null);
});

test("rowAt reads the reproduction number and counts from the series", () => {
  const r = { zone: "beni", as_of: "2026-10-20", R_T_median: 0.97, R_T_lower: 0.58,
              R_T_upper: 1.46, cases: 553, deaths: 365, cases_last_7: 71,
              cases_last_14: 144, cases_last_28: 261, last_case_date: "2026-10-20",
              forecast_median: 107, p_ge_10: 1, p_rt_above_one: 0.45 };
  const get = series({ rt: RT, cases: weekly([28, 39, 56, 51, 56, 61, 73, 71]),
                       deaths: weekly([20, 33, 40, 26, 25, 39, 54, 43]) });
  const o = T.rowAt(r, "zone", day("2026-10-15"), get, undefined);
  assert.equal(o.R_T_median, null);
  assert.equal(T.rowAt(r, "zone", day("2026-10-02"), get).R_T_upper, 1.5);
  assert.equal(o.cases, 482);
  assert.equal(o.deaths, 322);
  assert.equal(o.cases_last_14, 134);
  assert.equal(o.last_case_date, null);
  // With no saved cut-off, the estimates without a series stay at the cut-off.
  assert.equal(o.forecast_median, 107);
  assert.equal(o.p_rt_above_one, 0.45);
  assert.equal(r.cases, 553, "the row itself is not changed");
});

test("rowAt takes the estimates without a series from a saved cut-off", () => {
  const r = { zone: "beni", as_of: "2026-10-20", forecast_median: 107, p_rt_above_one: 0.45,
              p_ge_10: 1, cfr_median: 0.5 };
  const get = series({});
  const saved = { forecast_median: 86, p_rt_above_one: 0.6, p_ge_10: 0.9, cfr_median: null };
  const o = T.rowAt(r, "zone", CUT - 14, get, saved);
  assert.equal(o.forecast_median, 86);
  assert.equal(o.p_rt_above_one, 0.6, "the chance R > 1 is the saved cut-off's, not a median");
  assert.equal(o.p_ge_10, 0.9);
  assert.equal(o.cfr_median, null);
  // Saved cut-offs exist but none is this early: no estimate.
  const none = T.rowAt(r, "zone", CUT - 14, get, null);
  assert.equal(none.forecast_median, null);
  assert.equal(none.p_rt_above_one, null);
});

test("rowAt reads a province's reported treatment-centre figures", () => {
  const r = { province: "ituri", patch: "Ituri", pooled: 0, as_of: "2026-10-20", cases: 10 };
  const get = series({
    isolation: [{ date: "2026-09-10", count: 435 }, { date: "2026-10-01", count: 400 }],
    beds: [{ date: "2026-09-01", count: 1000 }, { date: "2026-09-20", count: 1015 }]
  });
  const o = T.rowAt(r, "province", day("2026-09-15"), get);
  assert.equal(o.isolation_reported, 435);
  assert.equal(o.isolation_reported_date, "2026-09-10");
  assert.equal(o.beds_reported, 1000);
  assert.equal(o.bed_use_reported, 0.435);
  // The beds were last counted 30 days before the patients: no share.
  const old = T.rowAt(r, "province", day("2026-10-05"), series({
    isolation: [{ date: "2026-10-01", count: 400 }], beds: [{ date: "2026-09-01", count: 1000 }]
  }));
  assert.equal(old.bed_use_reported, null);
  assert.equal(T.rowAt(r, "province", day("2026-08-01"), get).isolation_reported, null);
});

test("rowAt leaves a pooled province's counts blank before the cut-off", () => {
  const r = { province: "tshopo", patch: "Other", pooled: 1, as_of: "2026-10-20", cases: 76 };
  const get = series({ cases: weekly([1, 2, 3]) });
  assert.equal(T.rowAt(r, "province", CUT - 7, get).cases, null);
  assert.equal(T.rowAt({ ...r, pooled: 0 }, "province", CUT - 7, get).cases, 73);
});

test("vintageAt picks the latest saved cut-off on or before the date", () => {
  const v = ["2026-09-29", "2026-10-06", "2026-10-13"];
  assert.equal(T.vintageAt(v, day("2026-10-08")), "2026-10-06");
  assert.equal(T.vintageAt(v, day("2026-10-13")), "2026-10-13");
  assert.equal(T.vintageAt(v, day("2026-09-20")), null);
  assert.equal(T.vintageAt([], CUT), null);
});

test("the slider's span and clamp keep the date inside it", () => {
  const sp = T.span([CUT - 100, NaN, CUT - 40], CUT);
  assert.deepEqual(sp, { first: CUT - 100, end: CUT });
  assert.equal(T.span([CUT], CUT), null);
  assert.equal(T.span([], CUT), null);
  assert.equal(T.clampDay(CUT - 200, sp), CUT - 100);
  assert.equal(T.clampDay(CUT - 5, sp), CUT - 5);
  assert.equal(T.clampDay(CUT, sp), null);
  assert.equal(T.clampDay(CUT + 3, sp), null);
  assert.equal(T.clampDay(null, sp), null);
  assert.equal(T.clampDay(CUT - 5, null), null);
});

test("play steps a week on the cut-off's weekday and stops at the cut-off", () => {
  const sp = { first: CUT - 23, end: CUT };
  let d = T.playStart(sp);
  assert.equal(d, CUT - 21);
  const seen = [d];
  while (d !== null) { d = T.playNext(d, sp); seen.push(d); }
  assert.deepEqual(seen, [CUT - 21, CUT - 14, CUT - 7, null]);
  assert.equal(T.playNext(null, sp), null);
});

test("the range keeps the chosen weeks and the axis fits them", () => {
  const pts = [
    { date: "2026-06-01", lower_90: 1.5, upper_90: 3.4 },
    { date: "2026-09-30", lower_90: 0.8, upper_90: 1.3 },
    { date: "2026-10-20", lower_90: 0.58, upper_90: 1.46 }
  ];
  assert.equal(T.rangeFrom(CUT, 0), null);
  assert.equal(T.rangeFrom(CUT, 4), CUT - 28);
  assert.equal(T.rangeFrom(NaN, 4), null);
  assert.deepEqual(T.inRange(pts, null), pts);
  assert.deepEqual(T.inRange(pts, T.rangeFrom(CUT, 4)).map((p) => p.date),
    ["2026-09-30", "2026-10-20"]);
  assert.deepEqual(T.inRange(undefined, CUT), []);
  assert.deepEqual(T.rtAxis(pts), { bot: 0.5, top: 3.5, step: 0.5 });
  assert.deepEqual(T.rtAxis(pts.slice(1)), { bot: 0.5, top: 1.5, step: 0.25 });
  assert.deepEqual(T.rtAxis([{ lower_90: 0.2, upper_90: 9 }]), { bot: 0, top: 5, step: 1 });
  assert.deepEqual(T.rtAxis([]), { bot: 0.5, top: 1.5, step: 0.25 });
});
