// The dashboard map's date logic: the rows an area shows on a date before
// the cut-off, the date slider's span and steps, and the detail charts'
// period and axis. It has no DOM, so `docs/test/zone_map_timeline.test.js`
// runs it under Node. The page loads it before its own script.
(function (root) {
  "use strict";
  function num(x) { return typeof x === "number" && isFinite(x); }
  // ISO dates and day numbers (days since 1970-01-01, UTC).
  function dayOf(iso) { return Date.parse(iso + "T00:00:00Z") / 864e5; }
  function isoOf(d) { return new Date(d * 864e5).toISOString().slice(0, 10); }

  // A date-sorted series with a slot per day from its first point to its
  // last, so a point on a day, or the last point on or before it, is one
  // lookup.
  function indexSeries(pts) {
    var ix = { pts: pts, first: NaN, at: [], upto: [] };
    if (!pts.length) { return ix; }
    var days = pts.map(function (p) { return dayOf(p.date); });
    ix.first = days[0];
    for (var d = days[0], j = 0; d <= days[days.length - 1]; d++) {
      while (j + 1 < days.length && days[j + 1] <= d) { j++; }
      ix.upto.push(j);
      ix.at.push(days[j] === d ? pts[j] : null);
    }
    return ix;
  }
  var EMPTY = indexSeries([]);
  // The point of `ix` on day `d` with `exact`, otherwise the last on or
  // before it; null when there is none.
  function pointAt(ix, d, exact) {
    if (!ix || !ix.pts.length || !(d >= ix.first)) { return null; }
    var i = d - ix.first;
    if (i >= ix.upto.length) { return exact ? null : ix.pts[ix.pts.length - 1]; }
    return exact ? ix.at[i] : ix.pts[ix.upto[i]];
  }

  // The weekly counts end on the cut-off `end` and every seventh day
  // before it. `weekEnd` is the last such day on or before `d`.
  function weekEnd(d, end) { return end - 7 * Math.ceil((end - d) / 7); }
  function weeklyCovers(pts, end) {
    return pts.length > 0 && pts[pts.length - 1].date === isoOf(end);
  }
  // The count over the `k` weeks ending on week end `w`, or null when the
  // series does not reach that far back.
  function weeksSum(pts, w, k, end) {
    if (!weeklyCovers(pts, end) || w - 7 * (k - 1) < dayOf(pts[0].date)) { return null; }
    return pts.reduce(function (s, p) {
      var d = dayOf(p.date);
      return d <= w && d > w - 7 * k ? s + p.count : s;
    }, 0);
  }
  // The total to week end `w`: the total at the cut-off less the weeks
  // after `w`, or null before the series starts.
  function totalAt(total, pts, w, end) {
    if (!num(total) || !weeklyCovers(pts, end) || w < dayOf(pts[0].date) - 7) { return null; }
    return total - pts.reduce(function (s, p) {
      return dayOf(p.date) > w ? s + p.count : s;
    }, 0);
  }

  // The estimates with no time series, taken whole from a saved cut-off.
  var SAVED_FIELDS = ["forecast", "deaths_forecast", "share", "cfr", "ascertainment",
    "isolation", "beds", "bed_use", "isolation_forecast", "beds_forecast",
    "bed_use_forecast"].reduce(function (a, p) {
      return a.concat([p + "_median", p + "_lower", p + "_upper"]);
    }, ["p_ge_1", "p_ge_5", "p_ge_10", "p_ge_20", "p_rt_above_one"]);
  // The latest of the ascending ISO dates `vintages` on or before day `d`.
  function vintageAt(vintages, d) {
    var v = null;
    vintages.forEach(function (x) { if (dayOf(x) <= d) { v = x; } });
    return v;
  }

  // Area row `r` at level `lvl` on day `d`. `get(s)` gives the area's
  // indexed series `s`. The reproduction number, the counts and the
  // reported treatment-centre figures come from the series. `saved` is the
  // area's row at the latest saved cut-off on or before `d`, whose
  // estimates with no series replace the row's; null blanks them and
  // undefined keeps the cut-off's.
  function rowAt(r, lvl, d, get, saved) {
    var o = Object.assign({}, r), end = dayOf(r.as_of), w = weekEnd(d, end);
    var p = pointAt(get("rt"), d, true) || {};
    o.R_T_median = num(p.median) ? p.median : null;
    o.R_T_lower = num(p.lower_90) ? p.lower_90 : null;
    o.R_T_upper = num(p.upper_90) ? p.upper_90 : null;
    // A pooled province's counts are its own, but its series are the pool's.
    var own = !(lvl === "province" && r.pooled === 1);
    var cs = own ? (get("cases") || EMPTY).pts : [];
    o.cases = totalAt(r.cases, cs, w, end);
    o.deaths = totalAt(r.deaths, own ? (get("deaths") || EMPTY).pts : [], w, end);
    [7, 14, 28].forEach(function (n) {
      if ("cases_last_" + n in r) { o["cases_last_" + n] = weeksSum(cs, w, n / 7, end); }
    });
    o.last_case_date = null;
    if (lvl !== "zone") {
      // As at the cut-off: the patients over the last bed count on or
      // before their report and at most 14 days older.
      var iso = pointAt(get("isolation"), d), bed = pointAt(get("beds"), d);
      var bu = iso && pointAt(get("beds"), dayOf(iso.date));
      o.isolation_reported = iso ? iso.count : null;
      o.isolation_reported_date = iso ? iso.date : null;
      o.beds_reported = bed ? bed.count : null;
      o.beds_reported_date = bed ? bed.date : null;
      o.bed_use_reported = bu && bu.count > 0 &&
        dayOf(iso.date) - dayOf(bu.date) <= 14 ? iso.count / bu.count : null;
    }
    if (saved !== undefined) {
      SAVED_FIELDS.forEach(function (k) { o[k] = saved && num(saved[k]) ? saved[k] : null; });
    }
    return o;
  }

  // The slider's span: from the earliest of the day numbers `firsts` to the
  // cut-off `end`, or null when nothing precedes the cut-off.
  function span(firsts, end) {
    var first = Math.min.apply(null, firsts.filter(num).concat([Infinity]));
    return num(end) && first < end ? { first: first, end: end } : null;
  }
  // A day kept inside `sp`, or null (the cut-off) at or past its end.
  function clampDay(d, sp) {
    if (d === null || !sp || d >= sp.end) { return null; }
    return Math.max(sp.first, d);
  }
  // Play starts on the first day of the span on the cut-off's weekday and
  // steps a week, reaching the cut-off as null.
  function playStart(sp) { return sp.end - 7 * Math.floor((sp.end - sp.first) / 7); }
  function playNext(d, sp) { return d === null ? null : clampDay(d + 7, sp); }

  // The detail charts' period: the `weeks` weeks up to the cut-off day
  // `asOf`, or everything when `weeks` is 0.
  function rangeFrom(asOf, weeks) { return weeks && num(asOf) ? asOf - 7 * weeks : null; }
  function inRange(pts, from) {
    return from === null ? pts || [] :
      (pts || []).filter(function (p) { return dayOf(p.date) > from; });
  }
  // The reproduction-number chart's y-axis over the points shown: their
  // 90% intervals on half steps, always covering 0.5 to 1.5 and never
  // beyond 0 to 5, with ticks every 1, 0.5 or 0.25.
  function rtAxis(pts) {
    function ends(k) { return pts.map(function (p) { return p[k]; }).filter(num); }
    var top = Math.max(1.5, Math.ceil(Math.min(5, Math.max.apply(null, ends("upper_90"))) * 2) / 2);
    var bot = Math.min(0.5, Math.floor(Math.max(0, Math.min.apply(null, ends("lower_90"))) * 2) / 2);
    var step = top - bot > 3 ? 1 : top - bot > 1.5 ? 0.5 : 0.25;
    return { bot: bot, top: top, step: step };
  }

  // Where an observed count `obs` falls in its posterior predictive `q`
  // (median, lower_90/upper_90, lower_50/upper_50), as a probability
  // integral transform read off those quantiles: 0.025 below the 90%
  // interval, 0.15 below the 50% interval, 0.5 inside it, 0.85 above it and
  // 0.975 above the 90% interval. A count on a bound is inside it, so a
  // count the predictive puts all its mass on reads 0.5. Null without both.
  function fitClass(obs, q) {
    if (!num(obs) || !q || !num(q.lower_90) || !num(q.upper_90) ||
        !num(q.lower_50) || !num(q.upper_50)) { return null; }
    if (obs < q.lower_90) { return 0.025; }
    if (obs < q.lower_50) { return 0.15; }
    if (obs > q.upper_90) { return 0.975; }
    if (obs > q.upper_50) { return 0.85; }
    return 0.5;
  }
  // An area's observed weekly count and its predictive for the week ending
  // on week end `w`, from the indexed series `obsIx` and `fitIx`, or null
  // when either lacks that week.
  function fitAt(obsIx, fitIx, w) {
    var o = pointAt(obsIx, w, true), q = pointAt(fitIx, w, true);
    if (!o || !q || !num(o.count)) { return null; }
    return { date: q.date, obs: o.count, median: q.median, lower_90: q.lower_90,
             upper_90: q.upper_90, lower_50: q.lower_50, upper_50: q.upper_50,
             pit: fitClass(o.count, q) };
  }

  var T = {
    dayOf: dayOf, isoOf: isoOf, indexSeries: indexSeries, pointAt: pointAt,
    weekEnd: weekEnd, weeksSum: weeksSum, totalAt: totalAt,
    SAVED_FIELDS: SAVED_FIELDS, vintageAt: vintageAt, rowAt: rowAt,
    span: span, clampDay: clampDay, playStart: playStart, playNext: playNext,
    rangeFrom: rangeFrom, inRange: inRange, rtAxis: rtAxis,
    fitClass: fitClass, fitAt: fitAt
  };
  if (typeof module === "object" && module.exports) { module.exports = T; }
  else { root.BVDTime = T; }
})(this);
