/*
 * CFD-DEM Closure Calculator: physics core.
 * Particle-fluid closures for a sphere in a fluid and for a particle bed. SI units throughout.
 * Works in the browser (window.Closures) and in Node (module.exports).
 */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.Closures = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  var G = 9.81; // m/s^2

  /** Schiller & Naumann (1933): Cd = 24/Re (1 + 0.15 Re^0.687) for Re < 1000, then 0.44. */
  function cdSchillerNaumann(re) {
    if (!(re > 0)) return Infinity;
    return re < 1000 ? (24 / re) * (1 + 0.15 * Math.pow(re, 0.687)) : 0.44;
  }

  /** Dallavalle (1948): Cd = (0.63 + 4.8/sqrt(Re))^2, the single-particle law used by Di Felice (1994). */
  function cdDallavalle(re) {
    if (!(re > 0)) return Infinity;
    var a = 0.63 + 4.8 / Math.sqrt(re);
    return a * a;
  }

  /**
   * Terminal velocity of an isolated sphere from the force balance
   * (pi d^3/6)|rho_p - rho_f| g = 0.5 Cd(Re) rho_f v^2 (pi d^2/4), with Schiller-Naumann drag.
   * Cd(Re) v^2 increases monotonically with v, so bisection in log space is robust in every regime.
   */
  function terminalVelocity(d, rhoP, rhoF, mu) {
    var drho = Math.abs(rhoP - rhoF);
    if (drho === 0) return { v: 0, direction: "neutral", re: 0, stokes: 0 };
    var K = (4 / 3) * G * d * drho / rhoF;
    var lo = 1e-12, hi = 1e4;
    for (var i = 0; i < 200; i++) {
      var mid = Math.sqrt(lo * hi);
      if (cdSchillerNaumann(rhoF * mid * d / mu) * mid * mid - K > 0) hi = mid; else lo = mid;
    }
    var v = Math.sqrt(lo * hi);
    return {
      v: v,
      direction: rhoP > rhoF ? "settles" : "rises",
      re: rhoF * v * d / mu,
      stokes: G * d * d * drho / (18 * mu)
    };
  }

  /**
   * All readouts for one set of inputs.
   * p = { d [m], rp [kg/m3], rf [kg/m3], mu [Pa s], k [W/(m K)], cp [J/(kg K)],
   *       w [m/s] slip velocity of the isolated particle, us [m/s] superficial velocity, eps [-] bed voidage }
   */
  function compute(p) {
    var d = p.d, rp = p.rp, rf = p.rf, mu = p.mu, k = p.k, cp = p.cp, w = p.w, us = p.us, eps = p.eps;
    var r = {};
    r.pr = cp * mu / k;

    // Isolated particle
    r.re = rf * w * d / mu;
    r.cd = cdSchillerNaumann(r.re);
    var area = Math.PI * d * d / 4;
    r.fd = 0.5 * r.cd * rf * w * w * area;
    if (!isFinite(r.fd)) r.fd = 3 * Math.PI * mu * d * w; // Stokes drag as the w -> 0 limit
    r.tau = rp * d * d / (18 * mu);
    r.ar = rf * Math.abs(rp - rf) * G * d * d * d / (mu * mu);
    r.term = terminalVelocity(d, rp, rf, mu);
    r.nuRM = 2 + 0.6 * Math.sqrt(r.re) * Math.cbrt(r.pr); // Ranz & Marshall (1952)
    r.hRM = r.nuRM * k / d;

    // Particle bed
    r.res = rf * us * d / mu;
    r.ergV = 150 * mu * Math.pow(1 - eps, 2) * us / (Math.pow(eps, 3) * d * d); // Ergun (1952), viscous
    r.ergI = 1.75 * rf * (1 - eps) * us * us / (Math.pow(eps, 3) * d);           // Ergun (1952), inertial
    r.erg = r.ergV + r.ergI;
    r.remf = Math.sqrt(33.7 * 33.7 + 0.0408 * r.ar) - 33.7;                          // Wen & Yu (1966)
    r.umf = r.remf * mu / (rf * d);
    r.bedWeight = (1 - eps) * (rp - rf) * G;
    r.chi = r.res > 0 ? 3.7 - 0.65 * Math.exp(-Math.pow(1.5 - Math.log10(r.res), 2) / 2) : NaN; // Di Felice (1994)
    r.cdb = cdDallavalle(r.res);
    r.fDF = 0.5 * r.cdb * rf * area * us * us * Math.pow(eps, -r.chi);
    r.n = 6 * (1 - eps) / (Math.PI * Math.pow(d, 3));
    r.nF = r.n * r.fDF;
    // The drag summed over a fixed bed balances eps * (pressure gradient): the particles also carry the
    // pressure-gradient (buoyancy) force, which accounts for the remaining (1 - eps) fraction.
    r.epsErg = eps * r.erg;
    r.diff = r.epsErg > 0 ? (r.nF - r.epsErg) / r.epsErg : NaN;
    r.nuG = (7 - 10 * eps + 5 * eps * eps) * (1 + 0.7 * Math.pow(r.res, 0.2) * Math.cbrt(r.pr)) +
            (1.33 - 2.4 * eps + 1.2 * eps * eps) * Math.pow(r.res, 0.7) * Math.cbrt(r.pr); // Gunn (1978)
    r.hG = r.nuG * k / d;
    return r;
  }

  return { G: G, cdSchillerNaumann: cdSchillerNaumann, cdDallavalle: cdDallavalle, terminalVelocity: terminalVelocity, compute: compute };
});
