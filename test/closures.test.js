// Run with: node --test
const test = require("node:test");
const assert = require("node:assert/strict");
const C = require("../closures.js");

const close = (actual, expected, relTol, label) => {
  const rel = Math.abs(actual - expected) / Math.abs(expected);
  assert.ok(rel <= relTol, `${label}: got ${actual}, expected ${expected} (rel. error ${rel.toExponential(2)})`);
};

// Fixed bed of 2 mm spheres in air, eps = 0.40, superficial velocity 0.5 m/s.
const bed = { d: 2e-3, rp: 2000, rf: 1.1, mu: 1.9e-5, k: 0.028, cp: 1007, w: 0.5, us: 0.5, eps: 0.4 };

test("Ergun terms match the hand calculation", () => {
  const r = C.compute(bed);
  // 150 mu (1-eps)^2 U / (eps^3 d^2) = 150 * 1.9e-5 * 0.36 * 0.5 / (0.064 * 4e-6)
  close(r.ergV, 2003.906, 1e-6, "viscous term");
  // 1.75 rho (1-eps) U^2 / (eps^3 d) = 1.75 * 1.1 * 0.6 * 0.25 / (0.064 * 0.002)
  close(r.ergI, 2255.859, 1e-6, "inertial term");
});

test("Di Felice drag summed over the bed balances eps x Ergun within 5%", () => {
  const r = C.compute(bed);
  close(r.chi, 3.0720, 1e-3, "Di Felice exponent chi");
  close(r.nF, 1641.8, 1e-3, "n x F_DF");
  assert.ok(Math.abs(r.diff) < 0.05, `difference ${r.diff}`);
});

test("Gunn Nusselt number matches the hand calculation", () => {
  const r = C.compute(bed);
  close(r.nuG, 17.557, 1e-3, "Gunn Nu");
});

test("Terminal velocity reduces to Stokes law for a 10 micron glass bead in water", () => {
  const t = C.terminalVelocity(10e-6, 2500, 998, 1.0e-3);
  const stokes = 9.81 * 10e-6 ** 2 * (2500 - 998) / (18 * 1.0e-3);
  assert.ok(t.re < 1e-2, `Re_t = ${t.re}`);
  close(t.v, stokes, 5e-3, "Stokes limit");
});

test("Terminal velocity reaches the Newton limit (Cd = 0.44) for a 20 mm steel sphere in air", () => {
  const t = C.terminalVelocity(20e-3, 7800, 1.2, 1.81e-5);
  const newton = Math.sqrt((4 / 3) * 9.81 * 20e-3 * (7800 - 1.2) / (1.2 * 0.44));
  assert.ok(t.re > 1000, `Re_t = ${t.re}`);
  close(t.v, newton, 1e-6, "Newton limit");
});

test("Drag balances buoyant weight at the terminal velocity (intermediate regime)", () => {
  const d = 0.5e-3, rp = 2500, rf = 998, mu = 1.0e-3;
  const t = C.terminalVelocity(d, rp, rf, mu);
  const drag = 0.5 * C.cdSchillerNaumann(rf * t.v * d / mu) * rf * t.v ** 2 * Math.PI * d * d / 4;
  const weight = Math.PI * d ** 3 / 6 * (rp - rf) * 9.81;
  assert.ok(t.re > 1 && t.re < 1000, `Re_t = ${t.re}`);
  close(drag, weight, 1e-9, "force balance");
});

test("Wen-Yu minimum fluidisation velocity for 0.5 mm glass beads in water", () => {
  const r = C.compute({ d: 0.5e-3, rp: 2500, rf: 998, mu: 1.0e-3, k: 0.6, cp: 4182, w: 0.05, us: 0.002, eps: 0.42 });
  close(r.ar, 1838.1, 1e-3, "Archimedes number");
  close(r.umf, 2.194e-3, 2e-3, "U_mf");
});
