"""Water tank level control: P vs. PI vs. PID on a nonlinear tank model.

Model (Torricelli outflow):
    A * dh/dt = k_pump * u - c * sqrt(h) * leak

    h     tank level [m]
    u     pump command, 0..1 (saturated)
    leak  outflow multiplier (1.0 = normal, >1.0 = leak / open drain valve)

Scenario:
    - Start in steady state at h = 0.20 m
    - t = 50 s:  setpoint step 0.20 m -> 0.40 m
    - t = 900 s: disturbance, outflow +20 %

Run:
    python src/tank_control.py
Outputs:
    results/step_response.png
    results/metrics.csv
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # no GUI needed, works on every machine
import matplotlib.pyplot as plt
import numpy as np

# --------------------------------------------------------------------------
# Plant
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Tank:
    """Physical parameters of the tank."""

    area: float = 0.5  # cross-section A [m^2]
    pump_gain: float = 0.004  # inflow at u = 1 [m^3/s]
    outflow_coeff: float = 0.0035  # c in q_out = c * sqrt(h) [m^2.5/s]
    h_max: float = 1.0  # overflow level [m]

    def derivative(self, h: float, u: float, leak: float = 1.0) -> float:
        """dh/dt for level h, pump command u and outflow multiplier leak."""
        q_in = self.pump_gain * u
        q_out = self.outflow_coeff * np.sqrt(max(h, 0.0)) * leak
        return (q_in - q_out) / self.area

    def equilibrium_input(self, h: float) -> float:
        """Pump command that holds level h in steady state."""
        return self.outflow_coeff * np.sqrt(h) / self.pump_gain


# --------------------------------------------------------------------------
# Controller
# --------------------------------------------------------------------------


class PID:
    """PID controller with output limits and anti-windup.

    - Derivative acts on the measurement (no kick at setpoint steps).
    - Anti-windup: the integrator is frozen while the output is saturated
      and the error would push it further into saturation.
    """

    def __init__(
        self,
        kp: float,
        ki: float = 0.0,
        kd: float = 0.0,
        u_min: float = 0.0,
        u_max: float = 1.0,
        bias: float = 0.0,
    ) -> None:
        self.kp, self.ki, self.kd = kp, ki, kd
        self.u_min, self.u_max = u_min, u_max
        self.bias = bias
        self.integral = 0.0
        self._prev_measurement: float | None = None

    def step(self, setpoint: float, measurement: float, dt: float) -> float:
        error = setpoint - measurement

        if self._prev_measurement is None:
            d_term = 0.0
        else:
            d_term = -self.kd * (measurement - self._prev_measurement) / dt
        self._prev_measurement = measurement

        u_unsat = self.bias + self.kp * error + self.integral + d_term
        u = min(max(u_unsat, self.u_min), self.u_max)

        pushing_up = u_unsat > self.u_max and error > 0
        pushing_down = u_unsat < self.u_min and error < 0
        if not (pushing_up or pushing_down):
            self.integral += self.ki * error * dt
        return u


# --------------------------------------------------------------------------
# Simulation
# --------------------------------------------------------------------------

T_END = 1500.0  # simulation length [s]
DT = 0.5  # time step [s]
H_START = 0.20  # initial level [m]
H_TARGET = 0.40  # setpoint after the step [m]
T_STEP = 50.0  # time of the setpoint step [s]
T_DISTURBANCE = 900.0  # time of the outflow disturbance [s]
LEAK_FACTOR = 1.20  # outflow multiplier after the disturbance


def simulate(tank: Tank, controller: PID) -> dict[str, np.ndarray]:
    """Euler simulation of tank + controller. Returns time series."""
    n = int(T_END / DT) + 1
    t = np.linspace(0.0, T_END, n)
    h = np.zeros(n)
    u = np.zeros(n)
    sp = np.zeros(n)
    h[0] = H_START

    for k in range(n - 1):
        sp[k] = H_START if t[k] < T_STEP else H_TARGET
        leak = 1.0 if t[k] < T_DISTURBANCE else LEAK_FACTOR
        u[k] = controller.step(sp[k], h[k], DT)
        h[k + 1] = h[k] + DT * tank.derivative(h[k], u[k], leak)
        h[k + 1] = min(max(h[k + 1], 0.0), tank.h_max)

    sp[-1] = sp[-2]
    u[-1] = u[-2]
    return {"t": t, "h": h, "u": u, "sp": sp}


# --------------------------------------------------------------------------
# Metrics
# --------------------------------------------------------------------------


def step_metrics(res: dict[str, np.ndarray]) -> dict[str, float]:
    """Overshoot, rise time, settling time and steady-state error of the step."""
    t, h = res["t"], res["h"]
    mask = (t >= T_STEP) & (t < T_DISTURBANCE)
    ts, hs = t[mask], h[mask]
    step = H_TARGET - H_START

    overshoot = max(0.0, (hs.max() - H_TARGET) / step * 100.0)

    def first_crossing(level: float) -> float:
        idx = np.argmax(hs >= level)
        return float(ts[idx] - T_STEP) if hs[idx] >= level else float("nan")

    rise_time = first_crossing(H_START + 0.9 * step) - first_crossing(
        H_START + 0.1 * step
    )

    band = 0.02 * step
    outside = np.where(np.abs(hs - H_TARGET) > band)[0]
    if len(outside) == 0:
        settling = 0.0
    elif outside[-1] == len(hs) - 1:
        settling = float("nan")  # never settles inside the 2 % band
    else:
        settling = float(ts[outside[-1] + 1] - T_STEP)

    last = (t >= T_DISTURBANCE - 50.0) & (t < T_DISTURBANCE)
    sse = float(H_TARGET - h[last].mean())

    # Largest level drop caused by the disturbance
    after = t >= T_DISTURBANCE
    dip = float(H_TARGET - h[after].min())

    return {
        "overshoot_percent": overshoot,
        "rise_time_s": rise_time,
        "settling_time_s": settling,
        "steady_state_error_m": sse,
        "disturbance_dip_m": max(dip, 0.0),
    }


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------


def make_plot(results: dict[str, dict[str, np.ndarray]], path: Path) -> None:
    colors = {"P": "#d95f02", "PI": "#1b9e77", "PID": "#7570b3"}
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(9, 6.5), sharex=True, gridspec_kw={"height_ratios": [3, 1.6]}
    )

    first = next(iter(results.values()))
    ax1.plot(first["t"], first["sp"], "k--", lw=1.4, label="Setpoint")
    for name, res in results.items():
        ax1.plot(res["t"], res["h"], lw=2, color=colors[name], label=name)
        ax2.plot(res["t"], res["u"], lw=1.8, color=colors[name], label=name)

    for ax in (ax1, ax2):
        ax.axvline(T_STEP, color="grey", lw=0.8, ls=":")
        ax.axvline(T_DISTURBANCE, color="grey", lw=0.8, ls=":")
        ax.grid(alpha=0.3)

    ax1.annotate("Setpoint step", (T_STEP, 0.43), xytext=(T_STEP + 15, 0.44), fontsize=9)
    ax1.annotate(
        "Outflow +20 %", (T_DISTURBANCE, 0.43), xytext=(T_DISTURBANCE + 15, 0.44), fontsize=9
    )
    ax1.set_ylabel("Tank level h [m]")
    ax1.set_title("Water tank level control: P vs. PI vs. PID")
    ax1.set_ylim(0.15, 0.47)
    ax1.legend(loc="lower right", ncol=4)

    ax2.set_ylabel("Pump command u [-]")
    ax2.set_xlabel("Time t [s]")
    ax2.set_ylim(-0.05, 1.05)

    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def save_metrics(metrics: dict[str, dict[str, float]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "controller",
        "overshoot_percent",
        "rise_time_s",
        "settling_time_s",
        "steady_state_error_m",
        "disturbance_dip_m",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for name, m in metrics.items():
            writer.writerow({"controller": name, **{k: round(v, 4) for k, v in m.items()}})


def main() -> None:
    tank = Tank()
    u0 = tank.equilibrium_input(H_START)  # start exactly in steady state

    controllers = {
        "P": PID(kp=1.5, bias=u0),
        "PI": PID(kp=1.5, ki=1.5 / 150.0, bias=0.0),
        "PID": PID(kp=1.5, ki=1.5 / 150.0, kd=30.0, bias=0.0),
    }
    # PI/PID start with the integrator preloaded so they also begin in steady state
    controllers["PI"].integral = u0
    controllers["PID"].integral = u0

    results = {name: simulate(tank, ctrl) for name, ctrl in controllers.items()}
    metrics = {name: step_metrics(res) for name, res in results.items()}

    root = Path(__file__).resolve().parent.parent
    make_plot(results, root / "results" / "step_response.png")
    save_metrics(metrics, root / "results" / "metrics.csv")

    header = f"{'Ctrl':<5}{'Overshoot %':>13}{'Rise [s]':>10}{'Settle [s]':>12}{'SSE [m]':>10}{'Dip [m]':>10}"
    print(header)
    print("-" * len(header))
    for name, m in metrics.items():
        print(
            f"{name:<5}{m['overshoot_percent']:>13.2f}{m['rise_time_s']:>10.1f}"
            f"{m['settling_time_s']:>12.1f}{m['steady_state_error_m']:>10.4f}"
            f"{m['disturbance_dip_m']:>10.4f}"
        )
    print(f"\nSaved: {root / 'results' / 'step_response.png'}")


if __name__ == "__main__":
    main()
