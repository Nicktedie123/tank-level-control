"""Basic sanity tests. Run with:  python -m pytest"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from tank_control import PID, Tank  # noqa: E402


def test_equilibrium_is_stationary():
    tank = Tank()
    u = tank.equilibrium_input(0.3)
    assert abs(tank.derivative(0.3, u)) < 1e-12


def test_more_pump_raises_level():
    tank = Tank()
    u_eq = tank.equilibrium_input(0.3)
    assert tank.derivative(0.3, u_eq + 0.1) > 0


def test_controller_output_is_limited():
    ctrl = PID(kp=100.0)
    assert ctrl.step(1.0, 0.0, 0.5) == 1.0
    assert ctrl.step(0.0, 1.0, 0.5) == 0.0


def test_integrator_does_not_wind_up():
    ctrl = PID(kp=1.0, ki=1.0)
    for _ in range(1000):
        ctrl.step(10.0, 0.0, 0.5)  # output saturated the whole time
    assert ctrl.integral < 1.5
