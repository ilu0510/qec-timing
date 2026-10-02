"""Statistical and configuration guards for the full-circuit timing experiment."""

import importlib.util
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "stage5_validation", Path(__file__).resolve().parents[1] / "scripts/stage5_validation.py")
stage5 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stage5)


def test_confirmation_distinguishes_improvement_from_sampling_uncertainty():
    result = stage5.ratio_summary(dict(failures=100, shots=1000),
                                  dict(failures=10, shots=1000), confidence=0.99)
    assert result["status"] == "improved"
    assert result["gamma"] == 10
    assert result["gamma_low"] > 1
    tied = stage5.ratio_summary(dict(failures=100, shots=1000),
                                dict(failures=100, shots=1000), confidence=0.99)
    assert tied["status"] == "inconclusive"
    assert tied["gamma_low"] < 1 < tied["gamma_high"]


def test_zero_failures_do_not_prove_zero_error_probability():
    result = stage5.ratio_summary(dict(failures=0, shots=1000),
                                  dict(failures=0, shots=1000), confidence=0.99)
    assert result == dict(gamma=None, gamma_low=0.0, gamma_high=None, status="inconclusive")


@pytest.mark.parametrize("argv", [
    ["--rounds", "5", "5", "10"], ["--T", "nan"], ["--lam", "nan"],
    ["--distances", "4"], ["--q-d", "0"], ["--confidence", "1"],
])
def test_rejects_invalid_experiments(argv):
    with pytest.raises(SystemExit):
        stage5.parse_args(argv)
