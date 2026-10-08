"""Small in-memory example of quenched free-energy reconstruction."""

import numpy as np

from gauge_ite import TFIMSpec
from gauge_ite.studies import QuenchedFreeEnergyStudy


spec = TFIMSpec.rectangular(2, 3)
beta_values = np.linspace(0.0, 1.4, 29)
results = QuenchedFreeEnergyStudy(
    spec=spec,
    beta_values=beta_values,
    shots_per_basis=4096,
).run()

columns = [
    "beta",
    "ed_uniform_quenched_free_energy",
    "operational_ti_uniform_quenched_free_energy",
    "ed_annealed_free_energy",
]
print(results["summary"][columns].tail())
