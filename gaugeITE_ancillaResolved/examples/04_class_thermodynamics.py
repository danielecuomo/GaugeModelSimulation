"""Calculate exact thermodynamics for every gauge-flux class."""

from gauge_ite import TFIMSpec
from gauge_ite.thermodynamics import exact_class_thermodynamics


def main() -> None:
    spec = TFIMSpec.rectangular(rows=2, cols=3)
    table = exact_class_thermodynamics(
        spec,
        beta_values=(0.0, 0.25, 0.5, 0.75, 1.0),
    )
    columns = [
        "beta",
        "class_id",
        "flux_label",
        "ground_energy",
        "thermal_energy",
        "annealed_class_weight",
    ]
    print(table[columns].round(6).to_string(index=False))


if __name__ == "__main__":
    main()
