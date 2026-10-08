"""Smallest complete energy and gauge-class example."""

from gauge_ite import ITEProtocol, TFIMSpec
from gauge_ite.workflows import quickstart_report


def main() -> None:
    spec = TFIMSpec.rectangular(rows=2, cols=2)
    protocol = ITEProtocol(
        beta=0.75,
        steps=1,
        angle_mapping="exact_kraus",
        architecture="coherent_reuse",
        input_mode="purification",
    )
    report = quickstart_report(
        spec,
        protocol,
        beta_values=(0.0, 0.25, 0.5, 0.75),
    )

    print(f"Exact clean ground-state energy: {report['ground_energy']:.8f}")
    print(f"Circuit energy at beta=0.75:     {report['circuit_energy']:.8f}")
    print("\nClean thermodynamics")
    print(report["clean_thermodynamics"].round(6).to_string(index=False))
    print("\nGauge-class probabilities")
    print(report["class_probabilities"].round(6).to_string(index=False))


if __name__ == "__main__":
    main()
