"""Run a fresh trajectory and inspect its branch and class tables."""

from gauge_ite import GaugeITEExperiment, ITEProtocol, TFIMSpec


def main() -> None:
    spec = TFIMSpec.rectangular(rows=2, cols=2)
    protocol = ITEProtocol(
        beta=0.75,
        steps=2,
        architecture="fresh_trajectory",
        angle_mapping="exact_kraus",
        input_mode="purification",
    )
    result = GaugeITEExperiment(spec, protocol).run("trajectory_records")

    columns = [
        "ancilla_bitstring",
        "class_ids_by_step",
        "plaquette_fluxes_by_step",
        "persistent_signs",
        "count",
        "probability",
    ]
    print("Most common trajectory records")
    print(result.analysis["branches"][columns].head(12).to_string(index=False))
    print("\nFinal gauge classes")
    print(result.analysis["classes"].to_string(index=False))
    print("\nRun summary")
    print(result.analysis["summary"])


if __name__ == "__main__":
    main()
