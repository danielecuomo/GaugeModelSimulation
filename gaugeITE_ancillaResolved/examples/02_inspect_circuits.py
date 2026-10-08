"""Save logical and transpiled views of every circuit architecture."""

from gauge_ite import TFIMSpec
from gauge_ite.workflows import inspect_architectures


def main() -> None:
    # A 1x2 lattice keeps the diagrams readable. Change the dimensions after
    # first checking the small example.
    spec = TFIMSpec.rectangular(rows=1, cols=2)
    resources = inspect_architectures(
        spec,
        beta=0.75,
        steps=2,
        output_dir="artifacts/circuit_comparison",
    )
    print(resources.to_string(index=False))


if __name__ == "__main__":
    main()
