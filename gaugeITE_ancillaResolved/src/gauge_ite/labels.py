"""Central publication vocabulary used by tables, plots, and notebooks."""

from __future__ import annotations

from types import MappingProxyType


PUBLICATION_LABELS = MappingProxyType(
    {
        "uniform_class_prior": "uniform class prior",
        "ed_annealed_class_weights": "ED annealed class weights",
        "noiseless_circuit_result": "noiseless circuit result",
        "finite_shot_circuit_estimate": "finite-shot circuit estimate",
        "distance_to_ed_annealed": "distance to ED annealed class weights",
        "distance_to_uniform_prior": "distance to uniform class prior",
        "ed_gibbs_reference": "ED Gibbs reference",
        "ed_clean_sector": "ED clean-sector reference",
        "ed_annealed": "ED annealed reference",
        "ed_uniform_quenched": "ED uniform-quenched reference",
        "ed_circuit_weighted": (
            "ED class energies with noiseless-circuit class weights"
        ),
        "one_step_exact_elementary_angle": "one step, exact elementary angle",
        "one_step_linear_elementary_angle": "one step, linear elementary angle",
        "two_step_coherent_reuse": "two steps, coherent reuse",
        "two_step_persistent_fresh": "two steps, persistent fresh records",
        "operational_ti_quenched": "operational-energy TI estimate",
        "ed_uniform_quenched_free_energy": "ED uniform-quenched reference",
        "ed_annealed_free_energy": "ED annealed reference",
        "ed_quenched_ground_limit": "ED quenched ground-state limit",
    }
)


def publication_label(key: str) -> str:
    """Return one controlled manuscript-facing label."""

    try:
        return PUBLICATION_LABELS[key]
    except KeyError as error:
        raise KeyError(f"Unknown publication-label key: {key!r}") from error


__all__ = ["PUBLICATION_LABELS", "publication_label"]
