"""TFIM circuit construction with explicit multi-step architectures."""

from __future__ import annotations

from dataclasses import dataclass
from math import log2

from qiskit import ClassicalRegister, QuantumCircuit

from ..models import TFIMSpec
from ..protocols import ITEProtocol
from .primitive import append_term_interaction


@dataclass(frozen=True)
class RegisterLayout:
    system_qubits: tuple[int, ...]
    record_qubits_by_step: tuple[tuple[int, ...], ...]
    purification_qubits: tuple[int, ...]
    num_record_bits: int
    dynamic_record_register_name: str | None = None

    @property
    def record_qubits(self) -> tuple[int, ...]:
        return tuple(qubit for row in self.record_qubits_by_step for qubit in row)


@dataclass(frozen=True)
class CircuitBundle:
    circuit: QuantumCircuit
    spec: TFIMSpec
    protocol: ITEProtocol
    layout: RegisterLayout

    @property
    def interpretation(self) -> str:
        return self.protocol.interpretation


def _prepare_system_input(
    circuit: QuantumCircuit,
    spec: TFIMSpec,
    protocol: ITEProtocol,
    purification_qubits: tuple[int, ...],
) -> None:
    system = tuple(range(spec.num_system_qubits))
    if protocol.input_mode == "purification":
        # The legacy package pairs system qubit i with the reference qubit at
        # ``total_qubits - 1 - i``. Preserve that order exactly.
        for system_qubit, reference_qubit in zip(
            system, reversed(purification_qubits)
        ):
            circuit.h(system_qubit)
            circuit.cx(system_qubit, reference_qubit)
    elif protocol.input_mode == "basis":
        if len(protocol.basis_bits or ()) != spec.num_system_qubits:
            raise ValueError(
                f"basis_bits must contain {spec.num_system_qubits} entries."
            )
        for qubit, bit in enumerate(protocol.basis_bits or ()):
            if bit:
                circuit.x(qubit)
    circuit.barrier()


def build_tfim_circuit(spec: TFIMSpec, protocol: ITEProtocol) -> CircuitBundle:
    """Build the circuit while preserving the legacy interaction primitive.

    ``coherent_reuse`` exactly reproduces the original ancilla schedule.
    ``fresh_trajectory`` defers measurement but allocates fresh ancillas for
    every step. ``measure_reset_trajectory`` records and resets the same
    ancillas after every step.
    """

    n_system = spec.num_system_qubits
    n_terms = spec.num_terms
    n_reference = n_system if protocol.input_mode == "purification" else 0
    if protocol.architecture == "fresh_trajectory":
        n_ancilla = protocol.steps * n_terms
    else:
        n_ancilla = n_terms
    total_qubits = n_system + n_ancilla + n_reference

    if protocol.architecture == "measure_reset_trajectory":
        record_register = ClassicalRegister(protocol.steps * n_terms, "record")
        circuit = QuantumCircuit(total_qubits)
        circuit.add_register(record_register)
        dynamic_name = record_register.name
    else:
        circuit = QuantumCircuit(total_qubits)
        record_register = None
        dynamic_name = None

    system_qubits = tuple(range(n_system))
    ancilla_start = n_system
    purification_start = n_system + n_ancilla
    purification_qubits = tuple(
        range(purification_start, purification_start + n_reference)
    )
    _prepare_system_input(circuit, spec, protocol, purification_qubits)

    if protocol.architecture == "fresh_trajectory":
        record_rows = tuple(
            tuple(
                range(
                    ancilla_start + step * n_terms,
                    ancilla_start + (step + 1) * n_terms,
                )
            )
            for step in range(protocol.steps)
        )
    else:
        shared = tuple(range(ancilla_start, ancilla_start + n_terms))
        record_rows = tuple(shared for _ in range(protocol.steps))

    if protocol.architecture == "coherent_reuse":
        for ancilla in record_rows[0]:
            circuit.h(ancilla)
        circuit.barrier()
        for _step in range(protocol.steps):
            for term, ancilla in zip(spec.terms, record_rows[0]):
                append_term_interaction(
                    circuit,
                    term,
                    ancilla,
                    beta=protocol.beta,
                    steps=protocol.steps,
                    mapping=protocol.angle_mapping,
                )
            circuit.barrier()
        for ancilla in record_rows[0]:
            circuit.sx(ancilla)
        output_rows = (record_rows[0],)
        num_record_bits = n_terms
    else:
        for step, ancillas in enumerate(record_rows):
            for ancilla in ancillas:
                circuit.h(ancilla)
            circuit.barrier()
            for term, ancilla in zip(spec.terms, ancillas):
                append_term_interaction(
                    circuit,
                    term,
                    ancilla,
                    beta=protocol.beta,
                    steps=protocol.steps,
                    mapping=protocol.angle_mapping,
                )
            circuit.barrier()
            for local_index, ancilla in enumerate(ancillas):
                circuit.sx(ancilla)
                if record_register is not None:
                    cbit = step * n_terms + local_index
                    circuit.measure(ancilla, record_register[cbit])
                    if step + 1 < protocol.steps:
                        circuit.reset(ancilla)
            circuit.barrier()
        output_rows = record_rows
        num_record_bits = protocol.steps * n_terms

    layout = RegisterLayout(
        system_qubits=system_qubits,
        record_qubits_by_step=output_rows,
        purification_qubits=purification_qubits,
        num_record_bits=num_record_bits,
        dynamic_record_register_name=dynamic_name,
    )
    circuit.metadata = {
        "gauge_ite": True,
        "architecture": protocol.architecture,
        "interpretation": protocol.interpretation,
        "angle_mapping": protocol.angle_mapping,
        "steps": protocol.steps,
        "beta": protocol.beta,
        "input_mode": protocol.input_mode,
        "system_qubits": list(system_qubits),
        "record_qubits_by_step": [list(row) for row in output_rows],
        "purification_qubits": list(purification_qubits),
    }
    return CircuitBundle(circuit=circuit, spec=spec, protocol=protocol, layout=layout)


def add_final_measurements(
    bundle: CircuitBundle,
    basis: str = "ancilla",
) -> QuantumCircuit:
    """Return a copy with final record and optional system measurements."""

    basis = basis.lower()
    if basis not in {"ancilla", "z", "x"}:
        raise ValueError("basis must be 'ancilla', 'z', or 'x'.")
    if bundle.protocol.architecture == "measure_reset_trajectory":
        circuit = bundle.circuit.copy()
        include_system = basis in {"z", "x"}
        if basis == "x":
            for qubit in bundle.layout.system_qubits:
                circuit.h(qubit)
        if include_system:
            system_measurement = ClassicalRegister(
                len(bundle.layout.system_qubits), "system"
            )
            circuit.add_register(system_measurement)
            for index, qubit in enumerate(bundle.layout.system_qubits):
                circuit.measure(qubit, system_measurement[index])
        circuit.metadata = {
            **(circuit.metadata or {}),
            "measurement_basis": basis,
            "num_record_bits": bundle.layout.num_record_bits,
            "num_system_bits": (
                len(bundle.layout.system_qubits) if include_system else 0
            ),
        }
        return circuit

    circuit = bundle.circuit.copy()
    include_system = basis in {"z", "x"}
    if basis == "x":
        for qubit in bundle.layout.system_qubits:
            circuit.h(qubit)

    total_bits = bundle.layout.num_record_bits
    if include_system:
        total_bits += len(bundle.layout.system_qubits)
    measurement = ClassicalRegister(total_bits, "measurement")
    circuit.add_register(measurement)
    for index, qubit in enumerate(bundle.layout.record_qubits):
        circuit.measure(qubit, measurement[index])
    if include_system:
        offset = bundle.layout.num_record_bits
        for index, qubit in enumerate(bundle.layout.system_qubits):
            circuit.measure(qubit, measurement[offset + index])
    circuit.metadata = {
        **(circuit.metadata or {}),
        "measurement_basis": basis,
        "num_record_bits": bundle.layout.num_record_bits,
        "num_system_bits": (
            len(bundle.layout.system_qubits) if include_system else 0
        ),
    }
    return circuit


def circuit_resource_row(spec: TFIMSpec, protocol: ITEProtocol) -> dict:
    bundle = build_tfim_circuit(spec, protocol)
    dense_gib = 16 * (2**bundle.circuit.num_qubits) / 2**30
    return {
        "lattice": f"{spec.lattice.rows}x{spec.lattice.cols}",
        "boundary": spec.lattice.boundary,
        "system_qubits": spec.num_system_qubits,
        "bond_terms": spec.num_bond_terms,
        "field_terms": spec.num_system_qubits,
        "cycle_rank": spec.lattice.cycle_rank,
        "gauge_classes": spec.lattice.num_classes,
        "architecture": protocol.architecture,
        "steps": protocol.steps,
        "input_mode": protocol.input_mode,
        "record_bits": bundle.layout.num_record_bits,
        "total_circuit_qubits": bundle.circuit.num_qubits,
        "dense_statevector_GiB": dense_gib,
        "dense_statevector_log2_amplitudes": int(log2(2**bundle.circuit.num_qubits)),
    }
