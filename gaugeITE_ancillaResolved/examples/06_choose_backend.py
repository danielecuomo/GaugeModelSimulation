"""The same experiment on local Aer or IBM hardware.

The checked-in example remains on Aer.  IBM credentials are configured once
with ``gauge-ite configure ibm`` and are never pasted into Python code.
"""

from gauge_ite import ExecutionConfig, GaugeITEExperiment, ITEProtocol, TFIMSpec


PROVIDER = "aer"
BACKEND = "auto"
SUBMIT_HARDWARE_JOB = False


spec = TFIMSpec.rectangular(1, 2)
protocol = ITEProtocol(beta=0.5, steps=1, architecture="coherent_reuse")
execution = ExecutionConfig(
    provider=PROVIDER,
    backend=BACKEND,
    profile="default" if PROVIDER == "ibm" else None,
    shots=1024,
    submit=SUBMIT_HARDWARE_JOB,
)
result = GaugeITEExperiment(spec, protocol, execution=execution).run(
    "class_probabilities"
)
print(result.classes.to_string(index=False))
