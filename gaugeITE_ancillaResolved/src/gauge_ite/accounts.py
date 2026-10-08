"""Secure, notebook-free IBM account configuration helpers."""

from __future__ import annotations

from getpass import getpass
from typing import Any


def _runtime_service_class():
    try:
        from qiskit_ibm_runtime import QiskitRuntimeService
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise ImportError(
            "Install `gauge-ite[ibm]` before configuring an IBM account."
        ) from exc
    return QiskitRuntimeService


def _redact(value: Any) -> Any:
    if isinstance(value, dict):
        return {
            key: ("<stored securely>" if "token" in str(key).lower() else _redact(item))
            for key, item in value.items()
        }
    if isinstance(value, (list, tuple)):
        return [_redact(item) for item in value]
    return value


def configure_ibm_account(
    *,
    token: str | None = None,
    profile: str = "default",
    instance: str | None = None,
    channel: str = "ibm_quantum_platform",
    overwrite: bool = False,
    set_as_default: bool = True,
) -> dict[str, object]:
    """Save an IBM Runtime account locally and return only safe metadata.

    If ``token`` is omitted, it is requested with a hidden terminal prompt.
    API tokens should never be typed into a notebook or committed script.
    """

    secret = token if token is not None else getpass("IBM Quantum API token: ")
    if not secret.strip():
        raise ValueError("An IBM Quantum API token is required.")
    service_class = _runtime_service_class()
    kwargs: dict[str, object] = {
        "token": secret.strip(),
        "channel": channel,
        "name": profile,
        "overwrite": overwrite,
        "set_as_default": set_as_default,
    }
    if instance:
        kwargs["instance"] = instance
    service_class.save_account(**kwargs)
    return {
        "profile": profile,
        "channel": channel,
        "instance": instance,
        "token": "<stored securely>",
        "set_as_default": set_as_default,
    }


def ibm_account_status(profile: str | None = None) -> dict[str, object]:
    """Report whether a saved IBM profile exists without revealing secrets."""

    service_class = _runtime_service_class()
    accounts = service_class.saved_accounts(name=profile) if profile else service_class.saved_accounts()
    if profile:
        found = bool(accounts)
        selected = accounts
    else:
        found = bool(accounts)
        selected = accounts
    return {
        "profile": profile,
        "configured": found,
        "saved_accounts": _redact(selected),
    }


__all__ = ["configure_ibm_account", "ibm_account_status"]
