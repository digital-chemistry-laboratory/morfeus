"""Run in command line mode."""

from functools import wraps
from typing import Any, Callable

import fire

import morfeus.buried_volume
import morfeus.cone_angle
import morfeus.conformer
import morfeus.dispersion
import morfeus.local_force
import morfeus.multiwfn
import morfeus.pyramidalization
import morfeus.sasa
import morfeus.solid_angle
import morfeus.sterimol
import morfeus.visible_volume
import morfeus.xtb


def _wrap_cli(func: Callable[..., Any]) -> Callable[..., Any]:
    """Wrap CLI functions to ensure that they have a ``__name__`` attribute.

    Args:
        func: CLI function to wrap

    Returns:
        A wrapped version of ``func`` that sets ``__name__`` on the returned
        functools.partial object, as required by Fire in Python 3.14+.
    """

    @wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        result = func(*args, **kwargs)
        if not hasattr(result, "__name__"):
            result.__name__ = func.__name__
        return result

    return wrapper


def main() -> None:
    """Call Fire to access command line scripts."""
    fire.Fire(
        {
            "bite_angle": _wrap_cli(morfeus.bite_angle.cli),
            "buried_volume": _wrap_cli(morfeus.buried_volume.cli),
            "cone_angle": _wrap_cli(morfeus.cone_angle.cli),
            "conformer": _wrap_cli(morfeus.conformer.cli),
            "dispersion": _wrap_cli(morfeus.dispersion.cli),
            "local_force": _wrap_cli(morfeus.local_force.cli),
            "pyramidalization": _wrap_cli(morfeus.pyramidalization.cli),
            "multiwfn": _wrap_cli(morfeus.multiwfn.cli),
            "sasa": _wrap_cli(morfeus.sasa.cli),
            "solid_angle": _wrap_cli(morfeus.solid_angle.cli),
            "sterimol": _wrap_cli(morfeus.sterimol.cli),
            "visible_volume": _wrap_cli(morfeus.visible_volume.cli),
            "xtb": _wrap_cli(morfeus.xtb.cli),
        }
    )


if __name__ == "__main__":
    main()
