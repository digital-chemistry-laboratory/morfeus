"""Run in command line mode."""
import fire

import morfeus.buried_volume
import morfeus.cone_angle
import morfeus.conformer
import morfeus.dispersion
import morfeus.local_force
import morfeus.pyramidalization
import morfeus.sasa
import morfeus.sterimol
import morfeus.visible_volume
import morfeus.xtb


def main() -> None:
    """Call Fire to access command line scripts."""
    fire.Fire(
        {
            "buried_volume": morfeus.buried_volume.cli,
            "cone_angle": morfeus.cone_angle.cli,
            "conformer": morfeus.conformer.cli,
            "dispersion": morfeus.dispersion.cli,
            "local_force": morfeus.local_force.cli,
            "pyramidalization": morfeus.pyramidalization.cli,
            "sasa": morfeus.sasa.cli,
            "sterimol": morfeus.sterimol.cli,
            "visible_volume": morfeus.visible_volume.cli,
            "xtb": morfeus.xtb.cli,
        }
    )


if __name__ == "__main__":
    main()
