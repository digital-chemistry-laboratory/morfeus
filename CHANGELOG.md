# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - 2026-01-08

### Fixed
- Compatibility issue with Python 3.14 where command-line interface failed due to `fire` library's requirement for `__name__` attribute on `functools.partial` objects

## [0.8.0] - 2025-11-11

### Added
- Method `XTB.get_dipole_moment` to return the molecular dipole moment (in deybe or a.u.)
- Method `XTB.get_chemical_potential` to calculate the chemical potential
- Method `XTB.get_electronegativity` to calculate the electronegativity
- Method `XTB.get_hardness` to calculate the hardness
- HOMO and LUMO energy can be returned either in Eh of eV depending on the given unit argument
- Synonyms "minus" and "plus" for respectively "nucleophilicity" and "electrophilicity" in the Fukui coefficients varieties
- GFNn-xTB parametrisation (1 or 2) can be given either in `int` or `str` in the `XTB` object
- Method `XTB.get_softness` to calculate the softness
- Possibility to specifiy the number of parallel processes for the xtb runs
- Possibility to set up environment variables for xtb runs
- Support energy units Eh, eV, kcal/mol or kJ/mol
- Methods `XTB.get_atom_dipoles` and `XTB.get_atom_dipole_moments` to return the atomic dipole vectors and dipole moments
- Methods `XTB.get_fod_population` and `XTB.get_nfod` to return the atomic fractional occupation density (FOD) population and the integration over all space of the FOD
- Methods `XTB.get_atom_polarizabilities` and `XTB.get_molecular_polarizability` to return the atomic and molecular polarizabilities
- Method `XTB.get_homo_lumo_gap` to return HOMO-LUMO gap
- Implementation of descriptors calculation with PTB
- Methods `XTB.get_solvation_energy`, `XTB.get_solvation_h_bond_correction` and `XTB.get_atomic_h_bond_corrections` to get respectively the solvation free energy, the hydrogen bonding correction to the solvation energy, and the atomic hydrogen bonding corrections to the solvation free energy
- Possibility to return the CM5 atomic charges instead of Mulliken with GFN1-xTB
- Method `XTB.get_energy` to return the total energy
- Method `XTB.get_fermi_level` to return the Fermi level
- Method `XTB.get_covcn` to return the atomic covalent coordination numbers
- Method `XTB.get_s_pop`, `XTB.get_p_pop` and `XTB.get_d_pop` to return the population partitioned to the s/p/d shells

### Changed
- Require Python>=3.10
- Run xtb for descriptors calculations through the command line instead of the deprecated xtb-python API (xtb-python is still used for the conformers optimmisation via QCEngine)
- Arguments for number of conformers in all morfeus methods renamed into `n_conformers` for consistency. `num_children` also renamed in `n_children` in `conformers_from_ob_ga` method
- Switch to `True` by default for the correction term in the calculations of IP, EA, Fukui local electrophilicity and global descriptors
- `XTB` attribute `version` renamed into `method` to reflect better the choice between GFN2-xTB, GFN1-xTB, and PTB

### Removed
- test `test_homo_index` due to dependency to the deprecated xtb-python API
- Atomic volumes attribute from SASA due to problems in the implementation (see discussion in [issue #76](https://github.com/digital-chemistry-laboratory/morfeus/issues/76))

### Fixed
- Update `pkg_resources.parse_version` to `packaging.version.parse` due to deprecation of pkg_resources

## [0.7.2] - 2022-08-23

### Fixed 
- Calculation of p values of point clouds for Dispersion now works
- Updated PyVista calls to be consistent with changes in their API
- Added `py.typed` for compatibility with other typed code

## [0.7.1] - 2022-06-27

### Fixed
- Fixed type errors related to ConformerEnsemble and Python 3.8 

## [0.7.0] - 2022-06-22

### Added
- Added `SolidAngle` for solid angle calculations.

## [0.6.0] - 2022-03-28

### Added
- Add capacity to Boltzmann average 1D arrays in `ConformerEnsemble`
- `ConeAngle` now uses libconeangle as default with the internal algorithm as backup.
- Added `BiteAngle` for bite angle calculations.

### Fixed
- `XTB.get_fukui`, affecting the varieties "electrophilicity" (wrong sign), "nucleophilicity" (wrong sign), "radical" (wrong number) and "local_electrophilicty" (wrong number)
- Error when `ConformerEnsemble` was initiated with `connectivity_matrix=None`
- Interface to be compatible with version 3 of dftd4.

### Removed
- D3Grimme calculator due to removal from dftd4. The internal D3Calculator remains.

## [0.5.5] - 2021-10-07

### Fixed
- Fixed bug with `ConformerEnsemble.from_crest` with only one structure in the ensemble
- Fixed bug with `ConformerEnsemble.from_ob_ga` when generating RDKit mol 

## [0.5.4] - 2021-07-13

### Fixed
- Fixed bug with `io.write_xyz` and multiple structures

## [0.5.3] - 2021-04-29

### Fixed 
- Typing of external dependencies to allow conda-forge install

## [0.5.2] - 2021-04-28

### Added
- First public release. Changelog use starts here.

### Deprecated
- `BuriedVolume.percent_buried_volume`. Use `BuriedVolume.fraction_buried_volume` instead.`percent_buried_volume` will be reintroduced later with the proper meaning.

