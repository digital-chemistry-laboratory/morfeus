# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Add capacity to Boltzmann average 1D arrays in `ConformerEnsemble`
- `ConeAngle` now uses libconeangle as default with the internal algorithm as backup.
- Added `BiteAngle` for bite angle calculations.

### Fixed
- Fixed ``XTB.get_fukui``, affecting the varieties "electrophilicity" (wrong sign), "nucleophilicity" (wrong sign), "radical" (wrong number) and "local_electrophilicty" (wrong number)
- Fixed crash when `ConformerEnsemble` was initiated with `connectivity_matrix=None`

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

