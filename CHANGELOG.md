# phac-nml/aborator: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2025-07-28

### Added
- Numerous pytest workflow tests. [PR #18](https://github.com/phac-nml/arborator/pull/18)
- Warnings for currently unused parameters, incorrect parameters. Exceptions for parameter errors. [PR #18](https://github.com/phac-nml/arborator/pull/18)
- Enforcing that thresholds must be strictly decreasing. [PR #18](https://github.com/phac-nml/arborator/pull/18)
- Generally improved program robustness. [PR #18](https://github.com/phac-nml/arborator/pull/18)

### Fixed
- GitHub Workflow CI repository name. [PR #18](https://github.com/phac-nml/arborator/pull/18)
- Corrected and consolidated variable names in README, config.json examples. [PR #18](https://github.com/phac-nml/arborator/pull/18)
- Conflicts when running multiple runs in a row in the same Python execution. [PR #18](https://github.com/phac-nml/arborator/pull/18)
- A bug with hiding metadata. [PR #18](https://github.com/phac-nml/arborator/pull/18)
- Consolidated true/false boolean/string handling. [PR #18](https://github.com/phac-nml/arborator/pull/18)
- A bug where min_members was always 2. [PR #18](https://github.com/phac-nml/arborator/pull/18)

## [1.0.6] - 2025-05-23

### Added
- GitHub Workflows CI. [PR #15](https://github.com/phac-nml/arborator/pull/15)
- Very simple pytest-workflows testing. [PR #15](https://github.com/phac-nml/arborator/pull/15)

### Fixed
- Reorganized test data and updated the README to reflect this. [PR #15](https://github.com/phac-nml/arborator/pull/15)
- Updated many versions in `setup.py`. [PR #15](https://github.com/phac-nml/arborator/pull/15)

### Removed
- An example run's "results" directory that was included in the source. [PR #15](https://github.com/phac-nml/arborator/pull/15)
- An empty `tests.py` file in the root directory. [PR #15](https://github.com/phac-nml/arborator/pull/15)

[1.0.6]: https://github.com/phac-nml/arborator/releases/tag/1.0.6
[1.1.0]: https://github.com/phac-nml/arborator/releases/tag/1.1.0
