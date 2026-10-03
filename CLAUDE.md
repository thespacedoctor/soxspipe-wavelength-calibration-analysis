# CLAUDE.md

This file gives Claude Code guidance for the skyline wavelength calibration anchoring investigation, part of the soxspipe project.

## Project layout

- `input/`: source material for the project. Do not write into this folder.
- `output/`: deliverables and project-adjacent files.
- FITS files (`*.fits`) are too large for version control and are ignored by git.

## Agent skills

### Issue tracker

Issues live in Linear, team `thespacedoctor`, project `soxspipe`. See `docs/agents/issue-tracker.md`.

### Triage labels

Triage states are Linear workflow statuses (`Triage`, `needs-info`, `ready-for-agent`, `ready-for-dave`, `wontfix`); `bug` and `enhancement` are labels. See `docs/agents/triage-labels.md`.

### Domain docs

Single-context: one `CONTEXT.md` and `docs/adr/` at the repo root. See `docs/agents/domain.md`.
