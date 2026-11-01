# Repository administration status — 2026-10-26

**Partial rollout, verified against live settings and GitHub's ruleset API.** Nine new active rulesets were saved across five public repositories. Four repositories have both new branch and tag rules plus future-release immutability. The existing Red-team-Challenge rule was preserved. This is not a claim that all 18 repositories are fully locked.

The earlier [content evidence index](PUBLIC_RESEARCH_INDEX.md) now links to this receipt. Its content-audit findings and the dated baseline manifests are preserved. This receipt reports administration and archive integration, not new scientific results.

## Applied controls

Each new branch ruleset targets the default branch, blocks deletion and force pushes, and has an empty bypass list. Ordinary commits remain possible for authorized writers. Each new tag ruleset targets all tags and blocks updates, deletion and force pushes, with no bypass actors. New tag creation remains allowed.

Future-release immutability was enabled and verified again after reloading Settings for Publications, Quantum-Classical-public, quantum-classical-provisional and Artificial-intelligence. It does not retroactively make older releases immutable. Repository owners retain the ability to administer rules; these rules govern normal Git operations under the current configuration.

| Public repository | Default branch rule | Tag rule | Future release immutability | Zenodo auto-archive |
| --- | --- | --- | --- | --- |
| `Artificial-intelligence` | [25894272](https://github.com/Parveen117/Artificial-intelligence/rules/25894272) | [25894317](https://github.com/Parveen117/Artificial-intelligence/rules/25894317) | Enabled; reloaded | Enabled |
| `quantum-classical-provisional` | [25893756](https://github.com/Parveen117/quantum-classical-provisional/rules/25893756) | [25894113](https://github.com/Parveen117/quantum-classical-provisional/rules/25894113) | Enabled; reloaded | Enabled |
| `Cryptography-Recognition-Ledger-Systems` | [25894526](https://github.com/Parveen117/Cryptography-Recognition-Ledger-Systems/rules/25894526) | Pending | Pending verification/application | Enabled |
| `IEL-Information-Entropy-Ledger` | Pending | Pending | Pending verification/application | Enabled |
| `Atomic-validation-public` | Pending | Pending | Pending verification/application | Enabled |
| `EMK-material-response-public` | Pending | Pending | Pending verification/application | Enabled |
| `Quantum-Classical-public` | [25891525](https://github.com/Parveen117/Quantum-Classical-public/rules/25891525) | [25893094](https://github.com/Parveen117/Quantum-Classical-public/rules/25893094) | Enabled; reloaded | Enabled |
| `ECL-pulic` | Pending | Pending | Pending verification/application | Enabled |
| `recognition-seam-calculus` | Pending | Pending | Pending verification/application | Enabled |
| `Publications` | [25890100](https://github.com/Parveen117/Publications/rules/25890100) | [25890492](https://github.com/Parveen117/Publications/rules/25890492) | Enabled; reloaded | Enabled |
| `Recognition-Kernel-Framework` | Pending | Pending | Pending verification/application | Enabled |
| `Thermodynamics-Reproducibility` | Pending | Pending | Pending verification/application | Enabled |
| `Red-team-Challenge` | [20662589](https://github.com/Parveen117/Red-team-Challenge/rules/20662589) (existing) | Pending | Pending verification/application | Enabled |
| `RNKE-GATE` | Pending | Pending | Pending verification/application | Enabled |
| `Bio-tech` | Pending | Pending | Pending verification/application | Enabled |
| `Energy` | Pending | Pending | Pending verification/application | Enabled |
| `Parveen117` | Pending | Pending | Pending verification/application | Off; profile repo |
| `extra-ideas` | Pending | Pending | Pending verification/application | Enabled |

## Zenodo update

- Authenticated GitHub integration was synchronized and reloaded. All 17 research repositories are enabled; the profile repository remains off.
- `quantum-classical-provisional` was disabled despite having an earlier DOI. Automatic archiving is now enabled.
- [Publications v3.3, DOI 10.5281/zenodo.22164767](https://zenodo.org/records/22164767) now has a dated provenance note linking to the evidence index and repository provenance at commit `8101092ded133333a6c96a9fd5883631149cd600`.
- A before/after comparison of the public record API found only the `notes` metadata field changed. The original title, authors, version, publication date, rights metadata and file descriptors/checksum were preserved. No replacement archive or new DOI was created in that earlier administration pass. The subsequent collection deposit is recorded below.

## Published collection snapshot — 26 October 2026

[Version 2026.10.26, DOI 10.5281/zenodo.23671480](https://zenodo.org/records/23671480) archives **17 public research repositories and 3,252 tracked files** at pinned commits. All seven published attachments were downloaded and matched byte for byte; every source file, Git commit object and reconstructed root tree passed verification. See [COLLECTION_ARCHIVE.md](COLLECTION_ARCHIVE.md) for exact scope, commits, checksums and the verification receipt.

This independently versioned collection preserves current source snapshots and original rights/credits. It leaves the earlier per-repository DOI families unchanged. Scientific suites were not rerun for this archival operation; repository administration remains a separate, partially completed rollout.

## Remaining administration

The settings above are the verified result. Repeated browser runtime and input timeouts interrupted the remainder of the rollout. The uncompleted rows are not treated as protected.

Human access was inspected for Publications, Artificial-intelligence and quantum-classical-provisional: each showed zero collaborators and no invited collaborators. Other access audits, app/deploy-key write privileges, and remaining workflow account settings are pending. Public contributors and forks are not evidence of canonical write access.

Complete the pending branch, tag and future-release controls using the verified templates below, preserving existing stronger rules and required checks. Do not require a sole owner to approve their own pull request. Keep the working Zenodo integration and user-authorized tools.

- [Default branch history template](administration/templates/branch-history-ruleset.json)
- [All-tag snapshot template](administration/templates/tag-snapshot-ruleset.json)
- [Machine-readable receipt, live rule IDs and archive checksum](administration/2026-10-26/status.json)
- [Release metadata and rights checklist](REPOSITORY_MAINTENANCE.md)

The next scientific or repository-specific release still requires accurate release-specific citation metadata, resolution of the recorded license discrepancies, and the repository's existing release checks. A provenance note documents traceability; it does not turn a historical archive into a new code snapshot or establish scientific correctness or legal priority.
