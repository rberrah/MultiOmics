# Security policy

MultiOmics is research software. It is not a medical device, clinical decision-support system, or repository for patient data.

## Supported source lines

- `1.1.x`: current maintained research-software line.
- `1.0.0 / Research v1`: frozen reproducibility baseline; retained for provenance rather than ongoing feature development.

## Reporting a vulnerability

Do not include credentials, patient information, controlled data, or exploit details in a public issue.

For a security vulnerability, prefer GitHub's private vulnerability-reporting / Security Advisory channel when it is available for this repository. If that channel is unavailable, contact the repository maintainer through the GitHub profile before publishing technical exploit details.

Ordinary reproducibility, scientific-method, documentation, or non-sensitive bug reports may use public issues.

## Data boundary

Public CI, issues and pull requests must use public or synthetic data only. Confidential clinical, PK, genomic or other controlled data must remain in authorized storage and must not be uploaded to this repository merely to reproduce a software issue.

## Scope

A security fix does not imply scientific, statistical, biological, PK/PD or clinical validation. Scientific limitations remain tracked separately in the release documentation.
