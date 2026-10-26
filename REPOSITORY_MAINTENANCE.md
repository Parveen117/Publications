# Repository maintenance and archival checklist

The 2026-10-26 content audit covers the 18 public repositories in [the evidence index](PUBLIC_RESEARCH_INDEX.md). Private repositories and unpublished material are outside this index.

## Repository administration

The applied portion of the 2026-10-26 administration pass is recorded in [ADMINISTRATION_STATUS.md](ADMINISTRATION_STATUS.md), including live rule IDs, verified future-release settings, the Zenodo integration update, and pending work.

The list below describes target settings. Consult the administration status receipt for the subset actually applied and verified; this list is not a blanket completion claim. Verify them in GitHub Settings using the owner account:

1. **Collaborators:** retain owner control; remove unneeded human write collaborators and pending write invitations. Do not confuse commit authors or public contributors with current collaborators.
2. **Apps and deploy keys:** inspect which integrations can write. Keep only explicitly authorized integrations; preserve the working Zenodo integration. Check workflow token permissions and whether Actions can create/approve pull requests.
3. **Branches:** use active rulesets to prevent force pushes and deletion. Restrict canonical updates to the owner or an explicitly authorized integration. Preserve existing required checks. Configure a workable owner-controlled merge path before enforcing review rules; a sole owner cannot approve their own pull request.
4. **Tags:** protect published release and freeze tags against updates and deletion. Preserve historical tags and their exact targets.
5. **Releases:** enable release immutability in Settings → General → Releases. This applies to future releases; it does not retroactively lock old releases. Prepare all intended assets before publishing an immutable release.
6. **Verify:** record the applied rule IDs, enforcement state, permitted bypass actors and an access audit privately. `CODEOWNERS` alone is not evidence that these settings are enabled.

Public repositories remain readable and forkable. These controls protect the canonical repository; they do not prevent copies. Reuse rights remain those of the applicable existing licenses and notices.

## Before the next Zenodo release

- Link each new scientific claim to its exact files, assumptions, reproduction command and current review status. Archive only work that has passed that repository's existing release gates.
- Preserve existing archived versions. Use a new release for new content or corrections; never move an old release tag to the new commit.
- Review citation titles, version numbers, author/contributor credit and rights metadata against the actual release. Several existing citation files describe earlier packages.
- `Publications/.zenodo.json` currently describes the historical native-seam paper and fixes a July 2026 publication date. It overrides `CITATION.cff` for Zenodo. Prepare metadata appropriate to the actual next release rather than letting that historical metadata label a new portfolio snapshot.
- Reconcile existing rights discrepancies before publishing another archive. The observed Artificial-intelligence, Atomic-validation-public and Recognition-Kernel-Framework Zenodo records display CC-BY-4.0 while current repository notices reserve rights. This audit does not resolve the legal scope or change either existing notice.
- Confirm the repository is enabled in Zenodo, then publish a release at the verified commit. A DOI badge or a historical archive does not prove that the integration is currently enabled.
- After Zenodo processes it, verify the version DOI, concept DOI, record creation time, file checksum, archive content and any Software Heritage archival status. Do not report an archive or timestamp receipt as complete while it is pending.
- Add the verification receipt in a later commit. Do not edit a dated snapshot to make it describe newer work.

## Source guidance

- [GitHub rulesets](https://docs.github.com/en/repositories/configuring-branches-and-merges-in-your-repository/managing-rulesets/available-rules-for-rulesets)
- [GitHub release immutability](https://docs.github.com/en/code-security/how-tos/secure-your-supply-chain/establish-provenance-and-integrity/prevent-release-changes)
- [Zenodo metadata precedence](https://help.zenodo.org/docs/github/describe-software/zenodo-json/)
- [Archive a GitHub release](https://help.zenodo.org/docs/github/archive-software/github-upload/)
- [Zenodo version records](https://help.zenodo.org/docs/deposit/manage-versions/)
