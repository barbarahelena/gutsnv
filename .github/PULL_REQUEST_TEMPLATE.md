<!--
# gutsnv pull request

Please fill in the checklist below (delete whatever is not relevant).
PRs should be made against the dev branch, unless you're preparing a release.
-->

## PR checklist

- [ ] This comment contains a description of changes (with reason).
- [ ] If you've fixed a bug or added code that should be tested, add tests!
- [ ] Make sure your code lints (`nextflow lint .`).
- [ ] Ensure the test suite passes (`nf-test test tests/default.nf.test`).
- [ ] If parameters changed: `nextflow_schema.json` is regenerated (`python .github/scripts/params_to_schema.py main.nf nextflow_schema.json`).
- [ ] Usage documentation in `docs/usage.md` is updated.
- [ ] Output documentation in `docs/output.md` is updated.
- [ ] `CHANGELOG.md` is updated.
- [ ] `README.md` is updated (including new tool citations and authors/contributors).
