# PR Checklist

- [ ] Ths PR is up-to-date with the top of all sub-component repositories except for those sub-components which are the subject of this PR. Please consult the ufs-weather-model [wiki](https://github.com/ufs-community/ufs-weather-model/wiki/Making-code-changes-in-the-UFS-weather-model-and-its-subcomponents) if you are unsure how to do this.

- [ ] An Issue describing the work contained in this PR has been created either in the subcomponent(s) or in the ufs-weather-model. The Issue location will depend on which repository is more relevant to the changes in contained in the PR. 

- [ ] If new or updated input data is required by this PR, it is clearly stated in the text of the PR.

Once the PR is opened, ufs-weather-model code managers will use the information provided to add the applicable labels (``No Baseline Change``, ``New Input Data Reqd`` etc) and place it in the Commit Queue. Once the PR is in the Commit Queue, it is the PR owner's responsiblity to keep the PR up-to-date with the develop branch of ufs-weather-model. 


## Instructions: all subsequent sections of text should be filled in as appropriate.

## Description

Provide a detailed description of what this PR does. What bug does it fix, or what feature does it add? Is a change of answers expected from this PR? Are any library updates included in this PR (modulefiles etc.)?

### Issue(s) addressed

Link the issues to be closed with this PR, whether in this repository, or in another repository.
(Remember, issues must always be created before starting work on a PR branch!) 
- fixes #<issue_number>
- fixes noaa-emc/fv3atm/issues/<issue_number>

## Testing

How were these changes tested?
What compilers / HPCs was it tested with?
Are the changes covered by regression tests? (If not, why? Do new tests need to be added?)
Have regression tests and unit tests (utests) been run? On which platforms and with which compilers? (Note that unit tests can only be run on tier-1 platforms)

## Dependencies

If testing this branch requires non-default branches in other repositories, list them.
Those branches should have matching names (ideally)

Do PRs in upstream repositories need to be merged first?
If so add the "waiting for other repos" label and list the upstream PRs
- waiting on noaa-emc/nems/pull/<pr_number>
- waiting on noaa-emc/fv3atm/pull/<pr_number>
