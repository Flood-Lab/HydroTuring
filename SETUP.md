# One-time setup

Run these from the repository root. They need the GitHub CLI
(`brew install gh`) and one `gh auth login`.

## 1. Create the organisation

Organisations cannot be created from the CLI. Do this once in a browser:

    https://github.com/organizations/plan  ->  choose Free  ->  name: HydroTuring

## 2. Create the repository and push

    cd ~/Documents/Python/HydroTuring
    gh repo create HydroTuring/HydroTuring \
        --public \
        --source=. \
        --remote=origin \
        --description "Does your AI hydrologic model conserve mass, energy and momentum?" \
        --homepage "https://hydroturing.github.io/HydroTuring/" \
        --push

## 3. Turn on Pages

    gh api -X POST repos/HydroTuring/HydroTuring/pages \
        -f 'build_type=workflow'

The site goes live at https://hydroturing.github.io/HydroTuring/ once the
`pages` workflow finishes. Check it with:

    gh run watch

## 4. Topics, so people can find it

    gh repo edit HydroTuring/HydroTuring --add-topic \
        hydrology,benchmark,machine-learning,physics-informed,water-balance,earth-science

## 5. Labels the contribution flow expects

    gh label create accepted        -c 2C7A57 -d "Proposal approved, the probe is yours" -R HydroTuring/HydroTuring
    gh label create probe-proposal  -c 0E6C78 -d "A proposed new probe"                  -R HydroTuring/HydroTuring
    gh label create good-first-probe -c 7BC4A4 -d "A good place to start"                -R HydroTuring/HydroTuring
    gh label create help-wanted     -c A8501C -d "Unclaimed, we want this"               -R HydroTuring/HydroTuring

## 6. Check it worked

    gh workflow list -R HydroTuring/HydroTuring   # probe, model, pages
    gh run list      -R HydroTuring/HydroTuring

The `probe` workflow runs the acceptance gate. If it is green, the benchmark
demonstrably discriminates on GitHub's runners and not only on your laptop.

## Afterwards

Ask me to file the fourteen roadmap probes as labelled issues. A markdown
list gets read once; issues get claimed.
