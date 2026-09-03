# One-time setup

Run these from the repository root. They need the GitHub CLI
(`brew install gh`) and one `gh auth login`.

## 1. Create the repository and push

The repository lives in the existing `Flood-Lab` organisation, so there is no
organisation to create first.


    cd ~/Documents/Python/HydroTuring
    gh repo create Flood-Lab/HydroTuring \
        --public \
        --source=. \
        --remote=origin \
        --description "Does your AI hydrologic model conserve mass, energy and momentum?" \
        --homepage "https://flood-lab.github.io/HydroTuring/" \
        --push

## 2. Turn on Pages

    gh api -X POST repos/Flood-Lab/HydroTuring/pages \
        -f 'build_type=workflow'

The site goes live at https://flood-lab.github.io/HydroTuring/ once the
`pages` workflow finishes. Check it with:

    gh run watch

## 3. Topics, so people can find it

    gh repo edit Flood-Lab/HydroTuring --add-topic \
        hydrology,benchmark,machine-learning,physics-informed,water-balance,earth-science

## 4. Labels the contribution flow expects

    gh label create accepted        -c 2C7A57 -d "Proposal approved, the probe is yours" -R Flood-Lab/HydroTuring
    gh label create probe-proposal  -c 0E6C78 -d "A proposed new probe"                  -R Flood-Lab/HydroTuring
    gh label create good-first-probe -c 7BC4A4 -d "A good place to start"                -R Flood-Lab/HydroTuring
    gh label create help-wanted     -c A8501C -d "Unclaimed, we want this"               -R Flood-Lab/HydroTuring

## 5. Check it worked

    gh workflow list -R Flood-Lab/HydroTuring   # probe, model, pages
    gh run list      -R Flood-Lab/HydroTuring

The `probe` workflow runs the acceptance gate. If it is green, the benchmark
demonstrably discriminates on GitHub's runners and not only on your laptop.

## Afterwards

Ask me to file the fourteen roadmap probes as labelled issues. A markdown
list gets read once; issues get claimed.
