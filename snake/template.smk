# Master rule, starts up a brand-new project from a repo clone.
rule new_project:
    output:
        touch("build/new_project.flag"),
    input:
        [
            ancient(flagfile)
            for flagfile in [
                "build/initialize_project_config.flag",
                "build/initialize_git_from_template.flag",
                "build/initialize_files_from_template.flag",
            ]
        ],


# Re-arrange some files that were used for template dev.
rule initialize_files_from_template:
    output:
        touch("build/initialize_files_from_template.flag"),
    shell:
        "ln -fs doc/notes.md README.md"


# Configure git so it's not pointing at the template.
rule initialize_git_from_template:
    output:
        touch("build/initialize_git_from_template.flag"),
    input:
        ancient("build/initialize_files_from_template.flag"),
    shell:
        dd(
            """
        git remote rename origin template-origin
        git checkout -b main
        git add -A
        git commit -m "Initial project commit."
        git commit --amend
        """
        )


rule initialize_project_config:
    output:
        touch("build/initialize_project_config.flag"),
    input:
        local_config_files=[ancient(p) for p in ["env_local", "snake/local.smk"]],
    shell:
        dd(
            """
        git config --local filter.dropoutput_ipynb.clean scripts/ipynb_output_filter.py
        git config --local filter.dropoutput_ipynb.smudge cat
        git config --local diff.daff-csv.command "daff.py diff --git"
        git config --local merge.daff-csv.name "daff.py tabular merge"
        git config --local merge.daff-csv.driver "daff.py merge --output %A %O %A %B"
        git config --local alias.l "log --graph --pretty=format:'%Cred%h%Creset -%C(yellow)%d%Creset %s %Cgreen(%cr)' --abbrev-commit"
        git config --local alias.wdiff "diff --color-words"
        git config --local alias.s "status -uno"
        git config --local alias.ls-other "ls-files --other --exclude-standard"
        git config --local alias.grep-other "!git ls-other | git grep -f -"
        git config --local alias.grep-hist '!f() {{ patt=$1; shift; git grep $patt $(git rev-list --all) $@; }}; f'
        git config --local pull.ff only
        echo 'Add local configuration to {input.local_config_files}'
        echo 'Or by creating/relinking profile/default/'
        echo 'Remember to symlink data directories to the correct fs. (e.g. raw/, ref/, data/, etc.)'
        """
        )


rule build_empty_local_snakefile:
    output:
        touch("snake/local.smk"),


rule build_empty_local_env:
    output:
        touch("env_local"),
