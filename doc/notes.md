# sMWAS: Strain-level correlations with host phenotypes

TODO:

- [ ] Implement a deconvolution evaluation tool that measures
    - how well the entropy of the metagenotype and the community are associated
    - whether high-coverage strains have lower entropy
    - whether comparisons between subjects look more different than within
- [ ] TODO
- [x] Export de-identified subject metadata from ucfmt2 database.
    - see ucfmt2 project repository: lib.project_data.load_subject_table(con=con) where con=sqlite3.connect('sdata/database.db')
