rule count_ref_geno_positions:
    output:
        "data/species/sp-{species}/midasdb_{dbv}.gtpro.geno.npositions.tsv",
    input:
        script="scripts/count_ref_geno_positions.py",
        geno="data/species/sp-{species}/midasdb_{dbv}.gtpro.mgtp.nc",
    params:
        ambig_thresh=0.01,
    conda:
        "conda/sfacts.yaml"
    shell:
        "{input.script} {input.geno} {params.ambig_thresh} {output}"


# FIXME: This is generic to either v15 or v20,
# # NOTE (2024-06-10): This is now no longer true. The script name should also probably be updated.
# # but points at v15 only for
# # npositions data because it should be the same as v20 and I don't want to re-run it.
rule collect_metadata_for_uhgg_ref_strains_new:
    output:
        meta="data/species/sp-{species}/midasdb_{dbv}.gene{centroid}.strain_meta.tsv",
    input:
        script="scripts/extract_metadata_midasdb_v15.py",
        meta="ref/midasdb_uhgg_{dbv}/metadata/genomes-all_metadata.tsv",
        genome_to_species="ref/midasdb_uhgg_{dbv}/genomes.tsv",
        pos="data/species/sp-{species}/midasdb_v20.gtpro.geno.npositions.tsv",
        genes="data/species/sp-{species}/midasdb.gene{centroid}_{dbv}.uhgg-strain_gene.tsv",
    shell:
        "{input.script} {input.meta} {input.genome_to_species} {input.pos} {input.genes} {wildcards.species} {output}"


rule collect_filtering_metadata_for_uhgg_ref_strains_new:
    output:
        "data/species/sp-{species}/midasdb_{dbv}.gene{centroid}.strain_meta-complete90-contam5-pos{pos}.tsv",
    input:
        script="scripts/filter_ref_strains_v15.py",
        meta="data/species/sp-{species}/midasdb_{dbv}.gene{centroid}.strain_meta.tsv",
    params:
        min_completeness=90 / 100,
        max_contamination=5 / 100,
        min_positions=lambda w: int(w.pos),
    shell:
        "{input.script} {input.meta} {params.min_completeness} {params.max_contamination} {params.min_positions} {output}"


# NOTE: Split from `compile_spgc_to_ref_strain_report_new`:
rule calculate_gene_prevalence_in_ref_genomes:
    output:
        "{stem}/midasdb.gene{centroid}_{dbv}.{unit}-strain_gene.prevalence.tsv",
    input:
        script="scripts/strain_gene_to_prevalence.py",
        gene="{stem}/midasdb.gene{centroid}_{dbv}.{unit}-strain_gene.tsv",
        filt="{stem}/midasdb_{dbv}.gene{centroid}.strain_meta-complete90-contam5-pos0.tsv",
    params:
        pseudo=0,
    shell:
        "{input.script} {input.gene} {input.filt} {params.pseudo} {output}"

# NOTE: This rule takes the super long filename and turns it into a much shorter one for, e.g., notebooks.
rule alias_sfacts_outputs:
    output:
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.sfacts-fit.{stemB}",
    input:
        source=lambda w: (
            "data/group/{w.group}/species/sp-{w.species}/{w.stemA}.gtpro.{sfacts_params}.{w.stemB}".format(
                w=w,
                sfacts_params=get_sfacts_stem(config, w.species, w.group),
            )
        ),
    shell:
        alias_recipe


localrules:
    alias_sfacts_outputs,


# NOTE: This ruleorder section is a place to clear up ambiguity about whether the pipeline should be run on the full sfacts spec (yes)
# or on previously aliased files (no).
# When I get errors like this:
#
# AmbiguousRuleException:
# Rules export_sfacts_comm and alias_sfacts_outputs are ambiguous for the file data/group/hmp2/species/sp-102506/r.proc.gtpro.sfacts-fit.comm.tsv.
# Consider starting rule output with a unique prefix, constrain your wildcards, or use the ruleorder directive.
#
# I can just add that rule to the end of the list.
ruleorder: alias_sfacts_outputs > export_sfacts_comm

# That forces the aliasing to the end of the computation.
