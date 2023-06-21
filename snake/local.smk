rule link_local_data_directories:
    output:
        directory(config["local_data_dirs"]),
    input:
        ancient(
            [f"{config['local_data_root']}/{dir}" for dir in config["local_data_dirs"]]
        ),
    params:
        root=config["local_data_root"],
    shell:
        """
        for dir in {output}
        do
            ln -s "{params.root}/$dir"
        done
        """


rule link_secure_local_data_directories:
    output:
        directory(config["secure_local_data_dirs"]),
    input:
        ancient(
            [
                f"{config['secure_local_data_root']}/{dir}"
                for dir in config["secure_local_data_dirs"]
            ]
        ),
    params:
        root=config["secure_local_data_root"],
    shell:
        """
        for dir in {output}
        do
            ln -s "{params.root}/$dir"
        done
        """


rule link_gtpro_db:
    output:
        directory("ref/gtpro"),
    input:
        ancient("/pollard/data/gt-pro-db"),
    shell:
        alias_recipe_norelative


rule link_gtpro_snps_dict:
    output:
        "ref/gtpro.snp_dict.tsv",
    input:
        ancient("/pollard/data/gt-pro-db/variants_main.covered.hq.snp_dict.tsv"),
    shell:
        alias_recipe_norelative


rule link_midasdb_uhgg:
    output:
        directory("ref/midasdb_uhgg"),
    input:
        ancient("/pollard/scratch/bsmith/midasdb_uhgg/"),
    shell:
        alias_recipe_norelative


rule link_xjin_mgen_data:
    output:
        directory("raw/mgen/xjin"),
    input:
        "/pollard/data/internal/unpublished/CZBMI-Biofilm_BeadExperiment/biofilmBeadExpV2",
    shell:
        alias_recipe_norelative


rule linkage_xjin_genomes_data:
    output:
        directory("raw/genomes/xjin"),
    input:
        "/pollard/data/microbial_genomes/fischbachBiohubStrains/Hybrid_Closed_Genomes/",
    shell:
        alias_recipe_norelative


rule link_watson2022_raw_data:
    output:
        directory("raw/mgen/watson2022"),
    input:
        "/pollard/data/metagenomes/fmt_studies/watson2021_biorxiv/wms_dna",
    shell:
        alias_recipe_norelative


rule link_fmt_studies_raw_data:
    output:
        directory("raw/mgen/fmt_studies"),
    input:
        "/pollard/data/metagenomes/fmt_studies",
    shell:
        alias_recipe_norelative


rule link_eggnog_mapper_db:
    output:
        directory("ref/eggnog_mapper_db"),
    input:
        ancient("/pollard/scratch/vdubinki/eggnog_mapper_data"),
    shell:
        alias_recipe_norelative


localrules:
    link_eggnog_mapper_db,


rule link_dbcan_db:
    output:
        directory("ref/dbcan"),
    input:
        ancient("/pollard/scratch/vdubinki/dbcan/db/"),
    shell:
        alias_recipe_norelative


localrules:
    link_dbcan_db,


rule link_een_mgen_dir:
    output:
        directory("raw/een-mgen/mgen"),
    input:
        "/pollard/data/projects/bsmith/een-mgen/raw/2023-06-16_dropbox_gladstone_een-metagenomes",
    shell:
        alias_recipe_norelative


# rule link_een_gtpro_dir:
#     output:
#         directory("raw/een-mgen/gtpro"),
#     input:
#         "/pollard/data/projects/bsmith/een-mgen/raw/2023-06-01_aritra.mahapatra@tum.de/syncandshare.lrz.de/dl/fiL5YrfYpvTkm19DDBkrTe",
#     shell:
#         alias_recipe_norelative
#
#
# rule link_een_mgen_gtpro_results:
#     output:
#         "data/reads/{mgen}/{r}.proc.gtpro_raw.gz",
#     input:
#         lambda w: "raw/een-mgen/gtpro/" + config["een_mgen_local_src"][w.r][w.mgen],
#     shell:
#         alias_recipe
#
#
# ruleorder: link_een_mgen_gtpro_results > run_gtpro
