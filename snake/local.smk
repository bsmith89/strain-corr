rule link_local_data_directories:
    output:
        directory(config["local_data_dirs"]),
    input:
        [f"{config['local_data_root']}/{dir}" for dir in config["local_data_dirs"]],
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
        [
            f"{config['secure_local_data_root']}/{dir}"
            for dir in config["secure_local_data_dirs"]
        ],
    params:
        root=config["secure_local_data_root"],
    shell:
        """
        for dir in {output}
        do
            ln -s "{params.root}/$dir"
        done
        """
