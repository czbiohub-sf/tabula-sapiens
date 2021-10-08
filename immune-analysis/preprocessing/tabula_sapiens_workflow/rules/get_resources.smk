rule get_imgt_db:
    # https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/advanced/references#segment
    output:
        fasta="{}/db/10X/{}/fasta/regions.fa".format(
            workflow.basedir, config["vdj_ref_prefix"]
        ),
        ref=directory(
            "{}/db/10X/{}/".format(workflow.basedir, config["vdj_ref_prefix"])
        ),
        outdir=directory("{}/db/10X/".format(workflow.basedir)),
    params:
        species="Homo sapiens",
        cell_ranger=config["cell_ranger"],
        genome=config["vdj_ref_prefix"],
        scripts="{}/scripts".format(workflow.basedir),
    conda:
        "../envs/cellranger.yaml"
    shell:
        "mkdir -p {output.outdir} "
        "&& cd {output.outdir} "
        " && rm -rf {output.ref}"
        "&& bash {params.scripts}/cell_ranger_imgt_pull.sh {params.cell_ranger} "
        "{params.genome} {params.species}"


rule get_immcantation_image:
    """ pull the immcantation image using singularity

    """
    output:
        "{}/resources/repertoire/{}.sif".format(workflow.basedir, "immcantation"),
    threads: 1
    params:
        name="container_pull",
        docker_address="docker://immcantation/suite:4.1.0",
    resources:
        mem_mb=10000,
    shell:
        "module load system && "
        "mkdir -p $(dirname {output}) && cd $(dirname {output}) && "
        "img=$(basename {output}) && "
        "singularity build $img "
        "{params.docker_address}"
