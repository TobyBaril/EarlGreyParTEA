import os

SCRIPT_DIR = config["script_dir"]  # Path to EarlGrey scripts directory

DFAM_PARTS = [0, 1]


rule download_dfam:
    """
    Download DFAM libraries for RepeatMasker
    """
    output:
        directory("{outdir}/dfam_libraries/")
    cache: True
    conda:
        "envs/download_dfam.yaml"
    threads: max(1, min(16, workflow.cores))  # max 16 threads for aria2
    params:
        dfam_folder_url = "https://dfam.org/releases/current/families/FamDB/",
        dfam_parts_list = {" ".join([str(i) for i in DFAM_PARTS])}
    shell:
        """
        for i in ${dfam_parts_list}
        do
            filename="dfam39_full.${i}.h5.gz"
            url="${dfam_folder_url}/${filename}"
            outfile="{OUTDIR}/dfam_libraries/${filename}"

            echo "Downloading $url"
            aria2c \
                -s {threads} \
                -x {threads} \
                --optimize-concurrent-downloads \
                --check-integrity=true \
                --max-tries=3 \
                $url \
                -o $outfile

            echo "Checking checksum"
            wget -q ${url}.md5
            md5sum -c ${filename}.md5

            echo "Decompressing $filename"
            pigz -p {threads} -d $filename
        done

        echo "Done"
        """

rule configure_repeatmasker_library:
    input:
        dfam_lib_dir="{outdir}/dfam_libraries/"
    output:
        "{outdir}/.repeatmasker_configuration_done"
    conda:
        "envs/repeatmasker.yaml"
    threads: 1
    shell:
        """
        echo "Configuring RepeatMasker"
        {params.script_dir}/configure_repeatmasker.sh {input.dfam_lib_dir}
        """
