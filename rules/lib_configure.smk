import os

SCRIPT_DIR = config["script_dir"]  # Path to EarlGrey scripts directory
OUTDIR = config["output_dir"]

ARIA_THREADS = 16
DFAM_PARTS = [0, 1]


rule download_dfam:
    output:
        directory("{outdir}/dfam_libraries/")
    cache: True
    conda:
        "envs/download_dfam.yaml"
    threads: max(1, min(16, workflow.cores))  # max 16 threads for aria2
    params:
        dfam_folder_url = "https://dfam.org/releases/current/families/FamDB/",
        dfam_parts_list = {" ".join(DFAM_PARTS)}
    shell:
        """
        for i in "${dfam_parts_list[@]}"
        do
            filename="dfam39_full.${i}.h5.gz"
            url="${dfam_folder_url}/${filename}"
            outfile="{OUTDIR}/dfam_libraries/${filename}"

            echo "Downloading $url"
            aria2c -s {threads} -x {threads} $url -o $outfile

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
        directory("{outdir}/dfam_libraries/")
    output:
        directory("{outdir}/dfam_libraries/")
    cache: True
    conda:
        "envs/repeatmasker.yaml"
    threads: 1
    params:
        dfam_folder_url = "https://dfam.org/releases/current/families/FamDB/",
        dfam_parts_list = {" ".join(DFAM_PARTS)}
    shell:
        """
        DB_PATH=$1

        BIN_DIR=$(dirname $(which RepeatMasker))
        REPEATMASKER_SHARE_DIR=${BIN_DIR}/../share
        LIB_DIR=${REPEATMASKER_SHARE_DIR}/Libraries/
        FAMDB_DIR=${LIB_DIR}/famdb/
        MIN_INIT_PARTITION=${FAMDB_DIR}/min_init.0.h5

        echo "Creating symbolic links from ${FAMDB_DIR} to Dfam DB at ${DB_PATH}"
        for file in ${DB_PATH}/*.h5
        do
            link=${FAMDB_DIR}/$(basename $file)
            ln -s $file $link
        done


        echo "Backing up mini init partition"
        mv $MIN_INIT_PARTITION ${MIN_INIT_PARTITION}.bak

        echo "Configuring RepeatMasker"
        cd $repeat_masker_dir
        perl ./configure \
            -libdir $LIB_DIR \
            -trf_prgm /home/coen/micromamba/envs/earlgrey/bin/trf \
            -rmblast_dir /home/coen/micromamba/envs/earlgrey/bin \
            -hmmer_dir /home/coen/micromamba/envs/earlgrey/bin \
            -abblast_dir /home/coen/micromamba/envs/earlgrey/bin \
            -crossmatch_dir /home/coen/micromamba/envs/earlgrey/bin \
            -default_search_engine rmblast
        """
