rule ccs:
    """
    Create HiFi reads from SequelII IsoSeq subreads
    """
    input:
        "../../input/subreads/{pool}.subreads.bam"
    output:
        "output/ccs/{pool}/{pool}.hifi.bam"
    log:
        "log/ccs/{pool}/{pool}.ccs.log"
    params:
        outdir="output/ccs/{pool}",
        prefix="output/ccs/{pool}/{pool}"
    threads: 64
    conda:
        "isoseq"
    benchmark:
        "stats/ccs_lima/ccs_{pool}_benchmark.txt"
    shell:
        "mkdir -p {params.outdir} && ccs --report-file {params.prefix}.report.txt --metrics-json {params.prefix}.metrics.json.gz --hifi-summary-json {params.prefix}.hifi-summary.json --num-threads {threads} --log-level INFO --log-file {log} {input} {output}"


# lima
for idx, pool in enumerate(config['pools']):
    n2bc = config['name2bc'][idx]
    bcs = list(n2bc.values())
    rule lima:
        """
        Demultiplex the pools and remove primers + barcodes (TO DO: do the
        barcodes present also contain primers? --> Barcodes + universal primers
        from the downloaded isoseq file)
        """
        input:
            hifi=f"output/ccs/{pool}/{pool}.hifi.bam",
            biosamples=f"../../input/meta/biosamp_{pool}.csv",
            barcodes="../../input/meta/pb16_barcodes_uniprimers.fasta"
        output: 
            expand(f"output/lima/{pool}/{pool}.demux.hifi.{{bc}}.bam", bc=bcs )
        log:
            f"log/lima/{pool}.lima.log"
        params:
            pattern=f"output/lima/{pool}/{pool}.demux.hifi.bam"
        threads: 64
        conda:
            "isoseq"
        benchmark:
            f"stats/ccs_lima/lima_{pool}_benchmark.txt"
        shell:
            "lima"
            " --isoseq"
            " --peek-guess"
            " --min-length 100"
            " --num-threads {threads}"
            " --log-level INFO"
            " --log-file {log}"
            " --store-unbarcoded"
            " --split-bam"
            " --biosample-csv {input.biosamples}"
            " {input.hifi}"
            " {input.barcodes}"
            " {params.pattern}"

    rule renameit:
        input:
            lambda wc: f"output/lima/{pool}/{pool}.demux.hifi.{n2bc[wc.sample]}.bam"
        output:
            f"output/lima/{pool}/{pool}.demux.hifi.{{sample}}.renamed.bam"
        shell:
            "ln -sf $(readlink -f {input}) {output}"

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
