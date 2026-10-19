for bc, sample in zip(rename_in, rename_out):
    rule:
        input:
            bc
        output:
            sample
        shell:
            "mv -v {input} {output}"

#rule rename:
#    input:
#        "output/lima/pool1.demux.hifi.{bc}.bam"
#    output:
#        "output/lima/pool1.demux.hifi.{config[bc2sample][0][bc]}.bam"
#    shell:
#        "mv {input} {output}"
