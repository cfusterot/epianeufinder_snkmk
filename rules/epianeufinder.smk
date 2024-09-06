import glob
import os
from pathlib import Path

rule epianeufinder:
    input:
       path = lambda wc: expand("{dir}/{{sample}}/", dir=units.loc[wc.sample]['path'])
    output:
       "{}/{{sample}}/epianeufinder_results/Karyogram.png".format(OUTDIR)
    params:
       dir=directory("{}/{{sample}}/").format(OUTDIR),
       sample_id="{{sample}}",
       windowSize=config['epianeufinder']['windowSize'],
       blacklist=config['epianeufinder']['blacklist'],
       genome=config['epianeufinder']['genome'],
       reuse=config['epianeufinder']['reuse.existing']
    threads: get_resource("epianeufinder", "threads")
    resources:
        mem_mb=get_resource("epianeufinder", "mem_mb"),
        walltime=get_resource("epianeufinder", "walltime")
    log:
        err="{}/{{sample}}/epianeufinder.err".format(LOGDIR),
        out="{}/{{sample}}/epianeufinder.out".format(LOGDIR)
    script:
        "../scripts/epianeufinder.R"
