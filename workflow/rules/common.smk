__author__ = "Hanna Soderstrom"
__copyright__ = "Copyright 2023, Hanna Soderstrom"
__email__ = "hanna.soderstrom@gu.se"
__license__ = "GPL-3"

import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *

min_version("6.8.0")

### Set and validate config file

if not workflow.overwrite_configfiles:
    sys.exit("At least one config file must be passed using --configfile/--configfiles, by command line or a profile!")


validate(config, schema="../schemas/config.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file

units = pandas.read_table(config["units"], dtype=str).set_index(["sample", "type", "flowcell", "lane", "barcode"], drop=False).sort_index()
validate(units, schema="../schemas/units.schema.yaml")

### Set wildcard constraints


wildcard_constraints:
    barcode="[A-Z+]+",
    chr="[^_]+",
    flowcell="[A-Z0-9-]+",
    lane="L[0-9]+",
    sample="|".join(get_samples(samples)),
    type="N|T|R",

if config.get("trimmer_software", "None") == "fastp_pe":
    alignment_input = lambda wilcards: [
        "prealignment/fastp_pe/{sample}_{flowcell}_{lane}_{barcode}_{type}_fastq1.fastq.gz",
        "prealignment/fastp_pe/{sample}_{flowcell}_{lane}_{barcode}_{type}_fastq2.fastq.gz",
    ]
elif config.get("trimmer_software", "None") == "None":
    alignment_input = lambda wildcards: [
        get_fastq_file(units, wildcards, "fastq1"),
        get_fastq_file(units, wildcards, "fastq2"),
    ]


def compile_output_list(wildcards: snakemake.io.Wildcards):
    #print(wildcards.type)
    #print(samples)
    #print(units)
    output_files = [
        #"sentieon/dedup/{}_{}_DEDUP.bam".format(sample, t)
        "sentieon/qualcal/{}_{}_RECAL_DATA.TABLE".format(sample, t)
        for sample in get_samples(samples)
        for t in get_unit_types(units, sample)
    ]
    output_files.append(
        [
            #"sentieon/tnscope/{}.txt".format(sample)
            "sentieon/tnscope/{}_TNscope_tn_ML.vcf".format(sample)
            for sample in get_samples(samples)
        ]
    )
    output_files.append(
        [
            "sentieon/dnascope/{}_{}_DNAscope_modelfiltered.vcf".format(sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )

   # output_files.append(
   #     [
   #         "sentieon/tnscope/%s_TNscope_tn.vcf" % (sample)
   #         for sample in get_samples(samples)
   #     ]
   # )
    return output_files
