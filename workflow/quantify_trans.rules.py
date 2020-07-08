import pandas as pd
configfile: "configs/config_main.yaml"

samples = pd.read_csv(config["METAFILE"], sep = '\t', header = 0)['sample']
trimmed = config["TRIMMED"]
if trimmed:
    input_path = config["OUTPUTPATH"] + "/" + config["PROJECT"] + "/trim"
else:
    input_path = config["READSPATH"]
key = config["KEY"]
end = config["END"]
intermediate_path = config["OUTPUTPATH"] + "/" + config["PROJECT"] + "/trans"
final_path = config["FINALOUTPUT"] + "/" + config["PROJECT"] + "/trans"

rule end:
    input:
        report = final_path + "/report_quantify.html"

if end == "pair":
    rule getReads:
        output:
            forward = temp(intermediate_path + "/reads/{sample}_forward.fastq.gz"),
            reverse = temp(intermediate_path + "/reads/{sample}_reverse.fastq.gz")
        params:
            key = key,
            input_path = input_path
        run:
            shell("scp -i {params.key} {params.input_path}/{wildcards.sample}_*1*.f*q.gz {output.forward}"),
            shell("scp -i {params.key} {params.input_path}/{wildcards.sample}_*2*.f*q.gz {output.reverse}")

else:
    rule getReads:
        output:
            read = temp(intermediate_path + "/reads/{sample}.fastq.gz")
        params:
            key = key,
            input_path = input_path
        run:
            shell("scp -i {params.key} {params.input_path}/{wildcards.sample}.f*q.gz {output.read}")

rule indexTrans:
    input:
        trans = config["TRANS"]
    output:
        index = directory(intermediate_path + "/transcripts_index")
    shell:
        "salmon index -t {input} -i {output} --type quasi -k 31 -p {config[NCORE]}"

if end == "pair":
    rule quantify:
        input:
            forward = temp(intermediate_path + "/reads/{sample}_forward.fastq.gz"),
            reverse = temp(intermediate_path + "/reads/{sample}_reverse.fastq.gz"),
            index = directory(intermediate_path + "/transcripts_index")
        output:
            quant_dir = directory(final_path + "/quant/{sample}"),
            tpm = final_path + "/tpmFile/{sample}_tpm.tsv"
        shell:
            """
            salmon quant -i {input.index} -l A -1 {input.forward} -2 {input.reverse} -o {output.quant_dir} -p {config[NCORE]} --seqBias --useVBOpt --validateMappings 2> /dev/null
            awk 'NR==1{{next}}{{print $1"\\t"$4}}' {output.quant_dir}/quant.sf > {output.tpm}
            """
else:
    rule quantify:
        input:
            read = temp(intermediate_path + "/reads/{sample}.fastq.gz"),
            index = directory(intermediate_path + "/transcripts_index")
        output:
            quant_dir = directory(final_path + "/quant/{sample}"),
            tpm = final_path + "/tpmFile/{sample}_tpm.tsv"
        shell:
            """
            salmon quant -i {input.index} -l A -r {input.read} -o {output.quant_dir} -p {config[NCORE]} --seqBias --useVBOpt --validateMappings 2> /dev/null
            awk 'NR==1{{next}}{{print $1"\\t"$4}}' {output.quant_dir}/quant.sf > {output.tpm}
            """

rule summaryReport:
    input:
        quant_dir = directory(expand(final_path + "/quant/{sample}", sample = samples))
    output:
        report = final_path + "/report_quantify.html"
    shell:
        "multiqc {input} --filename {output} --quiet"
        
