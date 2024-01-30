# include the config file
configfile: "config.yaml"

# Define the ref based on the config file
# Sort of acts like a global variable so you don't need to always type the whole thing
REF = config['params']['ref_fa']


# define a function to return target files based on config settings
def run_all_input(wildcards):

    run_all_files = []
    # if mapping is set to true add fragment calculation metrics
    if config['modules']['mapping']:

        run_all_files.append("Align/{}_gene_count.txt".format(config['samples']['id']))
        run_all_files.append("Align/{}_transcript_count.txt".format(config['samples']['id']))
        run_all_files.append("Align/{}.uniq.bam_coverage_depth.txt".format(config['samples']['id']))
        run_all_files.append("Align/{}.raw.flagstat.txt".format(config['samples']['id']))
        run_all_files.append("Align/{}.uniq.flagstat.txt".format(config['samples']['id']))

 

    return run_all_files


# rule run all, the files above are the targets for snakemake
rule run_all:
    input:
        run_all_input
smk_path = config['params']['smk_path']
include: smk_path+"/map.smk"
