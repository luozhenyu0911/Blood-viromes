# include the config file
configfile: "config.yaml"

# define a function to return target files based on config settings
def run_all_input(wildcards):

    run_all_files = []
    # if mapping is set to true add fragment calculation metrics
    if config['modules']['mapping']:

        run_all_files.append("03_align2pg/{}.rm_pg.R1.fq.gz".format(config['samples']['id']))
        run_all_files.append("03_align2pg/{}.rm_pg.R2.fq.gz".format(config['samples']['id']))
        run_all_files.append('04_annotation/{}.krona.html'.format(config['samples']['id']))
        run_all_files.append('04_annotation/{}_viral_unclass.R1.fa'.format(config['samples']['id']))
        run_all_files.append('04_annotation/{}_viral_unclass.R2.fa'.format(config['samples']['id']))
        run_all_files.append("04_annotation/{}.k2_report.txt".format(config['samples']['id']))
        run_all_files.append("metrics/{}_seq_count.txt".format(config['samples']['id']))
        run_all_files.append("04_annotation/{}.beta_diversity.stats".format(config['samples']['id']))
        run_all_files.append("04_annotation/{}.alpha_beta_diversity.stats".format(config['samples']['id']))
        run_all_files.append("04_annotation/{}_virus_unclss.blastn.res.txt".format(config['samples']['id']))
        run_all_files.append("04_annotation/{}_virus_unclss.blastn.res.anno.txt".format(config['samples']['id']))
        
    return run_all_files


# rule run all, the files above are the targets for snakemake
rule run_all:
    input:
        run_all_input
		
smk_path = config['params']['smk_path']
include: smk_path+"/map.smk"
include: smk_path+"/annotation.smk"
include: smk_path+"/metrics.smk"
