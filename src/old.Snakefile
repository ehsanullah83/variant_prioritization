import os

#grab SAMPLE name from vcf input
SAMPLE = str(config['input_vcf']).split('.vcf')[0]
# logic to parse ped for families
# can either be  multiple families (a list) in the yaml as below
# need to grab the family names from it
if type(config['ped']) == list:
	PEDfile = [x.split('.ped')[0].split('/')[-1] for x in config['ped']]
	PEDbase = '/'.join(config['ped'][0].split('.ped')[0].split('/')[:-1]) + '/'
# otherwise just one family can be provided, which would be a str as below:
else:
	PEDfile = config['ped'].split('.ped')[0].split('/')[-1]
	PEDbase = '/'.join(config['ped'].split('.ped')[0].split('/')[:-1]) + '/'

if PEDbase == '/':
	PEDbase = ''

if type(config['family_name']) == str:
	family_name_con = config['family_name']
else:
	family_name_con = '|'.join(config['family_name'])

def pick():
	# modifies VEP call to either pick most canonical tx or
	# return all tx possibilities
	if config['pick'].upper() == 'CANONICAL':
		out = '--pick_allele_gene  --pick_order canonical, tsl, biotype, ccds, length '
	if config['pick'].upper() == 'ALL':
		out = '--flag_pick_allele_gene '
	return(out)

# import regions
REGIONS_file = config['regions']
if '/home/$USER' in REGIONS_file:
	REGIONS_file = os.environ['HOME'] + REGIONS_file.split('$USER')[-1]
REGIONS = open(REGIONS_file).readlines()
REGIONS = [r.strip() for r in REGIONS]
MT_CONTIGS="MT,GL000207.1,GL000226.1,GL000229.1,GL000231.1,GL000210.1,GL000239.1,GL000235.1,GL000201.1,GL000247.1,GL000245.1,GL000197.1,GL000203.1,GL000246.1,GL000249.1,GL000196.1,GL000248.1,GL000244.1,GL000238.1,GL000202.1,GL000234.1,GL000232.1,GL000206.1,GL000240.1,GL000236.1,GL000241.1,GL000243.1,GL000242.1,GL000230.1,GL000237.1,GL000233.1,GL000204.1,GL000198.1,GL000208.1,GL000191.1,GL000227.1,GL000228.1,GL000214.1,GL000221.1,GL000209.1,GL000218.1,GL000220.1,GL000213.1,GL000211.1,GL000199.1,GL000217.1,GL000216.1,GL000215.1,GL000205.1,GL000219.1,GL000224.1,GL000223.1,GL000195.1,GL000212.1,GL000222.1,GL000200.1,GL000193.1,GL000194.1,GL000225.1,GL000192.1,NC_007605"

# set global Snakemake wildcard constraints
wildcard_constraints:
	sample=SAMPLE,
	region = '|'.join(REGIONS),
#	region = '^[0-9a-zA-Z]+:\d+-\d+'
	family_name=family_name_con

localrules: query_gemini

if config['family_name'] == '':
	rule all:
		input:
			expand('{sample}.PED_{ped}.gemini.db', sample=SAMPLE, ped=PEDfile)
else:
	rule all:
		input:
			expand('sample_reports/{sample}.{family_name}.PED_{ped}.lenient{gemini_lenient}.SeeGEM.report.html', \
				sample=SAMPLE, \
				ped=PEDfile, \
				family_name=config['family_name'], \
				gemini_lenient = config['gemini_lenient']),

rule n_split_vcf:
	input:
		vcf = config['input_vcf']
	output:
		vcf = temp('temp/{sample}__{region}.vcf.gz'),
		index = temp('temp/{sample}__{region}.vcf.gz.tbi')
	shell:
		"""
		export REF_CACHE=/scratch/$SLURM_JOB_ID/
		module load {config[samtools_version]}
		if [[ {wildcards.region} != "MT_contigs" ]]; then
			bcftools view -r {wildcards.region} {input.vcf} | bgzip > {output.vcf}
		else
			bcftools view -r {MT_CONTIGS} {input.vcf} | bgzip > {output.vcf}
		fi
		tabix -p vcf {output.vcf}
		"""

rule vt_bgzip_and_tabix_vcf:
	input:
		'temp/{sample}__{region}.vcf.gz'
	output:
		vcf = temp('temp/vt.{sample}__{region}.vcf.gz'),
		index = temp('temp/vt.{sample}__{region}.vcf.gz.tbi')
	shell:
		"""
		export REF_CACHE=/scratch/$SLURM_JOB_ID/
		module load {config[samtools_version]}
		module load {config[vt_version]}
		zcat {input} \
			| sed 's/ID=AD,Number=./ID=AD,Number=R/' \
			| vt decompose -s - \
			| vt normalize -r {config[ref_genome]} - \
			| bgzip -c > {output.vcf}
		tabix -f -p vcf {output.vcf}
		"""

rule annovar:
	input:
		vcf = 'temp/vt.{sample}__{region}.vcf.gz',
		index = 'temp/vt.{sample}__{region}.vcf.gz.tbi'
	output:
		avinput = temp('temp/{sample}__{region}.avinput'),
	threads: 2
	shell:
		"""
		module unload python
		module load {config[annovar_version]}
		module load {config[InterVar_version]}
		module load {config[R_version]}
		convert2annovar.pl -format vcf4old {input.vcf} -includeinfo --outfile {output.avinput}
		table_annovar.pl {output.avinput} \
			$ANNOVAR_DATA/hg19 \
			-buildver hg19 \
			-remove \
			-out {output.avinput}.INTERVAR \
			--protocol refGene,esp6500siv2_all,1000g2015aug_all,avsnp147,dbnsfp33a,clinvar_20170905,exac03,dbscsnv11,dbnsfp31a_interpro,rmsk,ensGene,knownGene,refGeneWithVer,popfreq_max_20150413,gnomad_exome,gnomad_genome,spidex,avsnp150 \
			-operation  g,f,f,f,f,f,f,f,f,r,g,g,g,f,f,f,f,f \
			--argument '-hgvs',,,,,,,,,,,,'-splicing 50 -hgvs',,,,, \
			--polish -nastring . \
		    --thread {threads} \
			--otherinfo
		"""

rule intervar:
	input:
		avinput = 'temp/{sample}__{region}.avinput',
	output:
		temp('temp/{sample}__{region}.annovar_intervar.out')
	shell:
		"""
		module unload python
		module load {config[InterVar_version]}
		module load {config[R_version]}
		InterVar \
			-i {input.avinput} \
			--input_type=AVinput \
			-d $ANNOVAR_DATA/hg19 \
			-o {input.avinput}.INTERVAR \
			--skip_annovar

		cut -f -139 {input.avinput}.INTERVAR.hg19_multianno.txt | \
			sed "1 s/"Otherinfo"/"CHROM\\\tPOS\\\tID\\\tREF\\\tALT\\\tQUAL\\\tFILTER"/" - > \
			{input.avinput}.INTERVAR.hg19_multianno.txt.CUT

		Rscript {config[intervar_Rscript_path]} \
			{input.avinput}.INTERVAR.hg19_multianno.txt.intervar \
			{input.avinput}.INTERVAR.hg19_multianno.txt.CUT \
			{output}
		"""

rule merge_annovar_intervar:
	input:
		expand('temp/{sample}__{region}.annovar_intervar.out', sample=SAMPLE, region=REGIONS)
	output:
		temp('temp/annovar_intervar_annotation.txt.gz')
	shell:
		"""
		module load {config[samtools_version]}
		head -n 1 {input} | grep "^CHROM" | uniq > temp/annovar_intervar_header
		cat {input} | grep -v "^CHROM" | sort -k1,1 -k2,2n > {output}TEMP
		cat temp/annovar_intervar_header {output}TEMP | bgzip -f > {output}
		tabix -b 2 -e 2 -S 1 {output}
		"""


# annotate with VEP
# two paths here, set by config['pick'] in config.yaml and the pick() function:
# 1. 'canonical' only returns one consequence per variant.
# 		Use this for clinical returns, as no risk of using odd tx with
# 		high consequence as choice in gemini
# 2. 'all' will do multiple annotations with VEP and gemini will use
# 	the most serious consequence.
# 		Use for more research returns, to increase chances of finding
# 		interesting variants at cost of lower specificity.
rule VEP_annotate:
	input:
		'temp/vt.{sample}__{region}.vcf.gz'
	output:
		vcf = temp('temp/{sample}__{region}.SORTED.VT.VEP.vcf.gz'),
		index = temp('temp/{sample}__{region}.SORTED.VT.VEP.vcf.gz.tbi')
	params:
		pick = pick()
	shell:
		"""
		module load {config[VEP_version]}
		module load {config[samtools_version]}
		vep -i {input} --offline \
			--cache --dir_cache $VEPCACHEDIR \
			--fasta $VEPCACHEDIR/GRCh37.fa --species human --assembly GRCh37  \
			--format vcf \
			--output_file {output.vcf} \
			--plugin Grantham \
			--plugin dbscSNV,$VEPCACHEDIR/dbscSNV1.1.txt.gz \
            --plugin GeneSplicer,$GS/bin/genesplicer,$GS/human,context=200 \
            --plugin SpliceRegion \
			--plugin MaxEntScan,/data/OGVFB/resources/MaxEntScan \
			--plugin CADD,/fdb/CADD/1.3/prescored/whole_genome_SNVs.tsv.gz,/fdb/CADD/1.3/prescored/InDels.tsv.gz \
			--canonical \
			--ccds \
			--total_length \
			--hgvs \
			--sift b \
			--polyphen b \
			--symbol \
			--numbers \
			--biotype \
			--total_length \
			--pubmed \
			--domains \
			--gene_phenotype \
			{params.pick} \
            --fields Allele,Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,CANONICAL,DOMAINS,CLIN_SIG,PICK,Grantham,MaxEntScan,HGVSc,HGVSp,PUBMED,Phenotypes,CADD_RAW,CADD_PHRED,GeneSplicer,SpliceRegion,ada_score,rf_score \
			--vcf --compress_output bgzip --force_overwrite
		# tabix
		tabix -f -p vcf {output.vcf}
		"""

# annotate with vcfanno
rule vcfanno_annotate:
	input:
		annovar_intervar = 'temp/annovar_intervar_annotation.txt.gz',
		vcf = 'temp/{sample}__{region}.SORTED.VT.VEP.vcf.gz',
		index = 'temp/{sample}__{region}.SORTED.VT.VEP.vcf.gz'
	output:
		new_conf = temp('temp/{sample}__{region}_vcfanno.conf'),
		vcf = temp('temp/{sample}__{region}.SORTED.VT.VEP.VCFANNO.vcf.gz'),
		index = temp('temp/{sample}__{region}.SORTED.VT.VEP.VCFANNO.vcf.gz.tbi')
	threads: 4
	shell:
		"""
		module load {config[samtools_version]}
		module load {config[vcfanno_version]}
		# copy conf to local dir, then edit to put in the path
		# to the annovar_intervar this Snakemake creates
		cp {config[vcfanno_conf]} {output.new_conf}
		annovar_intervar_path=`echo {input.annovar_intervar} | sed 's:/:\\\\\/:g'`
		sed -i "s/REPLACE_ME_PROGRAMATICALLY/$annovar_intervar_path/g" {output.new_conf}
		vcfanno -p {threads} -lua {config[vcfanno_lua]} {output.new_conf} {input.vcf} | bgzip > {output.vcf}
		tabix -f -p vcf {output.vcf}
		"""

# fix number=A issue
# Since I decompose variants into multiple lines,
# I can make number=1 so vcf2b will keep
rule tweak_header:
	input:
		'temp/{sample}__{region}.SORTED.VT.VEP.VCFANNO.vcf.gz'
	output:
		vcf = 'temp/{sample}__{region}.SORTED.VT.VEP.VCFANNO.NUMFIX.vcf.gz',
		index = 'temp/{sample}__{region}.SORTED.VT.VEP.VCFANNO.NUMFIX.vcf.gz.tbi'
	shell:
		"""
		module load {config[samtools_version]}
		zcat {input} | sed 's/Number=A/Number=1/g' | bgzip -c > {output.vcf}
		tabix -f -p vcf {output.vcf}
		"""


# merge vcfs into one again
rule merge_vcf:
	input:
		vcf = expand('temp/{{sample}}__{region}.SORTED.VT.VEP.VCFANNO.NUMFIX.vcf.gz', region=REGIONS),
		index = expand('temp/{{sample}}__{region}.SORTED.VT.VEP.VCFANNO.NUMFIX.vcf.gz.tbi', region=REGIONS)
	output:
		vcf = ('{sample}.SORTED.VT.VEP.VCFANNO.vcf.gz'),
		index = ('{sample}.SORTED.VT.VEP.VCFANNO.vcf.gz.tbi')
	threads: 16
	shell:
		"""
		export REF_CACHE=/scratch/$SLURM_JOB_ID/
		module load {config[samtools_version]}
		bcftools concat --threads {threads} {input.vcf} | bcftools sort -T /scratch/$SLURM_JOB_ID/ -m 100G -O z -o {output.vcf}
		tabix -f -p vcf {output.vcf}
		"""

#if config['re_sort'] == 'TRUE':
	# ensure that the concat didn't mess up the order
#	rule sort_tabix:
#		input:
#			vcf = 'temp/{sample}.SORTED.VT.VEP.VCFANNO.vcf.gz'
#		output:
#			vcf = protected('{sample}.RESORTED.VT.VEP.VCFANNO.vcf.gz'),
#			index = '{sample}.RESORTED.VT.VEP.VCFANNO.vcf.gz.tbi'
#		shell:
#			"""
#			export REF_CACHE=/scratch/$SLURM_JOB_ID/
#			/home/mcgaugheyd/bin/gsort_linux_amd64 --memory 60000 {input.vcf} /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/GRCh37_gatk_order.genome | bgzip -c > {output.vcf}
#			tabix -f -p vcf {output.vcf}
#			"""
#else:
#	rule faux_sort_tabix:
#		input:
#			vcf = 'temp/{sample}.SORTED.VT.VEP.VCFANNO.vcf.gz',
#			index = 'temp/{sample}.SORTED.VT.VEP.VCFANNO.vcf.gz.tbi'
#		output:
#			vcf = protected('{sample}.RESORTED.VT.VEP.VCFANNO.vcf.gz'),
#			index = '{sample}.RESORTED.VT.VEP.VCFANNO.vcf.gz.tbi'
#		shell:
#			"""
#			module load {config[samtools_version]}
#			mv {input.vcf} {output.vcf}
#			mv {input.index} {output.vcf}
#			"""

rule peddy_QC:
	input:
		'{sample}.SORTED.VT.VEP.VCFANNO.vcf.gz'
	output:
		ped = '{sample}_PEDDY.ped_check.csv',
		het = '{sample}_PEDDY.het_check.csv',
		sex = '{sample}_PEDDY.sex_check.csv'
	threads: 4
	shell:
		"""
		module load {config[peddy_version]}
		peddy -p {threads} {input} {config[ped]} --prefix {wildcards.sample}_PEDDY
		"""


# create gemini database
rule make_gemini_db:
	input:
		vcf = '{sample}.SORTED.VT.VEP.VCFANNO.vcf.gz',
#		index = 'temp/{sample}.SORTED.VT.VEP.VCFANNO.vcf.gz.tbi'
	output:
		'{sample}.PED_{ped}.gemini.db'
	shell:
		"""
		module load {config[vcf2db_version]}
		echo {wildcards.ped}.ped
		vcf2db.py {input.vcf} {PEDbase}{wildcards.ped}.ped {output}
		"""

# now write report for each family given in the yaml family_name section
rule query_gemini:
	input:
		db = '{sample}.PED_{ped}.gemini.db',
		peddy_ped = '{sample}_PEDDY.ped_check.csv',
		peddy_het = '{sample}_PEDDY.het_check.csv',
		peddy_sex = '{sample}_PEDDY.sex_check.csv'
	output:
		report_name = '{sample}.{family_name}.PED_{ped}.lenient{gemini_lenient}.SeeGEM.report.html',
#		report_path = 'sample_reports/{sample}.{family_name}.PED_{ped}.lenient{gemini_lenient}.SeeGEM.report.html'
	resources: res=1
	run:
		report_name = output.report_name
	#	report_path = output.report_path
		if config["output_raw"].upper() == 'NO':
			shell("mkdir -p sample_reports; \
					module load {config[R_version]}; \
					module load {config[gemini_version]}; \
					module load {config[pandoc_version]}; \
					Rscript {config[SeeGEM_script]} \
						{input.db} \
						{wildcards.family_name} \
						{output.report_name} \
						{wildcards.sample}_PEDDY \
						{config[aaf_change]} \
						{config[gemini_lenient]}")
#			shell("mv " + report_name + ' ' + report_path)
		else:
			raw_name = report_name.replace('.html','.tsv')
			#raw_path = report_path.replace('.html','.tsv')
			shell("mkdir -p sample_reports; \
					module load {config[R_version]}; \
					module load {config[gemini_version]}; \
					module load {config[pandoc_version]}; \
					Rscript {config[SeeGEM_script]} \
						{input.db} \
						{wildcards.family_name} \
						{output.report_name} \
						{wildcards.sample}_PEDDY \
						{config[aaf_change]} \
						{config[gemini_lenient]} " + raw_name)
#			shell("mv " + report_name + ' ' + report_path)
#			shell("mv " + raw_name + ' ' + raw_path)

localrules: move_reports
rule move_reports:
	input:
		'{sample}.{family_name}.PED_{ped}.lenient{gemini_lenient}.SeeGEM.report.html'
	output:
		'sample_reports/{sample}.{family_name}.PED_{ped}.lenient{gemini_lenient}.SeeGEM.report.html'
	shell:
		"""
		mv {input} {output}
		"""
