import os
import tempfile
import subprocess
import os.path
import json
import argparse


#parser = argparse.ArgumentParser(description="Location of boost include folder")
#parser.add_argument("--boost",help="Location of boost install/include folder. \
g#	Most likely {anaconda install}/envs/\{proj_env\}/include")

# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)


boost=config['boost']


rule all:
	input:
		samtobed="samToBed",
		vit="viterbi",

rule stb:
	input:
		sTb="SamToBed.cpp",
	output:
		samtobed="samToBed",
	shell:"""

		g++ -O2 {input.sTb} -o {output.samtobed}
"""


rule vit:
	input:
		vitt="viterbi.cpp",
	output:
		samtobed="viterbi",
	params:
		proj_env={boost},

	shell:"""
g++ -O2 -W -I {params.proj_env} {input.vitt} -o {output.samtobed}
"""


	
