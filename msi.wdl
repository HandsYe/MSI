## user: yehuitree@gmail.com
## tags: v0.0.3
workflow msipipeline{
	String result_dir
	String tumor_name
	String? normal_name
	File ref_fasta
	File ref_dic
	File ref_fai
	File tumor_fq1
	File tumor_fq2
	File? normal_fq1
	File? normal_fq2
	File? bed_file

#####################################################
	call fastqc_multqc{
		input:
			result_dir = result_dir,
			tumor_fq1 = tumor_fq1,
			tumor_fq2 = tumor_fq2,
			normal_fq1 = normal_fq1,
			normal_fq2 = normal_fq2
}

	call alignment{
		input:
			result_dir = result_dir,
			ref_fasta = ref_fasta,
			tumor_fq1 = tumor_fq1,
			tumor_fq2 = tumor_fq2,
      normal_fq1 = normal_fq1,
      normal_fq2 = normal_fq2,
			tumor_name = tumor_name,
			normal_name = normal_name			
}

	call MSI{
		input:
			result_dir = result_dir,
			ref_fasta = ref_fasta,
			ref_fai = ref_fai,
			ref_dic = ref_dic,
			tumor_name = tumor_name,
			tumor_bam = alignment.tb,
			tumor_bam_bai = alignment.tbb,
			normal_bam = alignment.nb,
			normal_bam_bai = alignment.nbb,
			normal_name = normal_name,
			bed_file = bed_file
}
}
#
#####################################################
task fastqc_multqc{
	String result_dir
	String QC_dir = "${result_dir}/1_QC"
	File tumor_fq1
	File tumor_fq2
	File? normal_fq1
  File? normal_fq2

	command {
			if [[ ! -f "${QC_dir}/SUCCESS" ]]; then
				rm -rf ${QC_dir}
				mkdir ${QC_dir}
				fastqc ${tumor_fq1} ${tumor_fq2} -o ${QC_dir}

				if [[ ! -z "${normal_fq1}" ]]; then
					fastqc ${normal_fq1} ${normal_fq2} -o ${QC_dir}
				fi
				
				cd ${QC_dir}
				multiqc .
				touch "${QC_dir}/SUCCESS"
			fi
		}
	output {
		File results = stdout()
		}
}

task alignment{
	String result_dir
	String Fastq2Bam_dir = "${result_dir}/2_Fastq2Bam"
	String tumor_name
	String? normal_name
	File ref_fasta
	File ref_fai = "${ref_fasta}.fai"
  File ref_dict = "${ref_fasta}.dict"
  File ref_amb = "${ref_fasta}.amb"
  File ref_ann = "${ref_fasta}.ann"
  File ref_bwt = "${ref_fasta}.bwt"
  File ref_pac = "${ref_fasta}.pac"
  File ref_sa = "${ref_fasta}.sa"
	File tumor_fq1
  File tumor_fq2
  File? normal_fq1
  File? normal_fq2

	command {
			if [[ ! -f "${Fastq2Bam_dir}/SUCCESS" ]]; then	
				rm -rf  ${Fastq2Bam_dir}
				mkdir ${Fastq2Bam_dir}
				bwa mem -t 12 -R "@RG\tID:${tumor_name}\tSM:${tumor_name}\tPL:illumina\tCN:BDLS" ${ref_fasta} ${tumor_fq1} \
				${tumor_fq2} | samtools view -1 -bS - -o ${tumor_name}.bam
				sambamba markdup -t 12 ${tumor_name}.bam ${tumor_name}.markDup.bam 
				sambamba sort -t 12 ${tumor_name}.markDup.bam -o ${tumor_name}.markDup.sort.bam
				cp ${tumor_name}.markDup.sort.bam* ${Fastq2Bam_dir}
		
				if [[ ! -z "${normal_name}" ]]; then
					bwa mem -t 12 -R "@RG\tID:${normal_name}\tSM:${normal_name}\tPL:illumina\tCN:BDLS" ${ref_fasta} ${normal_fq1} \
					${normal_fq2} | samtools view -1 -bS - -o ${normal_name}.bam
					sambamba markdup -t 12 ${normal_name}.bam ${normal_name}.markDup.bam
					sambamba sort -t 12 ${normal_name}.markDup.bam -o ${normal_name}.markDup.sort.bam
					cp ${normal_name}.markDup.sort.bam* ${Fastq2Bam_dir}
				fi
				touch "${Fastq2Bam_dir}/SUCCESS"
				touch "${Fastq2Bam_dir}/.markDup.sort.bam.bai"
				touch "${Fastq2Bam_dir}/.markDup.sort.bam"
			fi
		}
	output {
		File tb = "${Fastq2Bam_dir}/${tumor_name}.markDup.sort.bam"
		File tbb = "${Fastq2Bam_dir}/${tumor_name}.markDup.sort.bam.bai"
		File nb = "${Fastq2Bam_dir}/${normal_name}.markDup.sort.bam"
		File nbb = "${Fastq2Bam_dir}/${normal_name}.markDup.sort.bam.bai"
	}
}

task MSI{
	String result_dir
	String MSI_dir = "${result_dir}/3_MSI"
	String tumor_name
	String? normal_name
	File tumor_bam
	File tumor_bam_bai
	File ref_fasta
	File ref_fai
	File ref_dic
	File normal_bam
	File normal_bam_bai
	File? bed_file

	command {
		if [[ ! -f "${MSI_dir}/SUCCESS" ]];then
			rm -rf ${MSI_dir}
			mkdir ${MSI_dir}
			msisensor scan -d ${ref_fasta} -o b37.list

			if [[! -z "${bed_file}" ]];then
				if [[ ! -z "${normal_name}" ]];then
					msisensor msi -d b37.list -t ${tumor_bam} -n ${normal_bam} -e ${bed_file} -o ${tumor_name}
				else
					msisensor msi -d b37.list -t ${tumor_bam} -e ${bed_file}  -o ${tumor_name}
				fi
			else
				if [[ ! -z "${normal_name}" ]];then
					msisensor msi -d b37.list -t ${tumor_bam} -n ${normal_bam} -o ${tumor_name}
				else
					msisensor msi -d b37.list -t ${tumor_bam} -o ${tumor_name}
				fi
			fi 

			cp ${tumor_name}* ${MSI_dir}
			touch "${MSI_dir}/SUCCESS"
		fi
	}
	output {
		File results = stdout()
	}
}
