task PullSRACodes {
    input {
        String accession
        Int cpu
        String mem
    }

    command <<< 
        set -eo pipefail
        fasterq-dump ${accession} -e ${cpu}
    >>>

    output {
        Array[File]+ fastq_files = ["${accession}_1.fastq", "${accession}_2.fastq"]
    }

    runtime {
        memory: "${mem}"
        cpu: "${cpu}"
        docker: "cameronnguyen/sratoolkit:latest"
    }
}

task qcFiles {
    input {
        Array[File] fastq_file_pair
        String accession
        Int cpu
        String mem
        Int read_len_req
    }

    command <<<
        fastqc -t ${cpu} -o . ${fastq_file_pair[0]} ${fastq_file_pair[1]}
        fastp -5 -c --overlap_diff_limit 0 -p --trim_poly_x -w ${cpu} --length_required ${read_len_req} -i ${fastq_file_pair[0]} -I ${fastq_file_pair[1]} \
            -o ${accession}_1.trim.fastq -O ${accession}_2.trim.fastq -j ${accession}.fastp.json -h ${accession}.fastp.html
        fastqc -t ${cpu} *trim.fastq
        multiqc .
    >>>

    output {
        Array[File] fastq_html_preqc = ["${accession}_1_fastqc.html","${accession}_2_fastqc.html"]
        Array[File] fastq_html_postqc = ["${accession}_1.trim_fastqc.html","${accession}_2.trim_fastqc.html"]
        Array[File] fastq_postqc = ["${accession}_1.trim.fastq","${accession}_2.trim.fastq"]
        File multiqc_report = "multiqc_report.html"
    }

    runtime {
        memory: "${mem}"
        cpu: "${cpu}"
        docker: "cameronnguyen/qualitycontrol:latest"
    }

}

task callVariants {
    input {
        Array[File] c_fastq_file_pair
        Int cpu
        String mem
        String accession
    }

    command <<<
        bowtie2 --phred33 -q -p ${cpu} -x /usr/local/work/bt2/SARS2 -1 ${c_fastq_file_pair[0]} -2 ${c_fastq_file_pair[1]} -S ${accession}'.output.sam'
        samtools view -S -b ${accession}'.output.sam' > ${accession}'.output.bam'
        samtools sort ${accession}'.output.bam' > ${accession}'.sorted.bam'
        samtools index -o ${accession}'.sorted.bam.bai' ${accession}'.sorted.bam'
        bcftools mpileup -Ou -f /usr/local/work/bt2/NC_045512.2_SARS2.fasta ${accession}'.sorted.bam' | bcftools call -mv -Oz -o ${accession}'.calls.vcf'
        /usr/src/jdk-24.0.1/bin/java -jar /usr/src/snpEff/snpEff.jar -c /usr/local/work/snpEff.config -v SARS2 ${accession}'.calls.vcf' > ${accession}'.ann.vcf'
    >>>

    output {
        Array[File] bam_array = ["${accession}.sorted.bam","${accession}.sorted.bam.bai"]
        File VCF = "${accession}.ann.vcf"
    }
    
    runtime {
        memory: "${mem}"
        cpu: "${cpu}"
        docker: "cameronnguyen/variant_calling:latest"
    }
}

workflow VariantCalling {
    input {
        Array[String] sraCodes
        Int cpu
        String mem
        Int read_len_req
    }

    scatter (sraCode in sraCodes) {
        call PullSRACodes {
            input:
                accession = sraCode,
                cpu = cpu,
                mem = mem
        }
    }

    scatter (index in range(length(sraCodes))) {
        call qcFiles {
            input:
                fastq_file_pair = PullSRACodes.fastq_files[index],
                read_len_req = read_len_req,
                cpu = cpu,
                mem = mem,
                accession = sraCodes[index]
        }
    }

    scatter (index in range(length(sraCodes))) {
        call callVariants {
            input:
                c_fastq_file_pair = qcFiles.fastq_postqc[index],
                cpu = cpu,
                mem = mem,
                accession = sraCodes[index]
        }
    }

    output {
        Array[Array[File]] preqc_summaries = qcFiles.fastq_html_preqc
        Array[Array[File]] postqc_summaries = qcFiles.fastq_html_postqc
        Array[Array[File]] cleaned_fastqs = qcFiles.fastq_postqc
        Array[Array[File]] BAMs_BAIs = callVariants.bam_array
        Array[File] mqc_report = qcFiles.multiqc_report
        Array[File] VCFs = callVariants.VCF
    }
}
