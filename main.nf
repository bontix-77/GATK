#!/usr/bin/env nextflow
/*
 * Pipeline parameters
 */

// Primary input

/*
 * Generate BAM index file
 */
process SAMTOOLS_INDEX {

    container 'staphb/samtools:latest'

    publishDir params.outdir, mode: 'symlink'

    input:
    path input_bam

    output:
     tuple path(input_bam), path ("*.bai"), emit: bam_and_index


    script:
    """
    samtools index "${input_bam}"
    """
}
process GATK_varian {

    container 'broadinstitute/gatk:latest'
    publishDir params.outdir, mode: 'symlink'

    input:
    tuple path (file),path (bai)
     tuple path(ref_fasta), path(ref_fai), path(ref_dict)
     val (intervals)
    output:
    path "${file}.g.vcf", emit: GVCF
    path "${file}.g.vcf.idx" , emit: idx

    script:

    """
gatk HaplotypeCaller \
        -R ${ref_fasta}\
        -I ${file}\
        -O "${file}.g.vcf" \
        -L ${intervals}\
        -ERC GVCF
        """
}
process GATK_JOINTGENOTYPING {
       container 'broadinstitute/gatk:latest'
    publishDir params.outdir, mode: 'symlink'

    input:
    path all_gvcfs
    path all_idx
    tuple val (ref), val(fai), val(dict)
    path (intervals)
    val cohort

    output:
       path "${cohort}.joint.vcf"     , emit: vcf
        path "${cohort}.joint.vcf.idx" , emit: idx

     script:
    def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
    """
    gatk GenomicsDBImport \
        ${gvcfs_line} \
        -L ${intervals}\
        --genomicsdb-workspace-path ${cohort}_gdb

    gatk GenotypeGVCFs \
        -R "${ref}" \
        -V gendb://${cohort}_gdb \
        -L "${intervals}"\
        -O "${cohort}.joint.vcf"
    """
}
workflow {

    file_bam = channel.fromPath("${projectDir}/data/sample_bams.txt").splitText()

    // Create index file for input BAM file
    SAMTOOLS_INDEX(file_bam)
    reference = channel.value(tuple (file("${params.ref}/ref.fasta"),file("${params.ref}/ref.fasta.fai"),file("${params.ref}/ref.dict")))
    intervals= channel.value(("${params.ref}/intervals.bed"))
    GATK_varian(SAMTOOLS_INDEX.out.bam_and_index, reference,intervals)


    all_vcf = GATK_varian.out.GVCF.collect()
    all_idx=GATK_varian.out.idx.collect()
    cohort=channel.value("${params.cohort}")
    GATK_JOINTGENOTYPING(all_vcf,all_idx,reference,intervals,cohort)
}
