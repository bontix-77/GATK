#!/usr/bin/env nextflow
// ===================================================================================
// PROCESS 1: SAMTOOLS_INDEX
// Purpose: Generates an index (.bai) for a BAM file.
// ===================================================================================
process SAMTOOLS_INDEX {

    container 'staphb/samtools:latest'
    publishDir params.outdir, mode: 'symlink'

    input:
    path input_bam

    output:
    tuple path(input_bam), path("*.bai"), emit: bam_and_index

    script:
    """
    samtools index "${input_bam}"
    """
}

// ===================================================================================
// PROCESS 2: GATK_varian (HaplotypeCaller)
// Purpose: Identifies SNPs and indels in a single sample (GVCF mode).
// ===================================================================================
process GATK_varian {

    container 'broadinstitute/gatk:latest'
    publishDir params.outdir, mode: 'symlink'

    input:
    tuple path(file), path(bai)
    tuple path(ref_fasta), path(ref_fai), path(ref_dict)
    path intervals

    output:
    path "${file}.g.vcf",     emit: GVCF
    path "${file}.g.vcf.idx", emit: idx

    script:
    """
    # -ERC GVCF: Outputs intermediate "Genomic VCF" with genotype likelihoods,
    # which is required for the subsequent Joint Genotyping step.
    gatk HaplotypeCaller \
        -R ${ref_fasta} \
        -I ${file} \
        -O "${file}.g.vcf" \
        -L ${intervals} \
        -ERC GVCF
    """
}

// ===================================================================================
// PROCESS 3: GATK_JOINTGENOTYPING
// Purpose: Combines GVCFs from multiple samples into one finalized VCF.
// ===================================================================================
process GATK_JOINTGENOTYPING {
    
    container 'broadinstitute/gatk:latest'
    publishDir params.outdir, mode: 'symlink'

    input:
    path all_gvcfs
    path all_idx
    tuple val(ref), val(fai), val(dict)
    path intervals
    val cohort

    output:
    path "${cohort}.joint.vcf",     emit: vcf
    path "${cohort}.joint.vcf.idx", emit: idx

    script:
    /* ------------------------------------------------------------------------------
    EDUCATIONAL TAKEAWAY: Groovy Logic inside Scripts
    ------------------------------------------------------------------------------
    Bash cannot easily loop through a Nextflow object (the list of inputs). 
    We use Groovy (the language Nextflow is built on) to manipulate the list of 
    input files into a formatted string that the GATK command line expects.
    
    Code below transforms the list [file1, file2] into "-V file1 -V file2"
    */
    def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
    
    """
    gatk GenomicsDBImport \
        ${gvcfs_line} \
        -L ${intervals} \
        --genomicsdb-workspace-path ${cohort}_gdb

    gatk GenotypeGVCFs \
        -R "${ref}" \
        -V gendb://${cohort}_gdb \
        -L "${intervals}" \
        -O "${cohort}.joint.vcf"
    """
}

// ===================================================================================
// WORKFLOW DEFINITION
// ===================================================================================
workflow {

    // Load input files from text list
    file_bam = channel.fromPath("${projectDir}/prnp.txt").splitText()

    // Step 1: Index the BAMs
    SAMTOOLS_INDEX(file_bam)
    
    // /* ------------------------------------------------------------------------------
    // EDUCATIONAL TAKEAWAY: Value Channels
    // ------------------------------------------------------------------------------
    // Reference genomes are large files used by every process. 
    // If you used a standard queue channel, the reference would be "consumed" 
    // by the first process and unavailable for the second.
    
    // channel.value(...) ensures the file is available indefinitely for all 
    // parallel processes.
    // */
    reference = channel.value(tuple(file("${params.ref}/ref.fasta"), file("${params.ref}/ref.fasta.fai"), file("${params.ref}/ref.dict")))
    intervals = channel.value(file("${params.ref}/intervals.bed"))
    
    // Step 2: Run HaplotypeCaller (Scatter)
    // This runs in parallel for every item in the SAMTOOLS_INDEX output channel.
    GATK_varian(SAMTOOLS_INDEX.out.bam_and_index, reference, intervals)

    // ------------------------------------------------------------------------------
    // EDUCATIONAL TAKEAWAY: The .collect() Operator
    // ------------------------------------------------------------------------------
    // Without .collect(): The next process would trigger every time ONE file 
    // finishes (1 input = 1 output).
    
    // With .collect(): Nextflow blocks the channel until EVERY parallel instance 
    // of GATK_varian is finished. It then emits the data as a single List 
    // (List of all 50 samples, for example). 
    
    // This is essential for the Joint Genotyping step, which needs all files at once.
    // */
    all_vcf = GATK_varian.out.GVCF.collect()
    all_idx = GATK_varian.out.idx.collect()
    
    cohort = channel.value("${params.cohort}")
    
    // Step 3: Joint Genotyping (Gather)
    GATK_JOINTGENOTYPING(all_vcf, all_idx, reference, intervals, cohort)
}
