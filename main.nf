#!/usr/bin/env nextflow

process CONCAT_FASTQS {
    tag "${sample_id}"
    publishDir "${launchDir}/larry-results-${params.project_tag}/concat_fastqs/", mode: 'copy'
    
    input:
    tuple val(sample_id), path(r1_files), path(r2_files)
    
    output:
    tuple val(sample_id), path("${sample_id}_R1.fastq.gz"), path("${sample_id}_R2.fastq.gz"), emit: concatenated
    
    when:
    r1_files.size() > 1 || r2_files.size() > 1
    
    script:
    def r1_list = r1_files instanceof List ? r1_files.join(' ') : r1_files
    def r2_list = r2_files instanceof List ? r2_files.join(' ') : r2_files
    """
    # Concatenate R1 files
    cat ${r1_list} > ${sample_id}_R1.fastq.gz
    
    # Concatenate R2 files  
    cat ${r2_list} > ${sample_id}_R2.fastq.gz
    """
}

process FIND_LARRY_SEQS {

  publishDir "${launchDir}/larry-results-${params.project_tag}/before_qc/", mode: 'copy'

  input:
  tuple val(larry_samp), val(gex_samp), val(group), path(r1_file), path(r2_file)

  output:
  tuple val(larry_samp), val(gex_samp), val(group), path("*.pkl")

  script:
  """
  export PIGZ='-p 8' 
  python ${baseDir}/bin/find_larry_seqs.py ${params.fastqs_path} ${r1_file} ${r2_file} ${larry_samp}
  """
}

process LARRY_QC {

  publishDir "${launchDir}/larry-results-${params.project_tag}/after_qc/", mode: 'copy'

  input:
  tuple val(larry_samp), val(gex_samp), val(group), path(pkl)

  output:
  tuple val(larry_samp), val(gex_samp), val(group), path("*.pkl")
  path("*.pdf"), optional: true

  memory { pkl.size() < 300.MB ? ( 97.GB * ( pkl.size() / ( 1024 * 1024 * 1024 ) ) * task.attempt ) : 160.GB * task.attempt }

  script:
  """
  python ${baseDir}/bin/larry_qc.py ${pkl} ${larry_samp} ${params.check_whitelist ? '--check_whitelist' : ''} ${params.make_pdf ? '--make_pdf' : ''} --whitelist_csv ${baseDir}/data/larry_whitelist.csv
  """
}

process ASSIGN_CLONES {

  publishDir "${launchDir}/larry-results-${params.project_tag}/clones/", mode: 'copy', saveAs: {filename -> "${group}_${filename}"}

  input:
  tuple val(larry_samp), val(gex_samp), val(group), path(pkl)

  output:
  tuple val(larry_samp), val(gex_samp), val(group), path("*.pkl")

  script:
  """
  python ${baseDir}/bin/assign_clones.py ${pkl.join(",")} ${params.dispr_filter}
  """
}

process MATCH_GEX {

  publishDir "${launchDir}/larry-results-${params.project_tag}/clones/", mode: 'copy'

  input:
  tuple val(larry_samp), val(gex_samp), val(group), path(pkl)

  output:
  tuple val(larry_samp), path("*.h5ad"), path("*.csv")
  path("*.png"), optional: true

  memory = { 8.GB * Math.max(1, ((larry_samp.size() / 2) as int)) * task.attempt }

  script:
  """
  python ${baseDir}/bin/match_gex.py ${larry_samp.join(',')} ${params.sample_csv} ${params.ss_out} ${group} ${pkl} ${params.gex_path} ${params.plot_cumulative ? '--plot_cumulative' : ''}
  """
}

def extractSampleInfo(file_path) {
    def matcher = file_path.name =~ /^(.+)_S(\d+)_L(\d{3})_R([12])_/
    if (!matcher) {
        error "Unexpected filename format: ${file_path.name}"
    }
    def sample_id = matcher[0][1]
    def sample_num = matcher[0][2] as Integer
    def lane = matcher[0][3] as Integer
    def mate = matcher[0][4] as Integer
    
    return [sample_id, sample_num, lane, mate, file_path]
}


workflow CONCAT_SAMPLE_FASTQS {
    take:
    fastq_files // channel of fastq file paths
    
    main:
    // Extract sample information and group by sample ID
    grouped_samples = fastq_files
        .map { file -> extractSampleInfo(file) }
        .groupTuple(by: 0) // Group by sample_id (index 0)
        .map { sid, sample_nums, lanes, mates, files ->
            // zip the parallel lists: [lane, mate, file]
            def zipped = [lanes, mates, files].transpose()

            // split by mate, sort by lane, collect file paths
            def r1_files = zipped
                .findAll { it[1] as int == 1 }
                .sort    { a, b -> (a[0] as int) <=> (b[0] as int) }
                .collect { it[2] }

            def r2_files = zipped
                .findAll { it[1] as int == 2 }
                .sort    { a, b -> (a[0] as int) <=> (b[0] as int) }
                .collect { it[2] }

            tuple(sid, r1_files, r2_files)
        }

    // Process the samples
    grouped_samples
        .branch { sid, r1_files, r2_files ->
            multi_lane: r1_files.size() > 1 || r2_files.size() > 1
                return tuple(sid, r1_files, r2_files)
            single_lane: true
                return tuple(sid, r1_files[0], r2_files[0])
        }
        .set { samples_to_process }

    // Concatenate multi-lane samples
    CONCAT_FASTQS(samples_to_process.multi_lane)
    
    // Combine results: concatenated multi-lane + original single-lane
    final_reads = CONCAT_FASTQS.out.concatenated
        .mix(samples_to_process.single_lane)

    emit:
    reads = final_reads      // [ sample_id, R1.fastq.gz, R2.fastq.gz ]

}


workflow all {
    // Get the list of sample IDs to filter by
    samples_larry_ch = Channel
        .fromPath(params.sample_csv, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> (row.sample_larry as String).trim() }
        .filter { it }    // drop empties
        .unique()
        .collect()        // collect into a list

    fastq_all_ch = Channel
        .fromPath("${params.fastqs_path}/*_S*_L*_R{1,2}_*.fastq.gz", checkIfExists: true)
        .map { f -> tuple(f.name.split('_S',2)[0], f) }   // (sid, file)

    fastq_larry_ch = fastq_all_ch
        .filter { sid, fq -> 
            samples_larry_ch.val.contains(sid)
        }
        .map { sid, fq -> fq }   // drop sid, keep file

    CONCAT_SAMPLE_FASTQS(fastq_larry_ch)

    samples_all_ch = Channel
        .fromPath(params.sample_csv, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            tuple( (row.sample_larry as String).trim(),
                (row.sample_gex   as String).trim(),
                (row.group_id     as String).trim() )
        }

    merged_larry_fastqs = samples_all_ch
        .join( CONCAT_SAMPLE_FASTQS.out.reads )                 // (sample_larry, gex, group) ⋈ (sample_larry, R1, R2)
        .map { sid, gex, group, r1, r2 ->
            // Re‑order to desired 5‑tuple
            tuple(sid, gex, group, r1, r2)
        }


    FIND_LARRY_SEQS(merged_larry_fastqs)
    
    LARRY_QC(FIND_LARRY_SEQS.out)

    LARRY_QC.out[0].groupTuple(by: 2)
        .set {samples_clones}

    ASSIGN_CLONES(samples_clones)

    MATCH_GEX(ASSIGN_CLONES.out)

}


workflow from_qc {
    Channel.fromPath(params.barm_csv)
        .splitCsv(header: false)
        .map { row -> tuple(row[0], file(row[1])) }
        .set { barm_tuples }

    LARRY_QC(barm_tuples)

    if (params.combine_samples) {
        ASSIGN_CLONES(tuple(LARRY_QC.out[0].collect(){ it[0] }, LARRY_QC.out[0].collect(){ it[1] }))
    } else {
        // Use each output separately
        ASSIGN_CLONES(LARRY_QC.out[0])
    }

    MATCH_GEX(ASSIGN_CLONES.out)

}