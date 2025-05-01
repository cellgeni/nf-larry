#!/usr/bin/env nextflow

process FIND_LARRY_SEQS {

  publishDir "${launchDir}/larry-results-${params.project_tag}/before_qc/", mode: 'copy'

  input:
  tuple val(samp), path(r1_file), path(r2_file)

  output:
  tuple val(samp), path("*.pkl")

  script:
  """
  python ${baseDir}/bin/find_larry_seqs.py ${params.fastqs_path} ${r1_file} ${r2_file} ${samp}
  """
}

process LARRY_QC {

  publishDir "${launchDir}/larry-results-${params.project_tag}/after_qc/", mode: 'copy'

  input:
  tuple val(samp), path(pkl)

  output:
  tuple val(samp), path("*.pkl")
  path("*.pdf"), optional: true

  script:
  """
  python ${baseDir}/bin/larry_qc.py ${pkl} ${samp} ${params.check_whitelist ? '--check_whitelist' : ''} ${params.make_pdf ? '--make_pdf' : ''}
  """
}

process ASSIGN_CLONES {

  publishDir "${launchDir}/larry-results-${params.project_tag}/clones/", mode: 'copy', saveAs: {filename -> "${samp}_${filename}"}

  input:
  tuple val(samp), path(pkl)

  output:
  tuple val(samp), path("*.pkl")

  script:
  """
  python ${baseDir}/bin/assign_clones.py ${pkl.join(",")} "${params.combine_samples}" ${params.dispr_filter}
  """
}

process MATCH_GEX {

  publishDir "${launchDir}/larry-results-${params.project_tag}/clones/", mode: 'copy'

  input:
  tuple val(samp), path(pkl)

  output:
  tuple val(samp), path("*.h5ad"), path("*.csv")
  path("*.png"), optional: true

  script:
  """
  python ${baseDir}/bin/match_gex.py ${samp} ${pkl} ${params.sample_json} ${params.gex_path} ${params.combine_samples ? '--combine_samples' : ''} ${params.plot_cumulative ? '--plot_cumulative' : ''}
  """
}


workflow all {
 // Read JSON file and extract keys as a Groovy list
    def jsonFile = file(params.sample_json)
    def jsonText = jsonFile.text
    def sample_keys = new groovy.json.JsonSlurper().parseText(jsonText).keySet() as List

    // Now sample_keys contains: ['sample1_larry', 'sample2_larry', ...]
    println(sample_keys)

    // Example: filter file pairs using these keys
    ch_pairs = Channel.fromFilePairs("${params.fastqs_path}/(.+)_S1_L001_R[12]_001.fastq.gz", flat: true)
        .filter { pair -> 
            def (sample, files) = pair
            sample_keys.contains(sample)
        }
        .set { filtered_pairs }

    FIND_LARRY_SEQS(filtered_pairs)
    LARRY_QC(FIND_LARRY_SEQS.out)

    if (params.combine_samples) {
        ASSIGN_CLONES(tuple(LARRY_QC.out.collect(){ it[0] }, LARRY_QC.out.collect(){ it[1] }))
    } else {
        // Use each output separately
        ASSIGN_CLONES(LARRY_QC.out)
    }

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