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

  script:
  """
  python ${baseDir}/bin/larry_qc.py ${pkl} ${samp} "${params.check_whitelist}" ${params.plot_qc}
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
  python ${baseDir}/bin/match_gex.py ${samp} ${pkl} ${params.sample_list} ${params.gex_path} "${params.combine_samples}" "${params.plot_cumulative}"
  """
}


workflow {
    def samples = new groovy.json.JsonSlurper().parseText(file(params.sample_list).text).keySet() as List

    Channel.fromFilePairs("${params.fastqs_path}/*_S1_L001_R[12]_001.fastq.gz", flat: true)
        .filter { pair -> 
            def (sample, files) = pair
            samples.contains(sample)
        }
        .map { it -> it }
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