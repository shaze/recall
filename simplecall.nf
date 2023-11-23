
// (c) University of the Witwatersrand, Johannesburg
// Author Scott Hazelhurst
// Released under https://creativecommons.org/licenses/by-sa/4.0/


nextflow.enable.dsl=2

params.output = "sample"
params.has_bam = false
trf = file(params.tandem_example)

ref_seq=file(params.ref_seq)
ref_fai=file(params.ref_fai)
ref_mmi=file(params.ref_mmi)

params.chrom_prefix="chr"
params.input_type="bam"


autosomes = 1..22
params.chroms    = autosomes +  ['X','Y','M'] 


chroms=params.chroms
chroms = 15..22

// Do an initial calling of the variants
// Note that this is done per chromosome
// -- so for each sample of input we have multiple processes running
process deepcall {
  maxForks 80
  label 'deepvariant'
  cpus params.call_cpus
  memory params.call_mem
  errorStrategy 'finish'
  if (params.constraint)
     clusterOptions="--constraint=${params.constraint}"
  input:
     tuple val(the_id), file(bam), file(bai) 
     file ref_seq
     file ref_fai
     each chroms
  output:
    tuple val(chroms), \
	file("${vcf}.gvcf.gz"), file("${vcf}.gvcf.gz.tbi"),\
	file(bam), file(bai)
    publishDir "${params.output_dir}/vcf/$c", mode: 'copy', pattern: "*.vcf.gz"
    publishDir "${params.output_dir}/gvcf/$c", mode: 'copy', pattern: "*.gvcf.gz"    
  script:
     name = bam.simpleName
     if (['X','Y','M'].contains(chroms))
        c=chroms
     else
         c = "${chroms}".padLeft(2,"0")
         vcf = "${name}-${c}"
         tbi = "${vcf}.vcf.gz.tbi"
     """
      hostname
      /opt/deepvariant/bin/run_deepvariant  \
                --model_type WGS \
                --ref $ref_seq  \
                --reads  $bam  \
                --output_vcf=${vcf}.vcf.gz \
                --output_gvcf=${vcf}.gvcf.gz \
                --num_shards ${params.call_cpus} \
                --regions chr$chroms
      """
}


process GLnexus {
    input:
      tuple val(chrom), file(gvcfsetal)
    output:
       file("${c}.bcf")  
    publishDir "$params.out/vcf/"
    script:
    c = "${chroms}".padLeft(2,"0")
    """
    glnexus_cli –config DeepVariantWGS *vcf.gz > ${c}.bcf
    """
}



workflow {
    bf = Channel.fromFilePairs("${params.bam}*${params.input_type}*",size:2) \
	    { file -> file.simpleName }.\
	    map { b, f -> [b,f[0],f[1]] }
    deepcall(bf,ref_seq,ref_fai,chroms) //| groupTuple | GLnexus
}
