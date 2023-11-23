

nextflow.enable.dsl=2

gvcfs = Channel.fromPath(params.gvcfs,type:'dir')


params.build="Specify build"   // b37, b38, t2t
params.qual=20
params.depth=25
params.max_size = 80000000
params.output_dir="output"

max_size = params.max_size

ref=file(params.ref_seq)

include { GLnexus } from "./modules/glnexus.nf"




process filter {
   input:
      file(bcf) 
   output:
      file(filterbcf) 
   script:
      filterbcf="${bcf}.filter.bcf"
   """
     hostname
     bcftools filter -i 'QUAL>${params.qual} & FORMAT/DP>${params.depth}' $bcf  -o $filterbcf 
   """
}

process norm {
    input:
      tuple path(bcf), path(ref)
   output:
      file(normbcf) 
   publishDir params.output_dir
   script:
      normbcf="${bcf.simpleName}.bcf"
   """
     bcftools norm --fasta-ref  $ref $bcf  -o $normbcf 
   """
}

process index {
   input:
      file(bcf)  
   output:
      file(index)
   publishDir params.output_dir
   script:
      index="${bcf}.csi"
   """
     bcftools index $bcf
   """
}


workflow {
    main:
       gvcfs = Channel.fromPath(params.gvcfs,type:'dir').view()
       refc = Channel.fromPath(ref)
       GLnexus(gvcfs)
       filter(GLnexus.out.raw) | combine(refc) | norm | index 
}
