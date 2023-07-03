nextflow.enable.dsl=1

gvcfA = Channel.fromPath("/external/diskB/recall/gvcf/[012M]*",type:'dir')
gvcfS = Channel.fromPath("/external/diskB/recall/gvcf/[XY]",type:'dir')

params.qual=20
params.depth=25
params.ref_seq = "/dataB/aux/38/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta"
params.output_dir="output"

ref=file(params.ref_seq)



def memoryForChrom ( the_chrom ) {
    chrom =  the_chrom.getName()
    if (['X','01','02','03'].contains(chrom))
	{ mem = "900.GB" }
    else if (['14','15','16','17','18'].contains(chrom))
	{ mem =  "500.GB" }
    else if (['19','20','21','22'].contains(chrom))
	{ mem = "400.GB" }
    else if (['M'].contains(chrom))
	{ mem = "200.GB" }
    else if (['Y'].contains(chrom))
	{ mem =  "250.GB" }
    else
	{mem = "700.GB" }
    return mem
}


process GLnexusA {
    scratch '/local/scott/$vdir'
    memory { memoryForChrom(vdir) }
    cpus   40
    input:
      path(vdir) from gvcfA
    output:
       path(bcf) into vcfA_ch
       path(res) into timeA_ch
    publishDir params.output_dir, pattern: "*.t"
    script:
       bcf="dv${vdir}.bcf"
       res="dv${vdir}.t"
       """
       hostname 
       /usr/bin/time -f "Elapsed time: %e; MAXRSS=%M" -o $res\
            /opt/exp_soft/bioinf/GLnexus/1.4.1-220324/glnexus_cli  --config DeepVariantWGS $vdir/*vcf.gz > $bcf
       """
}



process GLnexusS {
    scratch '/local/scott/$vdir'
    memory { memoryForChrom(vdir) }
    cpus   40
    input:
      path(vdir) from gvcfS
    output:
       path(bcf) into vcfS_ch
       path(res) into timeS_ch
    publishDir params.output_dir, pattern: "*.t"
    script:
       bcf="dv${vdir}.bcf"
       res="dv${vdir}.t"
       """
       hostname 
       /usr/bin/time -f "Elapsed time: %e; MAXRSS=%M" -o $res\
            /opt/exp_soft/bioinf/GLnexus/1.4.1-220324/glnexus_cli  --config DeepVariantWGS $vdir/*vcf.gz > $bcf
       """
}


vcf_ch = vcfS_ch.mix(vcfA_ch)
time_ch = timeS_ch.mix(timeA_ch)


process filter {
   input:
      file(bcf) from vcf_ch
   output:
      file(filterbcf) into filter_ch
   script:
      filterbcf="${bcf}.filter.bcf"
   """
     hostname
     bcftools filter -i 'QUAL>${params.qual} & FORMAT/DP>${params.depth}' $bcf  -o $filterbcf 
   """
}

process norm {
    input:
      file(ref)
      file(bcf)  from filter_ch
   output:
      file(normbcf) into norm_ch
   publishDir params.output_dir
   script:
      normbcf="${bcf.simpleName}.bcf"
   """
     bcftools norm --fasta-ref  $ref $bcf  -o $normbcf 
   """
}

process index {
   input:
      file(bcf)  from norm_ch
   output:
      file(index)
   publishDir params.output_dir
   script:
      index="${bcf}.csi"
   """
     bcftools index $bcf
   """
}
