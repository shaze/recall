


chrom_size_b38 = ['01':248956422, '02':242193529, '03':198295559, '04':190214555,
	      '05':181538259, '06':170805979, '07':159345973, '08':145138636,
	      '09':138394717, '10':133797422, '11':135086622, '12':133275039,
	      '13':114364328, '14':107043718, '15':101991180, '16':90338345,
	      '17':83257441,  '18':80373285,  '19':50617616,  '20':64444167,
	      '21':47709983,  '22':50818468,  'X' :156040896, 'Y' :57227415,
	      'M':17000]

chrom_size_t2t = ['01': 248387328, '02': 242696752, '03': 201105948, '04': 193574945,
	      '05': 182045439, '06': 172126628, '07': 160567428, '08': 146259331,
	      '09': 150617247, '10': 134758134, '11': 135127769, '12': 133324548,
	      '13': 113566686, '14': 101161492, '15': 99753195, '16': 96330374,
	      '17': 84276897, '18': 80542538, '19': 61707364, '20': 66210255,
	      '21': 45090682, '22': 51324926, 'X': 154259566, 'Y': 62460029, 'M': 16569]

chrom_size = chrom_size_t2t

def memoryForChrom ( the_chrom ) {
    chrom =  the_chrom.getName()
    req   = Math.round(chrom_size[chrom]*2.96/1000000+246)
    return "${req}.G"
}



if (params.scratch=="false")
    scratch=false;
else if (params.scratch=="true")
    scratch==true;
else
    scratch=params.scratch;


process GLnexus {
    scratch scratch
    memory { memoryForChrom(vdir) }
    cpus   params.glnexus_cpus
    input:
      path(vdir)
    output:
       path(bcf), emit: raw
       path(res), emit: times
    publishDir params.output_dir, pattern: "*.t"
    script:
       bcf="dv${vdir}.bcf"
       res="dv${vdir}.t"
       """
       hostname 
       set limit descriptors 3600
       /usr/bin/time -f "Elapsed time: %e; MAXRSS=%M" -o $res\
            glnexus_cli  --config DeepVariantWGS $vdir/*vcf.gz > $bcf
       """
}

