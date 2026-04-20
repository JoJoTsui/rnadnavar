process ENSEMBLVEP_VEP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ensembl-vep:113.4--pl5321h2a3209d_0' :
        'biocontainers/ensembl-vep:113.4--pl5321h2a3209d_0' }"

    input:
    tuple val(meta), path(vcf), path(custom_extra_files)
    val   genome
    val   species
    val   cache_version
    path  cache
    tuple val(meta2), path(fasta)
    path  extra_files

    output:
    tuple val(meta), path("*.vcf.gz")       , optional:true, emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi")   , optional:true, emit: tbi
    tuple val(meta), path("*.tab.gz")       , optional:true, emit: tab
    tuple val(meta), path("*.json.gz")      , optional:true, emit: json
    path "*.html"                           , optional:true, emit: report
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def file_extension = args.contains("--vcf") ? 'vcf' : args.contains("--json")? 'json' : args.contains("--tab")? 'tab' : 'vcf'
    def compress_cmd = args.contains("--compress_output") ? '' : '--compress_output bgzip'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dir_cache = cache ? "\${PWD}/${cache}" : "/.vep"
    def reference = fasta ? "--fasta $fasta" : ""
    def create_index = file_extension == "vcf" ? "tabix ${args2} ${prefix}.${file_extension}.gz" : ""
    def plugin_code = 'package SubstrSafe;use strict;use warnings;use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);my $patched=0;sub _wrap_safe{my($pkg,$func_name)=@_;my $full_name="${pkg}::${func_name}";my $orig=do{no strict "refs";\\&{$full_name}};no strict "refs";no warnings "redefine";*{$full_name}=sub{my $result=eval{$orig->(@_)};if($@){my $var_info="unknown";eval{my($bvfoa,$feat,$bvfo,$bvf)=@_;$bvfo||=$bvfoa->base_variation_feature_overlap;my $vf=$bvfo->variation_feature;my $chr=$vf->seq_region_name||"?";my $start=$vf->start||0;my $alleles=$vf->allele_string||"?";my $tr=$bvfo->transcript;my $tr_id=$tr?$tr->stable_id:"?";$var_info="${chr}:${start} ${alleles} transcript=${tr_id}"};warn "SubstrSafe: CAUGHT error in ${func_name} at ${var_info}: $@";return 0}return $result}}sub new{my $class=shift;my $self=$class->SUPER::new(@_);unless($patched){require Bio::EnsEMBL::Variation::Utils::VariationEffect;my $pkg="Bio::EnsEMBL::Variation::Utils::VariationEffect";for my $func(qw(ref_eq_alt_sequence stop_retained stop_lost frameshift inframe_insertion inframe_deletion)){eval{_wrap_safe($pkg,$func)}}$patched=1}return $self}sub feature_types{return["Transcript"]}sub get_header_info{return{}}sub run{return{}}1;'
    """
    mkdir -p vep_plugins
    echo '${plugin_code}' > vep_plugins/SubstrSafe.pm

    vep \\
        -i $vcf \\
        -o ${prefix}.${file_extension}.gz \\
        $args \\
        $compress_cmd \\
        $reference \\
        --assembly $genome \\
        --species $species \\
        --cache \\
        --cache_version $cache_version \\
        --dir_cache $dir_cache \\
        --dir_plugins vep_plugins \\
        --plugin SubstrSafe \\
        --fork $task.cpus

    ${create_index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_extension = args.contains("--vcf") ? 'vcf' : args.contains("--json")? 'json' : args.contains("--tab")? 'tab' : 'vcf'
    def create_index = file_extension == "vcf" ? "touch ${prefix}.${file_extension}.gz.tbi" : ""
    """
    echo "" | gzip > ${prefix}.${file_extension}.gz
    ${create_index}
    touch ${prefix}_summary.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
