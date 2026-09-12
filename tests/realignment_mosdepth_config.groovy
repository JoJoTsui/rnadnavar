// Run with the installed Nextflow *-one.jar on the Java classpath:
// java -cp "$NEXTFLOW_JAR" groovy.ui.GroovyMain tests/realignment_mosdepth_config.groovy
// Evaluates the real config closures and the real module's script block; no jobs run.
import nextflow.config.parser.v1.ConfigParserV1
import nextflow.script.ProcessConfig
import nextflow.script.dsl.ProcessConfigBuilder

def root = new File('.').canonicalFile
def moduleText = new File(root, 'modules/nf-core/mosdepth/main.nf').text
def script = moduleText.split('    script:\\n', 2)[1].split('    stub:', 2)[0]
def configurations = [
    'conf/modules/quality_control/quality_control.config',
    'conf/modules/prepare_realignment/vcf_realignment.config'
]

for (wes in [false, true]) {
    for (realignment in [false, true]) {
        def bed = realignment || wes ? 'targets.bed' : []
        def params = [wes: wes, skip_tools: null, publish_dir_mode: 'copy',
                      save_output_as_bam: false, bam_csi_index: false, save_align_intermeds: false]
        def meta = [id: 'RNA', status: 2]
        def config = new ConfigObject()
        configurations.each { path ->
            config.merge(new ConfigParserV1().setBinding([params: params, meta: meta, bed: bed])
                .parse(new File(root, path)))
        }
        def name = 'NFCORE_RNADNAVAR:RNADNAVAR:' +
            (realignment ? 'RNA_REALIGNMENT_WORKFLOW:' : 'BAM_PROCESSING:') +
            'BAM_GATK_PREPROCESSING:BAM_MARKDUPLICATES:CRAM_QC_MOSDEPTH_SAMTOOLS:MOSDEPTH'
        def processConfig = new ProcessConfig([:])
        new ProcessConfigBuilder(processConfig).applyConfig(config.process, 'MOSDEPTH', 'MOSDEPTH', name)
        def ext = processConfig.ext.collectEntries { key, value ->
            [key, value instanceof Closure ? value.call() : value]
        }
        def binding = new Binding([task: [ext: ext, cpus: 2, process: name],
                                   meta: meta, bed: bed, fasta: 'reference.fa', bam: 'rna.cram',
                                   error: { String message -> throw new IllegalArgumentException(message) }])
        def command = new GroovyShell(binding).evaluate(script).toString()
        assert command.count('--by ') == 1 : command
        assert command.contains(bed ? '--by targets.bed' : '--by 500') : command
        if (!realignment && !wes) {
            assert ext.args == '-n --fast-mode --by 500' : 'First-pass WGS command must stay cache-compatible'
        }
        if (wes) assert ext.args == '' : 'Preserve WES depth semantics'
        println "PASS wes=${wes}, realignment=${realignment}: ${ext.args}, BED=${bed}"
    }
}
