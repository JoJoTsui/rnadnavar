//
// CHANNEL_VARIANT_CALLING_CREATE_CSV
//

workflow CHANNEL_VARIANT_CALLING_CREATE_CSV {
    take:
        vcf_to_csv // channel: [mandatory] meta, vcf
        csv_name

    main:
        // Creating csv files to restart from this step
        vcf_to_csv.collectFile(keepHeader: true, skip: 1,sort: true, storeDir: "${params.outdir}/csv"){ meta, vcf ->
            def patient       = meta.patient
            def sample        = meta.id
            def variantcaller = meta.variantcaller
            def status        = meta.status
            // Use different path structure for normalized files
            if (csv_name == "normalized") {
                vcf = "${params.outdir}/normalized/${variantcaller}/${meta.id}/${vcf.getName()}"
            } else if (csv_name == "annotated") {
                vcf = "${params.outdir}/annotation/${variantcaller}/${meta.id}/${vcf.getName()}"
            } else {
                vcf = "${params.outdir}/variant_calling/${variantcaller}/${meta.id}/${vcf.getName()}"
            }
            ["${csv_name}.csv", "patient,sample,status,variantcaller,vcf\n${patient},${sample},${status},${variantcaller},${vcf}\n"]
        }
}
