rule draw_stacked_ig_manhattans:
    input:
        iga = "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/filtered_meta.tsv.gz",
        igg = "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/filtered_meta.tsv.gz",
        igm = "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/filtered_meta.tsv.gz"
    output:
        "results/paper/figures/ig_manhattans.png"
    params:
        width = 12,
        height = 14,
        iga_ylim = [1, 1e-60],
        igg_ylim = [1, 1e-40],
        igm_ylim = [1, 1e-60]
    threads: 12
    resources:
        runtime = 20
    conda: env_path("global.yaml")
    script: script_path("paper/figures/draw_stacked_ig_manhattans.R")

rule candidate_ig_ig_coloc_plots:
    input:
        expand("results/coloc/iga_and_igm/{rsid_pair}/trimmed/lz_plots.png", rsid_pair = ["rs188468174_and_rs188468174", "rs6570530_and_rs7758383", "rs35969813_and_rs35629860", "rs2872516_and_rs4795397", "rs6885567_and_rs6879652", "rs10517538_and_rs1397934"]),
        "results/coloc/iga_and_igg/rs3803800_and_rs758641530/trimmed/lz_plots.png",
        "results/coloc/iga_and_igg/rs3803800_and_rs758641530/untrimmed/lz_plots.png",
        "results/coloc/raw_vs_trimmed_coloc_pp_h4_ig_pairs.png"
    output:
        expand("export/coloc/iga_and_igm/{rsid_pair}/trimmed/lz_plots.png", rsid_pair = ["rs188468174_and_rs188468174", "rs6570530_and_rs7758383", "rs35969813_and_rs35629860", "rs2872516_and_rs4795397", "rs6885567_and_rs6879652", "rs10517538_and_rs1397934"]),
        "export/coloc/iga_and_igg/rs3803800_and_rs758641530/trimmed/lz_plots.png",
        "export/coloc/iga_and_igg/rs3803800_and_rs758641530/untrimmed/lz_plots.png",
        "export/coloc/raw_vs_trimmed_coloc_pp_h4_ig_pairs.png"
    localrule: True
    run:
        for i,x in enumerate(input):
            y = output[i]
            shell(f"cp {x} {y}")

use rule candidate_ig_ig_coloc_plots as candidate_ig_lymphocyte_counts_coloc_plots with:
    input:
        expand("results/coloc/iga_and_lymphocyte-counts/{rsid}/lz_plots.png", rsid = ["rs3757387", "rs1003342", "rs2256609", "rs2286564", "rs6570530", "rs2872516", "rs3810277", "rs148076268", "rs3743"]),
        expand("results/coloc/igm_and_lymphocyte-counts/{rsid}/lz_plots.png", rsid = ["rs2476601", "rs6891054", "rs4795397", "rs735665", "rs9293515", "rs28671336", "rs6495122", "rs290243"]),
        "results/coloc/igg_and_lymphocyte-counts/rs1260326/lz_plots.png"
    output:
        expand("export/coloc/iga_and_lymphocyte-counts/{rsid}/lz_plots.png", rsid = ["rs3757387", "rs1003342", "rs2256609", "rs2286564", "rs6570530", "rs2872516", "rs3810277", "rs148076268", "rs3743"]),
        expand("export/coloc/igm_and_lymphocyte-counts/{rsid}/lz_plots.png", rsid = ["rs2476601", "rs6891054", "rs4795397", "rs735665", "rs9293515", "rs28671336", "rs6495122", "rs290243"]),
        "export/coloc/igg_and_lymphocyte-counts/rs1260326/lz_plots.png"

rule ig_ig_coloc_plots:
    input:
        expand("results/coloc/iga_and_igm/{rsids}/{filtering}/lz_plots.png",
               rsids = [f'{x['iga']}_and_{x['igm']}' for x in config['coloc']['iga_and_igm']['hits']],
               filtering = ['trimmed', 'raw']
              )

rule ig_lymphocyte_counts_coloc_plots:
    input:
        expand("results/coloc/iga_and_lymphocyte-counts/{rsid}/lz_plots.png",
        rsid = [x['iga'] for x in config['coloc']['iga_and_lymphocyte-counts']['hits']]
               ),
        expand("results/coloc/igm_and_lymphocyte-counts/{rsid}/lz_plots.png",
        rsid = [x['igm'] for x in config['coloc']['igm_and_lymphocyte-counts']['hits']]
               ),
        expand("results/coloc/igg_and_lymphocyte-counts/{rsid}/lz_plots.png",
        rsid = [x['igg'] for x in config['coloc']['igg_and_lymphocyte-counts']['hits']]
               )

