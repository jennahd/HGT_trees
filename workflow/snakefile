num=["0","1","2"]

rule all:
    input:
        #expand("{workdir}/ALL.sister_nested_taxonomy.tsv", workdir=config["workdir"]),
        expand("{workdir}/ALL.sister_nested_taxonomy.pdf", workdir=config["workdir"]),
        expand("{workdir}/7.subtree/{cluster}.{num}.subtree.align.trim.treefile.pdf", cluster=config["clusters"], workdir=config["workdir"], num=num)
        #expand("{workdir}/7.subtree/{cluster}.{num}.sister_nested_taxonomy.tsv", cluster=config["clusters"], workdir=config["workdir"], num=num)
    params:
        workdir=config["workdir"]
    run:
        shell("find {params.workdir}/5.select_subtree/*keep*.pdf -type f -empty -delete")
        shell("find {params.workdir}/6.subtree_faa -type f -empty -delete")
        shell("find {params.workdir}/7.subtree -type f -empty -delete")
        shell("find {params.workdir}/8.subtree_interproscan -type f -empty -delete")

###Find top blast hits of sequences of interest###

rule make_lists:
    input:
        "{workdir}/faa/{cluster}.faa"
    output:
        "{workdir}/faa/{cluster}.list"
    shell:
        "grep '>' {input} | cut -d '>' -f2 > {output}"

rule combine_faa:
    input:
        expand("{workdir}/faa/{cluster}.faa", cluster=config["clusters"], workdir=config["workdir"])
    output:
        "{workdir}/faa/all_proteins.faa"
    shell:
        "cat {input} > {output}"

rule diamond_blast:
    input:
        db="workflow/databases/nr.dmnd",
        faa=ancient("{workdir}/faa/all_proteins.faa")
    output:
        "{workdir}/1.diamond_blast/all_vs_nr_diamond_blast.tsv"
    threads: workflow.cores
    log:
        "{workdir}/logs/1.diamond_blast/all_vs_nr_diamond_blast.log"
    shell:
        "diamond blastp "
        "--threads {threads} "
        "--db {input.db} "
        "--out {output} "
        "--outfmt 6 qseqid qlen sseqid slen stitle length qcovhsp pident evalue bitscore "
        "--query {input.faa} "
        "--max-target-seqs 2000 "
        "--more-sensitive > {log}"

rule split_faa:
    input:
        blast=ancient("{workdir}/1.diamond_blast/all_vs_nr_diamond_blast.tsv"),
        list="{workdir}/faa/{cluster}.list"
    output:
        "{workdir}/1.diamond_blast/{cluster}.diamond_blast.tsv"
    run:
        try:
            shell("grep -w -f {input.list} {input.blast} > {output}")
        except:
            shell("touch {output}")

rule unique_hits:
    input:
        ancient("{workdir}/1.diamond_blast/{cluster}.diamond_blast.tsv")
    output:
        "{workdir}/1.diamond_blast/{cluster}.diamond_blast.uniq.list"
    run:
        try:
            shell("cut -f3 {input} | sort | uniq > {output}")
        except:
            shell("touch {output}")

rule retrive_nr_fasta:
    input:
        db="workflow/databases/nr.gz",
        list="{workdir}/1.diamond_blast/{cluster}.diamond_blast.uniq.list"
    output:
        "{workdir}/2.nr_faa/{cluster}.nr.faa"
    run:
        try:
            shell("seqtk subseq <(zcat {input.db}) {input.list} > {output}")
        except:
            shell("touch {output}")

###Remove sequence redundancy and get taxonomy of blast hits###

rule cd_hit:
    input:
        "{workdir}/2.nr_faa/{cluster}.nr.faa"
    output:
        "{workdir}/2.nr_faa/{cluster}.nr.cd-hit.faa"
    threads: 4
    params:
        cd_hit_perc_id = config["cd_hit_perc_id"]
    log:
        "{workdir}/logs/2.nr_faa/{cluster}.nr.cd-hit.log"
    run:
        try:
            shell("cd-hit -i {input} -o {output} -c {params.cd_hit_perc_id} -T {threads} -n 5 > {log}")
        except:
            shell("touch {output}")

rule get_accessions:
    input:
        "{workdir}/2.nr_faa/{cluster}.nr.cd-hit.faa"
    output:
        "{workdir}/2.nr_faa/{cluster}.nr.cd-hit.list"
    run:
        try:
            shell("grep '>' {input} | cut -d '>' -f2 | cut -d ' ' -f1 > {output}")
        except:
            shell("touch {output}")

rule get_taxids:
    input:
        "{workdir}/2.nr_faa/{cluster}.nr.cd-hit.list"
    output:
        taxids="{workdir}/2.nr_faa/{cluster}.nr.cd-hit.taxids.tsv",
        notaxids="{workdir}/2.nr_faa/{cluster}.nr.cd-hit.no-taxonomy"
    run:
        shell("touch {output.taxids} {output.notaxids}")
        shell("while read p ; do blastdbcmd -db /local/three/databases/ncbi_nr_v5/nr "
            "-dbtype prot -entry $p -target_only -outfmt '%a\t%T' >> {output.taxids}"
            " || echo $p >> {output.notaxids}; done <{input}")

rule add_taxonomy:
    input:
        "workflow/databases/fullnamelineage.dmp",
        "{workdir}/2.nr_faa/{cluster}.nr.cd-hit.taxids.tsv",
        "{workdir}/2.nr_faa/{cluster}.nr.cd-hit.faa"
    output:
        "{workdir}/2.nr_faa/{cluster}.nr.cd-hit.tax.faa"
    script:
        "scripts/add_taxonomy.py"

rule combine_files:
    input:
        cluster="{workdir}/faa/{cluster}.faa",
        nr="{workdir}/2.nr_faa/{cluster}.nr.cd-hit.tax.faa"
    output:
        "{workdir}/3.combined_faa/{cluster}.comb.faa"
    run:
        try:
            shell("cat {input.cluster} {input.nr} > {output}")
        except:
            shell("cat {input.cluster} > {output}")

rule cd_hit2:
    input:
        "{workdir}/3.combined_faa/{cluster}.comb.faa"
    output:
        "{workdir}/3.combined_faa/{cluster}.comb-99.faa"
    threads: 4
    log:
        "{workdir}/logs/3.combined_faa/{cluster}.comb-99.log"
    shell:
        "cat {input} > {output}"
        #"cd-hit -i {input} -o {output} -c 0.99 -T {threads} -n 5 > {log}"

###Initial alignment and tree, removal of both short and long-branching sequences###

rule inital_alignment:
    input:
        "{workdir}/3.combined_faa/{cluster}.comb-99.faa"
    output:
        "{workdir}/4.remove_long_branches/{cluster}.align"
    threads: 4
    shell:
        "mafft --anysymbol --thread {threads} --reorder --auto {input} > {output}"

rule inital_trimming:
    input:
        "{workdir}/4.remove_long_branches/{cluster}.align"
    output:
        "{workdir}/4.remove_long_branches/{cluster}.align.trim"
    shell:
        "trimal -in {input} -out {output} -gappyout"

rule remove_short_seqs:
    input:
        "{workdir}/4.remove_long_branches/{cluster}.align.trim",
        config["mapping_file"]
    output:
        "{workdir}/4.remove_long_branches/{cluster}.align.trim.len"
    params:
        "40"
    script:
        "scripts/remove_seqs_alignment_length.py"

rule inital_fasttree:
    input:
        "{workdir}/4.remove_long_branches/{cluster}.align.trim.len"
    output:
        "{workdir}/4.remove_long_branches/{cluster}.align.trim.tree"
    log:
        "{workdir}/logs/4.remove_long_branches/{cluster}.align.tree.log"
    shell:
        "FastTree -log {log} {input} > {output}"

rule remove_long_branches:
    input:
        map = config["mapping_file"],
        info = config["dataset_chacteristics"],
        tree = "{workdir}/4.remove_long_branches/{cluster}.align.trim.tree"
    output:
        out = "{workdir}/4.remove_long_branches/{cluster}.long-branches.tree.pdf",
        long = "{workdir}/4.remove_long_branches/{cluster}.long-branches.list",
        keep = "{workdir}/5.select_subtree/{cluster}.keep.list"
    params:
        cluster_name = "{cluster}"
    script:
        "scripts/remove_long-branches.py"

###Extract sequences retained and re-run tree before selection of focal groups of taxa of interest and the subtree###

rule extract_seqs:
    input:
        list = "{workdir}/5.select_subtree/{cluster}.keep.list",
        db = "{workdir}/3.combined_faa/{cluster}.comb-99.faa"
    output:
        "{workdir}/5.select_subtree/{cluster}.keep.faa"
    run:
        shell("touch {output}")
        shell("seqtk subseq {input.db} {input.list} > {output}")

rule alignment:
    input:
        "{workdir}/5.select_subtree/{cluster}.keep.faa"
    output:
        "{workdir}/5.select_subtree/{cluster}.keep.align"
    threads: 4
    run:
        try:
            shell("mafft --anysymbol --thread {threads} --reorder --auto {input} > {output}")
        except:
            shell("touch {output}")

rule trimming:
    input:
        "{workdir}/5.select_subtree/{cluster}.keep.align"
    output:
        "{workdir}/5.select_subtree/{cluster}.keep.align.trim"
    run:
        try:
            shell("trimal -in {input} -out {output} -gappyout")
        except:
            shell("touch {output}")

rule fasttree:
    input:
        "{workdir}/5.select_subtree/{cluster}.keep.align.trim"
    output:
        "{workdir}/5.select_subtree/{cluster}.keep.align.trim.tree"
    log:
        "{workdir}/logs/5.select_subtree/{cluster}.keep.align.trim.tree.log"
    run:
        try:
            shell("grep -c '>' {input}")
        except:
            shell("touch {output}")
        else:
            shell("FastTree -log {log} {input} > {output}")

rule show_focal_nodes:
    input:
        map = config["mapping_file"],
        info = config["dataset_chacteristics"],
        tree = "{workdir}/5.select_subtree/{cluster}.keep.align.trim.tree"
    output:
        overview = "{workdir}/5.select_subtree/{cluster}.focal_nodes.tree.pdf",
        taxa = "{workdir}/5.select_subtree/{cluster}.focal_nodes.taxa.tsv"
    params:
        cluster_name = "{cluster}",
        focal_taxonomy = config["focal_taxonomy"]
    script:
        "scripts/show_focal_nodes.py"

rule select_subtree:
    input:
        map = config["mapping_file"],
        info = config["dataset_chacteristics"],
        tree = "{workdir}/5.select_subtree/{cluster}.keep.align.trim.tree",
        taxa = "{workdir}/5.select_subtree/{cluster}.focal_nodes.taxa.tsv"
    output:
        out = "{workdir}/5.select_subtree/{cluster}.{num}.keep.subtree.tree.pdf",
        seqs = "{workdir}/6.subtree_faa/{cluster}.{num}.seqs.list",
        focal = "{workdir}/6.subtree_faa/{cluster}.{num}.focal.list",
        outgroup = "{workdir}/6.subtree_faa/{cluster}.{num}.outgroup.list",
        removed = "{workdir}/6.subtree_faa/{cluster}.{num}.removed.list",
        summary = "{workdir}/5.select_subtree/{cluster}.{num}.summary.tsv"
    params:
        cluster_name = "{cluster}",
        num = "{num}",
        min_taxa = config["min_taxa"],
        max_taxa = config["max_taxa"],
        focal_taxonomy = config["focal_taxonomy"]
    script:
        "scripts/select_subtree.py"

rule combine_subtree:
    input:
        expand("{workdir}/5.select_subtree/{cluster}.{num}.summary.tsv", cluster=config["clusters"], workdir=config["workdir"], num=num)
    output:
        "{workdir}/summary_subtree_selection.tsv"
    shell:
        "cat {input} > {output}"

rule extract_focal:
    input:
        list = "{workdir}/6.subtree_faa/{cluster}.{num}.focal.list",
        db = "{workdir}/3.combined_faa/{cluster}.comb-99.faa"
    output:
        "{workdir}/6.subtree_faa/{cluster}.{num}.focal.faa"
    shell:
        "seqtk subseq {input.db} {input.list} > {output}"

rule extract_subset:
    input:
        list = "{workdir}/6.subtree_faa/{cluster}.{num}.seqs.list",
        db = "{workdir}/3.combined_faa/{cluster}.comb-99.faa"
    output:
        "{workdir}/6.subtree_faa/{cluster}.{num}.seqs.faa"
    shell:
        "seqtk subseq {input.db} {input.list} > {output}"

rule extract_outgroup:
    input:
        list = "{workdir}/6.subtree_faa/{cluster}.{num}.outgroup.list",
        db = "{workdir}/3.combined_faa/{cluster}.comb-99.faa"
    output:
        "{workdir}/6.subtree_faa/{cluster}.{num}.outgroup.faa"
    run:
        shell("touch {output}")
        shell("seqtk subseq {input.db} {input.list} > {output}")

rule extract_combine:
    input:
        focal="{workdir}/6.subtree_faa/{cluster}.{num}.focal.faa",
        seqs="{workdir}/6.subtree_faa/{cluster}.{num}.seqs.faa",
        outgroup="{workdir}/6.subtree_faa/{cluster}.{num}.outgroup.faa"
    output:
        "{workdir}/6.subtree_faa/{cluster}.{num}.subtree.faa"
    params:
        num = "{num}"
    run:
        try:
            shell("grep -c '>' {input.focal}")
        except:
            shell("touch {output}")
        else:
            shell("sed -i 's/>/>OUTGROUP_{params.num}__/g' {input.outgroup}")
            shell("sed -i 's/>/>FOCAL_{params.num}__/g' {input.focal}")
            shell("cat {input.focal} {input.seqs} {input.outgroup} > {output} ")

rule get_phyla_list:
    input:
        expand("{workdir}/6.subtree_faa/{cluster}.{num}.subtree.faa", cluster=config["clusters"], workdir=config["workdir"], num=num)
    output:
        "{workdir}/7.subtree/all_subset_phyla.list"
    run:
        try:
            shell("grep '>' {input} | cut -d '|' -f3 | cut -d '@' -f1,2,3 | sort | uniq  -c | sort -rn |"
                "grep -v 'na' | grep -v 'environmental_samples' | grep '@[A-Z]' | grep -v 'unclassified'"
                "| grep -v 'incertae_sedis' | tr -s ' ' > {output}")
        except:
            shell("touch {output}")

rule subtree_alignment:
    input:
        faa = "{workdir}/6.subtree_faa/{cluster}.{num}.subtree.faa"
    output:
        "{workdir}/7.subtree/{cluster}.{num}.subtree.align"
    threads: 4
    run:
        try:
            shell("grep -c '>' {input}")
        except:
            shell("touch {output}")
        else:
            shell("mafft --anysymbol --thread {threads} --reorder --auto {input} > {output}")

rule subtree_trimming:
    input:
        "{workdir}/7.subtree/{cluster}.{num}.subtree.align"
    output:
        "{workdir}/7.subtree/{cluster}.{num}.subtree.align.trim"
    run:
        try:
            shell("trimal -in {input} -out {output} -gappyout")
        except:
            shell("touch {output}")

rule subtree:
    input:
        "{workdir}/7.subtree/{cluster}.{num}.subtree.align.trim"
    output:
        "{workdir}/7.subtree/{cluster}.{num}.subtree.align.trim.treefile"
    threads: 4
    log:
        "{workdir}/logs/7.subtree/{cluster}.{num}.subtree.align.trim.tree.log"
    run:
        try:
            shell("grep -c '>' {input}")
        except:
            shell("touch {output}")
        else:
            #if config["iqtree"] == "yes":
            #    shell("iqtree "
            #            "-s {input} "
            #            "-m LG"
            #            "-bb 1000 "
            #            "-nt {threads} ")
            #else:
            shell("FastTree -log {log} {input} > {output}")

rule interproscan_domains:
    input:
        "{workdir}/6.subtree_faa/{cluster}.{num}.subtree.faa"
    output:
        "{workdir}/8.subtree_interproscan/{cluster}.{num}.subtree.interproscan.tsv"
    threads: 4
    log:
        "{workdir}/logs/8.subtree_interproscan/{cluster}.{num}.subtree.interproscan.log"
    run:
        try:
            shell("grep -c '>' {input}")
        except:
            shell("touch {output}")
        else:
            shell("sed -i 's/*//g' {input}")
            shell("interproscan.sh -appl Pfam -cpu {threads} -f tsv -iprlookup --pathways -i {input} -o {output}")

rule colour_subtree:
    input:
        map = config["mapping_file"],
        info = config["dataset_chacteristics"],
        phyla = "{workdir}/7.subtree/all_subset_phyla.list",
        domains = "{workdir}/8.subtree_interproscan/{cluster}.{num}.subtree.interproscan.tsv",
        tree = "{workdir}/7.subtree/{cluster}.{num}.subtree.align.trim.treefile",
        faa = "{workdir}/6.subtree_faa/{cluster}.{num}.subtree.faa"
    output:
        out = "{workdir}/7.subtree/{cluster}.{num}.subtree.align.trim.treefile.pdf",
    params:
        cluster_name = "{cluster}",
        focal_taxonomy = config["focal_taxonomy"],
        num = "{num}"
    script:
        "scripts/colour_subtree.py"

rule HGT_sources:
    input:
        map = config["mapping_file"],
        tree = "{workdir}/7.subtree/{cluster}.{num}.subtree.align.trim.treefile"
    output:
        out = "{workdir}/7.subtree/{cluster}.{num}.sister_nested_taxonomy.tsv"
    params:
        cluster_name = "{cluster}",
        focal_taxonomy = config["focal_taxonomy"],
        num = "{num}"
    script:
        "scripts/find_sister_nested_taxonomy.py"

rule infer_HGT_sources:
    input:
        expand("{workdir}/7.subtree/{cluster}.{num}.sister_nested_taxonomy.tsv", cluster=config["clusters"], workdir=config["workdir"], num=num)
    output:
        "{workdir}/ALL.sister_nested_taxonomy.tsv"
    run:
        shell("cat {input} > {output}")

rule plot_HGT_sources:
    input:
        "{workdir}/ALL.sister_nested_taxonomy.tsv"
    output:
        "{workdir}/ALL.sister_nested_taxonomy.pdf"
    script:
        "scripts/plot_HGT_sources.R"
