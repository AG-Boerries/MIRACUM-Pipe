GEMINI 0.20.1 databases
=======================

clinvar 30.1.2017   bash pipeline uses 5.9.2017

cosmic  v68         bash pipeline uses v84?

dbsnp   b147

ExAC    r0.3

gnomad_exome  r2.0.1    bash pipeline version ?

ESP6500si   v2

1000g   phase3

fitcons v1.01

cadd13  as add-on

vista enhancers db  8.11.2013

HPRD interactions   ??


Exhaustive lists of fields stored in a GEMINI 0.20.1 sqlite database
====================================================================

::

    table_name          column_name                   type      
    variants            chrom                         VARCHAR(20)
    variants            start                         INTEGER   
    variants            end                           INTEGER   
    variants            vcf_id                        TEXT      
    variants            variant_id                    INTEGER   
    variants            anno_id                       INTEGER   
    variants            ref                           TEXT      
    variants            alt                           TEXT      
    variants            qual                          FLOAT     
    variants            filter                        TEXT      
    variants            type                          VARCHAR(20)
    variants            sub_type                      TEXT      
    variants            gts                           BLOB      
    variants            gt_types                      BLOB      
    variants            gt_phases                     BLOB      
    variants            gt_depths                     BLOB      
    variants            gt_ref_depths                 BLOB      
    variants            gt_alt_depths                 BLOB      
    variants            gt_alt_freqs                  BLOB      
    variants            gt_quals                      BLOB      
    variants            gt_copy_numbers               BLOB      
    variants            gt_phred_ll_homref            BLOB      
    variants            gt_phred_ll_het               BLOB      
    variants            gt_phred_ll_homalt            BLOB      
    variants            call_rate                     FLOAT     
    variants            max_aaf_all                   FLOAT     
    variants            in_dbsnp                      BOOLEAN   
    variants            rs_ids                        TEXT      
    variants            sv_cipos_start_left           INTEGER   
    variants            sv_cipos_end_left             INTEGER   
    variants            sv_cipos_start_right          INTEGER   
    variants            sv_cipos_end_right            INTEGER   
    variants            sv_length                     INTEGER   
    variants            sv_is_precise                 BOOLEAN   
    variants            sv_tool                       TEXT      
    variants            sv_evidence_type              TEXT      
    variants            sv_event_id                   TEXT      
    variants            sv_mate_id                    TEXT      
    variants            sv_strand                     TEXT      
    variants            in_omim                       BOOLEAN   
    variants            clinvar_sig                   TEXT      
    variants            clinvar_disease_name          TEXT      
    variants            clinvar_dbsource              TEXT      
    variants            clinvar_dbsource_id           TEXT      
    variants            clinvar_origin                TEXT      
    variants            clinvar_dsdb                  TEXT      
    variants            clinvar_dsdbid                TEXT      
    variants            clinvar_disease_acc           TEXT      
    variants            clinvar_in_locus_spec_db      BOOLEAN   
    variants            clinvar_on_diag_assay         BOOLEAN   
    variants            clinvar_causal_allele         TEXT      
    variants            clinvar_gene_phenotype        TEXT      
    variants            geno2mp_hpo_ct                INTEGER   
    variants            pfam_domain                   TEXT      
    variants            cyto_band                     TEXT      
    variants            rmsk                          TEXT      
    variants            in_cpg_island                 BOOLEAN   
    variants            in_segdup                     BOOLEAN   
    variants            is_conserved                  BOOLEAN   
    variants            gerp_bp_score                 FLOAT     
    variants            gerp_element_pval             FLOAT     
    variants            num_hom_ref                   INTEGER   
    variants            num_het                       INTEGER   
    variants            num_hom_alt                   INTEGER   
    variants            num_unknown                   INTEGER   
    variants            aaf                           FLOAT     
    variants            hwe                           FLOAT     
    variants            inbreeding_coeff              FLOAT     
    variants            pi                            FLOAT     
    variants            recomb_rate                   FLOAT     
    variants            gene                          VARCHAR(60)
    variants            transcript                    VARCHAR(60)
    variants            is_exonic                     BOOLEAN   
    variants            is_coding                     BOOLEAN   
    variants            is_splicing                   BOOLEAN   
    variants            is_lof                        BOOLEAN   
    variants            exon                          TEXT      
    variants            codon_change                  TEXT      
    variants            aa_change                     TEXT      
    variants            aa_length                     TEXT      
    variants            biotype                       TEXT      
    variants            impact                        VARCHAR(60)
    variants            impact_so                     TEXT      
    variants            impact_severity               VARCHAR(20)
    variants            polyphen_pred                 TEXT      
    variants            polyphen_score                FLOAT     
    variants            sift_pred                     TEXT      
    variants            sift_score                    FLOAT     
    variants            anc_allele                    TEXT      
    variants            rms_bq                        FLOAT     
    variants            cigar                         TEXT      
    variants            depth                         INTEGER   
    variants            strand_bias                   FLOAT     
    variants            rms_map_qual                  FLOAT     
    variants            in_hom_run                    INTEGER   
    variants            num_mapq_zero                 INTEGER   
    variants            num_alleles                   INTEGER   
    variants            num_reads_w_dels              FLOAT     
    variants            haplotype_score               FLOAT     
    variants            qual_depth                    FLOAT     
    variants            allele_count                  INTEGER   
    variants            allele_bal                    FLOAT     
    variants            in_hm2                        BOOLEAN   
    variants            in_hm3                        BOOLEAN   
    variants            is_somatic                    BOOLEAN   
    variants            somatic_score                 FLOAT     
    variants            in_esp                        BOOLEAN   
    variants            aaf_esp_ea                    FLOAT     
    variants            aaf_esp_aa                    FLOAT     
    variants            aaf_esp_all                   FLOAT     
    variants            exome_chip                    BOOLEAN   
    variants            in_1kg                        BOOLEAN   
    variants            aaf_1kg_amr                   FLOAT     
    variants            aaf_1kg_eas                   FLOAT     
    variants            aaf_1kg_sas                   FLOAT     
    variants            aaf_1kg_afr                   FLOAT     
    variants            aaf_1kg_eur                   FLOAT     
    variants            aaf_1kg_all                   FLOAT     
    variants            grc                           TEXT      
    variants            gms_illumina                  FLOAT     
    variants            gms_solid                     FLOAT     
    variants            gms_iontorrent                FLOAT     
    variants            in_cse                        BOOLEAN   
    variants            encode_tfbs                   TEXT      
    variants            encode_dnaseI_cell_count      INTEGER   
    variants            encode_dnaseI_cell_list       TEXT      
    variants            encode_consensus_gm12878      TEXT      
    variants            encode_consensus_h1hesc       TEXT      
    variants            encode_consensus_helas3       TEXT      
    variants            encode_consensus_hepg2        TEXT      
    variants            encode_consensus_huvec        TEXT      
    variants            encode_consensus_k562         TEXT      
    variants            vista_enhancers               TEXT      
    variants            cosmic_ids                    TEXT      
    variants            info                          BLOB      
    variants            cadd_raw                      FLOAT     
    variants            cadd_scaled                   FLOAT     
    variants            fitcons                       FLOAT     
    variants            in_exac                       BOOLEAN   
    variants            aaf_exac_all                  FLOAT     
    variants            aaf_adj_exac_all              FLOAT     
    variants            aaf_adj_exac_afr              FLOAT     
    variants            aaf_adj_exac_amr              FLOAT     
    variants            aaf_adj_exac_eas              FLOAT     
    variants            aaf_adj_exac_fin              FLOAT     
    variants            aaf_adj_exac_nfe              FLOAT     
    variants            aaf_adj_exac_oth              FLOAT     
    variants            aaf_adj_exac_sas              FLOAT     
    variants            exac_num_het                  INTEGER   
    variants            exac_num_hom_alt              INTEGER   
    variants            exac_num_chroms               INTEGER   
    variants            aaf_gnomad_all                FLOAT     
    variants            aaf_gnomad_afr                FLOAT     
    variants            aaf_gnomad_amr                FLOAT     
    variants            aaf_gnomad_asj                FLOAT     
    variants            aaf_gnomad_eas                FLOAT     
    variants            aaf_gnomad_fin                FLOAT     
    variants            aaf_gnomad_nfe                FLOAT     
    variants            aaf_gnomad_oth                FLOAT     
    variants            aaf_gnomad_sas                FLOAT     
    variants            gnomad_num_het                INTEGER   
    variants            gnomad_num_hom_alt            INTEGER   
    variants            gnomad_num_chroms             INTEGER   
    variant_impacts     variant_id                    INTEGER   
    variant_impacts     anno_id                       INTEGER   
    variant_impacts     gene                          VARCHAR(60)
    variant_impacts     transcript                    VARCHAR(60)
    variant_impacts     is_exonic                     BOOLEAN   
    variant_impacts     is_coding                     BOOLEAN   
    variant_impacts     is_lof                        BOOLEAN   
    variant_impacts     exon                          TEXT      
    variant_impacts     codon_change                  TEXT      
    variant_impacts     aa_change                     TEXT      
    variant_impacts     aa_length                     TEXT      
    variant_impacts     biotype                       TEXT      
    variant_impacts     impact                        VARCHAR(60)
    variant_impacts     impact_so                     TEXT      
    variant_impacts     impact_severity               VARCHAR(20)
    variant_impacts     polyphen_pred                 TEXT      
    variant_impacts     polyphen_score                FLOAT     
    variant_impacts     sift_pred                     TEXT      
    variant_impacts     sift_score                    FLOAT     
    samples             sample_id                     INTEGER   
    samples             family_id                     TEXT      
    samples             name                          VARCHAR(50)
    samples             paternal_id                   TEXT      
    samples             maternal_id                   TEXT      
    samples             sex                           TEXT      
    samples             phenotype                     TEXT      
    gene_detailed       uid                           INTEGER   
    gene_detailed       chrom                         VARCHAR(60)
    gene_detailed       gene                          VARCHAR(60)
    gene_detailed       is_hgnc                       BOOLEAN   
    gene_detailed       ensembl_gene_id               TEXT      
    gene_detailed       transcript                    VARCHAR(60)
    gene_detailed       biotype                       TEXT      
    gene_detailed       transcript_status             TEXT      
    gene_detailed       ccds_id                       VARCHAR(60)
    gene_detailed       hgnc_id                       TEXT      
    gene_detailed       entrez_id                     TEXT      
    gene_detailed       cds_length                    TEXT      
    gene_detailed       protein_length                TEXT      
    gene_detailed       transcript_start              TEXT      
    gene_detailed       transcript_end                TEXT      
    gene_detailed       strand                        TEXT      
    gene_detailed       synonym                       TEXT      
    gene_detailed       rvis_pct                      FLOAT     
    gene_detailed       mam_phenotype_id              TEXT      
    gene_summary        uid                           INTEGER   
    gene_summary        chrom                         VARCHAR(60)
    gene_summary        gene                          VARCHAR(60)
    gene_summary        is_hgnc                       BOOLEAN   
    gene_summary        ensembl_gene_id               TEXT      
    gene_summary        hgnc_id                       TEXT      
    gene_summary        transcript_min_start          INTEGER   
    gene_summary        transcript_max_end            INTEGER   
    gene_summary        strand                        TEXT      
    gene_summary        synonym                       TEXT      
    gene_summary        rvis_pct                      FLOAT     
    gene_summary        mam_phenotype_id              TEXT      
    gene_summary        in_cosmic_census              BOOLEAN

