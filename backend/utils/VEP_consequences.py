# Mapping of VEP consequence terms to their attributes
VEP_CONSEQUENCES = [
    {
        "term": "transcript_ablation",
        "rank": 1,
        "impact": "HIGH",
        "description": "A feature ablation whereby the deleted region includes a transcript feature",
        "display_term": "Transcript ablation"
    },
    {
        "term": "splice_acceptor_variant",
        "rank": 2,
        "impact": "HIGH",
        "description": "A splice variant that changes the 2 base region at the 3' end of an intron",
        "display_term": "Splice acceptor variant"
    },
    {
        "term": "splice_donor_variant",
        "rank": 3,
        "impact": "HIGH",
        "description": "A splice variant that changes the 2 base region at the 5' end of an intron",
        "display_term": "Splice donor variant"
    },
    {
        "term": "stop_gained",
        "rank": 4,
        "impact": "HIGH",
        "description": "A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript",
        "display_term": "Stop gained"
    },
    {
        "term": "frameshift_variant",
        "rank": 5,
        "impact": "HIGH",
        "description": "A sequence variant which causes a disruption of the translational reading frame",
        "display_term": "Frameshift variant"
    },
    {
        "term": "start_lost",
        "rank": 6,
        "impact": "HIGH",
        "description": "A codon variant that changes at least one base of the start codon",
        "display_term": "Start lost"
    },
    {
        "term": "transcript_amplification",
        "rank": 7,
        "impact": "HIGH",
        "description": "A feature amplification of a region containing a transcript",
        "display_term": "Transcript amplification"
    },
    {
        "term": "inframe_insertion",
        "rank": 8,
        "impact": "MODERATE",
        "description": "An inframe insertion of one or more bases into the coding sequence",
        "display_term": "Inframe insertion"
    },
    {
        "term": "inframe_deletion",
        "rank": 9,
        "impact": "MODERATE",
        "description": "An inframe deletion of one or more bases from the coding sequence",
        "display_term": "Inframe deletion"
    },
    {
        "term": "missense_variant",
        "rank": 10,
        "impact": "MODERATE",
        "description": "A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved",
        "display_term": "Missense variant"
    },
    {
        "term": "protein_altering_variant",
        "rank": 11,
        "impact": "MODERATE",
        "description": "A sequence variant which is predicted to change the protein encoded in the coding sequence",
        "display_term": "Protein altering variant"
    },
    {
        "term": "splice_region_variant",
        "rank": 12,
        "impact": "LOW",
        "description": "A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron",
        "display_term": "Splice region variant"
    },
    {
        "term": "incomplete_terminal_codon_variant",
        "rank": 13,
        "impact": "LOW",
        "description": "A sequence variant where at least one base of the final codon of an incompletely annotated transcript is changed",
        "display_term": "Incomplete terminal codon variant"
    },
    {
        "term": "start_retained_variant",
        "rank": 14,
        "impact": "LOW",
        "description": "A sequence variant where at least one base in the start codon is changed, but the start remains",
        "display_term": "Start retained variant"
    },
    {
        "term": "stop_retained_variant",
        "rank": 15,
        "impact": "LOW",
        "description": "A sequence variant where at least one base in the stop codon is changed, but the stop remains",
        "display_term": "Stop retained variant"
    },
    {
        "term": "synonymous_variant",
        "rank": 16,
        "impact": "LOW",
        "description": "A sequence variant where there is no resulting change to the encoded amino acid",
        "display_term": "Synonymous variant"
    },
    {
        "term": "coding_sequence_variant",
        "rank": 17,
        "impact": "MODIFIER",
        "description": "A sequence variant that changes the coding sequence",
        "display_term": "Coding sequence variant"
    },
    {
        "term": "mature_miRNA_variant",
        "rank": 18,
        "impact": "MODIFIER",
        "description": "A sequence variant located within the mature miRNA region",
        "display_term": "Mature miRNA variant"
    },
    {
        "term": "5_prime_UTR_variant",
        "rank": 19,
        "impact": "MODIFIER",
        "description": "A UTR variant of the 5' UTR",
        "display_term": "5 prime UTR variant"
    },
    {
        "term": "3_prime_UTR_variant",
        "rank": 20,
        "impact": "MODIFIER",
        "description": "A UTR variant of the 3' UTR",
        "display_term": "3 prime UTR variant"
    },
    {
        "term": "non_coding_transcript_exon_variant",
        "rank": 21,
        "impact": "MODIFIER",
        "description": "A sequence variant that changes non-coding exon sequence in a non-coding transcript",
        "display_term": "Non coding transcript exon variant"
    },
    {
        "term": "intron_variant",
        "rank": 22,
        "impact": "MODIFIER",
        "description": "A transcript variant occurring within an intron",
        "display_term": "Intron variant"
    },
    {
        "term": "NMD_transcript_variant",
        "rank": 23,
        "impact": "MODIFIER",
        "description": "A variant in a transcript that is the target of NMD",
        "display_term": "NMD transcript variant"
    },
    {
        "term": "non_coding_transcript_variant",
        "rank": 24,
        "impact": "MODIFIER",
        "description": "A transcript variant of a non coding transcript",
        "display_term": "Non coding transcript variant"
    },
    {
        "term": "upstream_gene_variant",
        "rank": 25,
        "impact": "MODIFIER",
        "description": "A sequence variant located upstream of a gene",
        "display_term": "Upstream gene variant"
    },
    {
        "term": "downstream_gene_variant",
        "rank": 26,
        "impact": "MODIFIER",
        "description": "A sequence variant located downstream of a gene",
        "display_term": "Downstream gene variant"
    },
    {
        "term": "TFBS_ablation",
        "rank": 27,
        "impact": "MODIFIER",
        "description": "A feature ablation whereby the deleted region includes a transcription factor binding site",
        "display_term": "TFBS ablation"
    },
    {
        "term": "TFBS_amplification",
        "rank": 28,
        "impact": "MODIFIER",
        "description": "A feature amplification of a region containing a transcription factor binding site",
        "display_term": "TFBS amplification"
    },
    {
        "term": "TF_binding_site_variant",
        "rank": 29,
        "impact": "MODIFIER",
        "description": "A sequence variant located within a transcription factor binding site",
        "display_term": "TF binding site variant"
    },
    {
        "term": "regulatory_region_ablation",
        "rank": 30,
        "impact": "MODIFIER",
        "description": "A feature ablation whereby the deleted region includes a regulatory region",
        "display_term": "Regulatory region ablation"
    },
    {
        "term": "regulatory_region_amplification",
        "rank": 31,
        "impact": "MODIFIER",
        "description": "A feature amplification of a region containing a regulatory region",
        "display_term": "Regulatory region amplification"
    },
    {
        "term": "feature_elongation",
        "rank": 32,
        "impact": "MODIFIER",
        "description": "A sequence variant that causes the extension of a genomic feature, with regard to the reference sequence",
        "display_term": "Feature elongation"
    },
    {
        "term": "regulatory_region_variant",
        "rank": 33,
        "impact": "MODIFIER",
        "description": "A sequence variant located within a regulatory region",
        "display_term": "Regulatory region variant"
    },
    {
        "term": "feature_truncation",
        "rank": 34,
        "impact": "MODIFIER",
        "description": "A sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence",
        "display_term": "Feature truncation"
    },
    {
        "term": "intergenic_variant",
        "rank": 35,
        "impact": "MODIFIER",
        "description": "A sequence variant located in the intergenic region, between genes",
        "display_term": "Intergenic variant"
    }
]

VEP_RANK_DICT = {entry["term"]: entry for entry in VEP_CONSEQUENCES}
