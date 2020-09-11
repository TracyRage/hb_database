library('Gmisc')
library('glue')
library('dplyr')

flow_chart <- function(enzyme_input, gene_input, accession_input,
         pmid_input, organism_input,
         cdd_input, ko_input, pfam_input, 
         outgroup_input, reviewed_input,
         unfiltered_input, filtered_input, 
         total_seqs_input) {
  seed <- boxGrob(glue("Curated seed (NCBI):",
                       "Enzyme: {enzyme}",
                       "Organism: {organism}",
                       "Gene: {gene}",
                       "Accession: {accession}",
                       "PMID: {pmid}",
                       enzyme = enzyme_input,
                       gene = gene_input,
                       accession = accession_input,
                       pmid = pmid_input,
                       organism = organism_input,
                       .sep = "\n"))
  
  seed_domain <- boxGrob(glue("CDD: {cdd}",
                              "Ko: {ko}",
                              "Pfam: {pfam}",
                              cdd = cdd_input,
                              ko = ko_input,
                              pfam = pfam_input,
                              .sep="\n"))
  
  blast_result <- boxGrob(glue("BLASTP results:",
                               "Database: UniProt",
                               "Outgroup: {outgroup}",
                               "Reviewed: {reviewed}",
                               outgroup = outgroup_input,
                               reviewed = reviewed_input,
                               .sep = "\n"))
  
  sp_result <- boxGrob(glue("UniProt query results:",
                            "Unfiltered: {unfiltered}",
                            "- fragments",
                            "- uncultured",
                            "Filtered: {filtered}",
                            unfiltered = unfiltered_input,
                            filtered = filtered_input,
                            .sep="\n"))
  
  tree <- boxGrob(glue("Tree features:",
                       "Total sequences: {total_seqs}",
                       "Clade exists if:",
                       "UFboot >= 95 & SH-aLRT >= 80",
                       total_seqs = total_seqs_input,
                       .sep="\n"))
  
  vert <- spreadVertical(seed = seed, 
                         seed_domain = seed_domain,
                         blast_res = blast_result,
                         phylo_tree = tree)
  treat <- alignVertical(vert$seed_domain, blast_result,
                         sp_result) %>% 
    spreadHorizontal()
  vert$blast_res <- NULL
  
  result_1 <- connectGrob(vert$seed, treat[[1]], type = "Z")
  result_2 <- connectGrob(vert$seed, vert$seed_domain, type="N")
  result_3 <- connectGrob(treat[[1]], vert$phylo_tree, type = "L")
  result_4 <- connectGrob(treat[[2]], vert$phylo_tree, type = "L")
  result_5 <- vert
  result_6 <- treat
  
  general_result <- list(result_1, result_2, result_3, 
                         result_4, result_5, result_6)
  
  return(general_result)
  
}


# Biphenyl_2,3_dioxygenase ------------------------------------------------

# flow_chart(enzyme_input = "Biphenyl 2,3 dioxygenase",
#            gene_input = "bphA1",
#            accession_input = "ABE37059",
#            organism_input = "Paraburkholderia xenovorans",
#            pmid_input = "17030797",
#            cdd_input = "239550",
#            ko_input = "K08689",
#            pfam_input = "PF00355.27",
#            outgroup_input = 35,
#            reviewed_input = 2,
#            unfiltered_input = 100,
#            filtered_input = 66,
#            total_seqs_input = 37+66)


# Dibenzofuran_dioxygenase ------------------------------------------------

# flow_chart(enzyme_input = "Dibenzofuran dioxygenase",
#            gene_input = 'dbfA1',
#            accession_input = 'ALS21084',
#            organism_input = 'Paenibacillus naphthalenovorans',
#            pmid_input = '26868401',
#            cdd_input = '239550',
#            ko_input = 'K14599',
#            pfam_input = 'PF00355.27',
#            outgroup_input = 39,
#            reviewed_input = 0,
#            unfiltered_input = 10,
#            filtered_input = 5,
#            total_seqs_input = 5+39
#            )


# Diphenyl_ether_dioxygenase ----------------------------------------------

# flow_chart(enzyme_input = 'Diphenylether dioxygenase',
#            gene_input = "dpeA1",
#            accession_input = "ANW37879.1",
#            organism_input = "Sphingobium phenoxybenzoativorans",
#            pmid_input = "Direct submission",
#            cdd_input = "239550",
#            ko_input = "N/A",
#            pfam_input = "PF00355.27",
#            outgroup_input = 38,
#            reviewed_input = 0,
#            unfiltered_input = 1,
#            filtered_input = 1,
#            total_seqs_input = 39)


# Ethylbenzene dioxygenase ----------------------------------------------

# flow_chart(enzyme_input = "Ethylbenzene dioxygenase",
#            gene_input = "ebdA1",
#            accession_input = "BAC92718.1",
#            organism_input = "Rhodococcus jostii",
#            pmid_input = "Direct Submission",
#            cdd_input = "239550",
#            ko_input = "K14748",
#            pfam_input = "PF00355.27",
#            outgroup_input = 38,
#            reviewed_input = 0,
#            unfiltered_input = 1,
#            filtered_input = 1,
#            total_seqs_input = 39)


# Ethylbenzene dehydrogenase ----------------------------------------------

# flow_chart(enzyme_input = "Ethylbenzene dehydrogenase",
#            gene_input = "ebdA1",
#            accession_input = "CAI07432",
#            organism_input = "Aromatoleum aromaticum",
#            pmid_input = "15551059",
#            cdd_input = "238218",
#            ko_input = "K00370",
#            pfam_input = "PF00384.23",
#            outgroup_input = 50,
#            reviewed_input = 0,
#            unfiltered_input = 3,
#            filtered_input = 3,
#            total_seqs_input = 53
#            )


# Naphthalene dioxygenase (nahAc) -----------------------------------------

# flow_chart(enzyme_input = "Naphthalene dioxygenase (nahAc)",
#            gene_input = "nahAc",
#            accession_input = "AAP44288",
#            organism_input = "Pseudomonas putida",
#            pmid_input = "15246534",
#            cdd_input = "239550",
#            ko_input = "K14579",
#            pfam_input = "PF00355.27",
#            outgroup_input = 38,
#            reviewed_input = 0,
#            unfiltered_input = 100,
#            filtered_input = 6,
#            total_seqs_input = 44)


# Naphthalene dioxygenase (narAa) -----------------------------------------

# flow_chart(enzyme_input = "Naphthalene dioxygenase (narAa)",
#            gene_input = "narAa",
#            accession_input = "ADM94827.1",
#            organism_input = "Rhodococcus",
#            pmid_input = "21503712",
#            cdd_input = "239550",
#            ko_input = "N/A",
#            pfam_input = "PF00355.27",
#            outgroup_input = 38,
#            reviewed_input = 0,
#            unfiltered_input = 19,
#            filtered_input = 12,
#            total_seqs_input = 50
#            )


# PAH dioxygenase (a) -----------------------------------------------------

# flow_chart(enzyme_input = "PAH dioxygenase (a)",
#            gene_input = "ahdA1a",
#            accession_input = "BAC65448.1",
#            organism_input = "Sphingomonas",
#            pmid_input = "12565867",
#            cdd_input = "239550",
#            ko_input = "K14748",
#            pfam_input = "PF00355.27",
#            outgroup_input = 37,
#            reviewed_input = 0,
#            unfiltered_input = 1,
#            filtered_input = 1,
#            total_seqs_input = 38
#            )


# PAH dioxygenase (b) -----------------------------------------------------

# flow_chart(enzyme_input = "PAH dioxygenase (b)",
#            gene_input = "ahdA1b",
#            accession_input = "BAC65446.1",
#            organism_input = "Sphingomonas",
#            pmid_input = "12565867",
#            cdd_input = "239550",
#            ko_input = "K14748",
#            pfam_input = "PF00355.27",
#            outgroup_input = 37,
#            reviewed_input = 0,
#            unfiltered_input = 2,
#            filtered_input = 2,
#            total_seqs_input = 39)


# PAH dioxygenase (c) -----------------------------------------------------

# flow_chart(enzyme_input = "PAH dioxygenase (c)",
#            gene_input = "ahdA1c",
#            accession_input = "AGZ63455.1",
#            organism_input = "Sphingobium",
#            pmid_input = "Direct Submission",
#            cdd_input = "239550",
#            ko_input = "K16319",
#            pfam_input = "PF00355.27",
#            outgroup_input = 38,
#            reviewed_input = 0,
#            unfiltered_input = 4,
#            filtered_input = 4,
#            total_seqs_input = 42)


# PAH dioxygenase (d) -----------------------------------------------------

# flow_chart(enzyme_input = "PAH dioxygenase (d)",
#            gene_input = "ahdA1d",
#            accession_input = "AGZ63462.1",
#            organism_input = "Sphingobium",
#            pmid_input = "Direct Submission",
#            cdd_input = "239550",
#            ko_input = "K16319",
#            pfam_input = "PF00355.27",
#            outgroup_input = 37,
#            reviewed_input = 0,
#            unfiltered_input = 2,
#            filtered_input = 2,
#            total_seqs_input = 39
#            )


# PAH dioxygenase (f) -----------------------------------------------------

# flow_chart(enzyme_input = "PAH dioxygenase (f)",
#            gene_input = "ahdA1f",
#            accession_input = "AGZ63449.1",
#            organism_input = "Sphingobium",
#            pmid_input = "Direct Submission",
#            cdd_input = "239550",
#            ko_input = "N/A",
#            pfam_input = "PF00355.27",
#            outgroup_input = 38,
#            reviewed_input = 0,
#            unfiltered_input = 1,
#            filtered_input = 1,
#            total_seqs_input = 39
#            )


# P-cymene-methyl hydroxylase ---------------------------------------------

# flow_chart(enzyme_input = "P-cymene-methyl hydroxylase",
#            gene_input = "cymA",
#            accession_input = "AAC45296.1",
#            organism_input = "Pseudomonas chlororaphis",
#            pmid_input = "9144566",
#            cdd_input = "238126",
#            ko_input = "N/A",
#            pfam_input = "PF00111.28",
#            outgroup_input = 50,
#            reviewed_input = 0,
#            unfiltered_input = 5,
#            filtered_input = 5,
#            total_seqs_input = 55
#            )


# Phenol 2-monooxygenase --------------------------------------------------

# flow_chart(enzyme_input = "Phenol 2-monooxygenase",
#            gene_input = "pheA1",
#            accession_input = "ABS30825.1",
#            organism_input = "Rhodococcus erythropolis",
#            pmid_input = "19787347",
#            cdd_input = "N/A",
#            ko_input = "K00483",
#            pfam_input = "PF11794.9",
#            outgroup_input = 12,
#            reviewed_input = 0,
#            unfiltered_input = 54,
#            filtered_input = 34,
#            total_seqs_input = 12+34)


# Styrene monooxygenase ---------------------------------------------------

# flow_chart(enzyme_input = "Styrene monooxygenase",
#            gene_input = "styA",
#            accession_input = "ABB03727.1",
#            organism_input = "Pseudomonas putida",
#            pmid_input = "Direct submission",
#            cdd_input = "N/A",
#            ko_input = "N/A",
#            pfam_input = "PF17885.2",
#            outgroup_input = 0,
#            reviewed_input = 2,
#            unfiltered_input = 85,
#            filtered_input = 82,
#            total_seqs_input = 84
#            )


# Toluene 2-monooxygenase (tbc) -------------------------------------------

# flow_chart(enzyme_input = "Toluene 2-monooxygenase",
#            gene_input = "tbc2A",
#            accession_input = "AAG40794.1",
#            organism_input = "Burkholderia cepacia",
#            pmid_input = " 11571188",
#            cdd_input = "N/A",
#            ko_input = "K15760",
#            pfam_input = "PF02332.19",
#            outgroup_input = 9,
#            reviewed_input = 0,
#            unfiltered_input = 1,
#            filtered_input = 1,
#            total_seqs_input = 10)


# Toluene 2-monooxygenase -------------------------------------------------

# flow_chart(enzyme_input = "Toluene 2-monooxygenase",
#            gene_input = "tomA1",
#            accession_input = "AAK07411.1",
#            organism_input = "Burkholderia cepacia",
#            pmid_input = "Direct Submission",
#            cdd_input = "153097",
#            ko_input = "K16242",
#            pfam_input = "PF02332.19",
#            outgroup_input = 9,
#            reviewed_input = 0,
#            unfiltered_input = 3,
#            filtered_input = 3,
#            total_seqs_input = 12)


# Toluene 3-monooxygenase -------------------------------------------------
# 
# flow_chart(enzyme_input = "Toluene 3-monooxygenase",
#            gene_input = "tbuA1",
#            accession_input = "AAB09618.1",
#            organism_input = "Ralstonia pickettii",
#            pmid_input = "7867951",
#            cdd_input = "N/A",
#            ko_input = "K15760",
#            pfam_input = "PF02332.19",
#            outgroup_input = 9,
#            reviewed_input = 0,
#            unfiltered_input = 3,
#            filtered_input = 3,
#            total_seqs_input = 12)


# Toluene 4-monooxygenase -------------------------------------------------

# flow_chart(enzyme_input = "Toluene 4-monooxygenase",
#            gene_input = "tmoA1",
#            accession_input = "QCT24447.1",
#            organism_input = "Pseudomonas oleovorans",
#            pmid_input = "Direct submission",
#            cdd_input = 'N/A',
#            ko_input = "K15760",
#            pfam_input = "PF02332.19",
#            outgroup_input = 8,
#            reviewed_input = 1,
#            unfiltered_input = 91,
#            filtered_input = 54,
#            total_seqs_input = 9+54)


# Toluene-benzene 2-monooxygenase -----------------------------------------

# flow_chart(enzyme_input = "Toluene/benzene 2-monooxygenase",
#            gene_input = "touA",
#            accession_input = "AAT40431.1",
#            organism_input = "Pseudomonas",
#            pmid_input = "15184119",
#            cdd_input = "153097",
#            ko_input = "K15760",
#            pfam_input = "PF02332.19",
#            outgroup_input = 9,
#            reviewed_input = 0,
#            unfiltered_input = 2,
#            filtered_input = 2,
#            total_seqs_input = 11)


# Toluene dioxygenase -----------------------------------------------------

# flow_chart(enzyme_input = "Toluene dioxygenase",
#            gene_input = "todC1",
#            accession_input = "BAC05504",
#            organism_input = "Thauera",
#            pmid_input = '15006757',
#            cdd_input = "N/A",
#            ko_input = "N/A",
#            pfam_input = "PF00848.20",
#            outgroup_input = 26,
#            reviewed_input = 0,
#            unfiltered_input = 15,
#            filtered_input = 2,
#            total_seqs_input = 26+2)
