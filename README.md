### Init conda environment


1. Download and install [Anaconda](https://www.anaconda.com/products/individual), set [Bioconda](https://bioconda.github.io/user/install.html) as well.

2. From `~/hb_database/hb_database_env.txt` create the working environment.

		conda env create -file hb_database_env.txt

3. Init conda environment:

		conda activate hb_database

### Install local databases

#### Create local SwissProt database
Download SwissProt [uniprot_sprot.fa](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz) and
create a local database.

		makeblastdb  \
		-in ~/db_uniprot/uniprot_sprot.fasta \
		-out ~/db_uniprot/ \
		-parse_seqids \
		-dbtype prot

#### Create local Cdd database for RPS-Blast
[Download](ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/) Cdd database.

		makeprofiledb \
		-title Cdd \
		-in Cdd.pn \
		-out db_cdd \
		-threshold 9.82 \
		-scale 100.0 \
		-dbtype rps \
		-index true


#### Create Pfam local database
[Download](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.1/) Pfam-A.hmm & Pfam-A.hmm.dat and create a local Pfam database.

		hmmpress -f data/seqs/db_pfam/Pfam-A.hmm

#### Create local ghostKoala database
Check [here](https://taylorreiter.github.io/2019-05-11-kofamscan/) for set up steps.

----------------------------------------------------------------

### Create database

#### Get sequences of curated alpha subunits
Get from NCBI the `.fasta` file for the curated alpha subunits.
> Only one alpha subunit per enzyme of interest.

		python 01_fetching_aa_seqs.py  \
		--input data/seqs/01_list_alpha_subunits.csv \
		--email <insert your email>

> Input table example:

enzyme | alpha_subunits | accession
------ | -------------- | ---------
Xylene_monooxygenase | xylM | NP_542887
Toluene_dioxygenase |todC1 | BAC05504
Toluene_2-monooxygenase | tomA1 | AAK07411.1

> Do not use symbols such as (",", "/", etc.) in your input table. It may cause various parsing failures.

#### Get Blast results for each gene sequence
`02_alpha_subunits_aa_sequences` splits your multiple sequence `02_alpha_subunits_aa_sequences.fa` file, and BLASTP each singular
`fasta` file.

		python 02_do_blast.py \
		--input analyses/02_alpha_subunits_aa_sequences.fa  \
		--database data/seqs/db_uniprot/swiss_prot.fa

> Check in the output xml files last/trailing entries; sometimes you can get eukaryotic genes

#### Get Blast hits (interest, outgroup, total)
Parse `.xml` files and get UniProt reviewed and unreviewed (outgroup) sequences.

		python 03_get_blast_hits.py \
		--xml_in analyses/blast_outputs/separated_blast/

> Nota Bene:

* `separated_hits/interest` (SwissProt reviewed sequences)
* `separated_hits/outgroup` (outgroup sequences)
* `separated_hits/merged` (reviewed + outgroup)

#### Get fasta file for each Blast hit
Get `.fasta` file for each hit.

		python 04_get_uniprot_fa.py \
		--input analyses/blast_outputs/separated_hits/

#### Get RPS-Blast (Cdd, Ko, Pfam) hits
`05_do_rps_blast.py` uses `.fa` files in `blast_outputs/separated_hits` and outputs `.cdd`, `.ko` and `.pfam` tables.

		python 05_do_rps_blast.py \
		--input_fa analyses/blast_outputs/separated_hits/ \
		--db_cdd data/seqs/db_cdd/db_cdd \
		--evalue 0.01 \
		--db_pfam data/seqs/db_pfam/ \
		--db_ko data/seqs/db_ko/profiles/ \
		--ko_list data/seqs/db_ko/ko_list \
		--ko_config data/seqs/db_ko/config.yml

> To enhance the performance of Ko annotation check `.yml` file.

#### Generate annotation tables for each alpha subunit
Merge `.cdd`, `.pfam`, `.ko` tables in `blast_outputs/separated_hits`.

		python 06_create_annotated_table.py /
		--fasta_input analyses/blast_outputs/separated_hits/

> Output example (ABE37059.1_Paraburkholderia_xenovorans_LB400.tsv):

UniProt | Cdd | Ko | Pfam
------- | --- | -- | ----
Q52438  | CDD:239550 | K08689 | PF00355.27
Q53122  | CDD:239550 | K08689 | PF00355.27

#### Access UniProt API and get sequences based on keyword and perform filtering
Perform UniProt direct query.

> Keyword implies exact gene name. Filtering implies deleting uncultured and fragments sequences.

		python 07_do_sw_query.py \
		--input_file data/seqs/01_list_alpha_subunits.csv

#### Get RPS-Blast (Cdd, Ko, Pfam) hits (2)
Get annotation for initial curated sequences (NBCI) `analyses/blast_outputs/separated_fasta/` and SwissProt query sequences `analyses/sp_query/`

		python 05_do_rps_blast.py \
		-in analyses/blast_outputs/separated_fasta/ \
		-cdd data/seqs/db_cdd/db_cdd \
		-e 0.01 \
		-pfam data/seqs/db_pfam/ \
		-ko data/seqs/db_ko/profiles/ \
		-kl data/seqs/db_ko/ko_list \
		-kc data/seqs/db_ko/config.yml

		python 06_create_annotated_table.py /
		--fasta_input analyses/blast_outputs/separated_fasta/

		---------------------------------------------

		python 05_do_rps_blast.py \
		-in analyses/sp_query/ \
		-cdd data/seqs/db_cdd/db_cdd \
		-e 0.01 \
		-pfam data/seqs/db_pfam/ \
		-ko data/seqs/db_ko/profiles/ \
		-kl data/seqs/db_ko/ko_list \
		-kc data/seqs/db_ko/config.yml

		python 06_create_annotated_table.py /
		-in analyses/sp_query/


#### Format all the analyzed fasta files for further processing

		python 08_do_format_seqs.py \
		--input_fa analyses/sp_query/

		python 08_do_format_seqs.py \
		--input_fa analyses/blast_outputs/

> I would recommend to put `.fa` of `analyses/sp_query/` in a separate directory

#### Generate a table with all the collected fasta files and number of sequences within
Keys: accession (NCBI accession #); sp_query (List of fasta files gathered after direct UniProt keyword query); filtered# (Number of SwissProt query
sequences after filtration of fragments, etc.); unfiltered# (Number of SwissProt raw keyword query sequences); ncbi_ref (Manually curated NCBI files);
sp_ref (Files resulted from RPS Blast); sp_ref# (Number of sequences in each RPS Blast file).

		python 09_prep_phylo.py \
		--initial_table data/seqs/01_list_alpha_subunits.csv \
		--ncbi_fa analyses/blast_outputs/separated_fasta/ \
		--sp_blast analyses/blast_outputs/separated_hits/merged/ \
		--sp_query analyses/sp_query/

> Output: `analyses/annotated_table.tsv`

#### Conclusion
* Curated NCBI sequences can be found at `analyses/blast_outputs/separated_fasta/`
* SwissProt review sequences can be found at `analyses/blast_outputs/separated_hits/relevant/`
* SwissProt outgroup sequences can be found at `analyses/blast_outputs/separated_hits/outgroup/`
* SwissProt outgroup + reviewed sequences can be found at `analyses/blast_outputs/separated_hits/merged/`
* SwissProt keyword query sequences (filtered) can be found at `analyses/sp_query/*.filtered`
* Table with fasta files to be used further in the phylogenetic analyses `analyses/annotated_table.tsv`

-------------------------------------------------------------------------------------------------------------

### Phylogenetic analysis

#### Prepare working environment
Create `analyses/phylo` directory and put there merged `.fasta` files from `analyses/blast_outputs/separated_fasta/`,
`analyses/blast_outputs/separated_hits` and `analyses/sp_query/`.

		python 10_phylo_env.py \
		--input_table analyses/annotated_table/annotated_table.tsv

#### Perform alignment, trimming and generate phylogenetic tree
Perform multiple alignment with MUSCLE. Use Trimal for alignment trimming (gap threshold of 0.05). Use IQ-TREE for phylogenetic tree rendering.
Use UFboot option for bootstrapping.

> In order to get accurate phylogenetic tree. Instead of GAMMA method use LG+R7 (default method in `11_generate_best_tree`)

		python 11_generate_best_tree.py \
		--input analyses/phylo

#### (Optional) Get references for each curate alpha subunits
Left join references titles for each alpha subunit in `analyses/annotated_table.tsv`

		python 12_get_references.py \
		--input_table analyses/annotated_table.tsv \
		--email <insert your email>

-------------------------------------------------------------------------------------------------------------

### Generate report

#### Annotate phylogenetic trees
In order to annotate phylogenetic trees in `analyses/phylo` use R & ggtree (check the R script in `/phylo/*.ipynb`).
> Use for jupyter notebook for `.ipynb`

#### Create pipeline flow chart
Check `/flow_charts` for flow chart creation.

#### Final reports
Check `/analyses/phylo/*.Rmd` for rendering final reports.

--------------------------------------------------------------------------------------------------------------

### Final results
Final files can be accessed in `analyses/phylo/*`.
Final `.html` reports can be accessed in `reports/*`.
