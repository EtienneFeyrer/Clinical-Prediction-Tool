# Annotation

Given is the following topic: Variant prioritization for rare disease (annotations used and why)

## Vocabulary

**Allele**  
A variant of a specific gene at a particular position in the genome.

**Annotation**  
Process of describing the structure and function of the components of a genome.

**Splicing**  
Process where mRNA is edited to remove introns (non-coding locations of DNA).  
A single gene can produce multiple transcripts by including or excluding different exons.

**Gene**  
A DNA template locus on the genome coding a functional product or protein.

**Transcript**  
A variant can have different transcripts (mRNAs), due to several options in splicing.

## NGS-Analysis Pipeline

### Preprocessing Steps

- **Adapter trimming:** `SeqPurge` from `ngs-bits`
- **Mapping:** `BWA mem`
- **Duplicate marking:** `samblaster`

### Variant Calling & Normalization

- **Variant calling:** `DeepVariant`
- **Splitting multi-allelic variants:** `VcfBreakMulti` from `ngs-bits`
- **Left-Alignment:** `VcfLeftAlign` from `ngs-bits`

### Variant-to-Gene Mapping

- **Variant annotation:** `VPE`

### Filtering Process

- **Allele Frequency Filtering:** `GnomAD`
- **Conservation Filtering:** `phyloP`
- **In-House Database Filter:** Phenotype from patient data
- **Quality Filtering:** Data provided by sequencing
- **Pathogenicity & Impact Filtering:** `SIFT` and `PolyPhen-2`
- **Inheritance & Phenotype Filtering:** `Human Phenotype Ontology`, `OMIM`, `TumorPortal`

## 1. Where is Annotation in Our Project

- **Raw Data Input**  
  Start with sequencing reads obtained from a sequencer.

- **Read Processing**  
  The raw reads are processed and aligned to a reference genome.  
  **Annotation Element:**  
  This step assigns each read a genomic coordinate, effectively annotating the raw data with positional information.

- **Variant Calling**  
  The aligned reads are transformed into a Variant Call Format (VCF) file.  
  **Annotation Element:**  
  The VCF file includes detailed information for each genomic position where a variant is detected (if present), along with quality scores and likelihood estimates.

- **Technical Filtering**  
  Technical filters (e.g., allele frequency, quality filtering) are applied to reduce the number of variants.  
  **Annotation Element:**  
  This filtering annotates the variants with technical criteria, refining the dataset into a manageable subset for further analysis.

- **Variant Effect Prediction**  
  Tools like the Variant Effect Predictor (VEP) are used to classify variants, specifically identifying those with no impact on protein function.  
  **Annotation Element:**  
  VEP annotates each variant with its potential biological impact, guiding the prioritization of variants that may have functional consequences.

- **Gene-Level Consolidation**  
  After variant effect prediction, the analysis yields both the variants and the genes they are thought to affect.  
  **Annotation Element:**  
  Aggregating variant-level annotations to the gene level provides a focused list of candidate genes for further clinical assessment.

- **Disease and Inheritance Filtering**  
  A specific filter is then applied using disease and inheritance information (e.g., OMIM for disease data and HPO for phenotype data).  
  **Annotation Element:**  
  This final step further annotates the candidate genes by integrating clinical databases and phenotype associations, prioritizing those that match the expected disease and inheritance patterns.

## 2. Filter and Database for Disease and Inheritance Information

- **Database:**  
  [Online Mendelian Inheritance in Man (OMIM)](https://omim.org/) provides curated disease-related gene information.

- **Filter:**  
  [Human Phenotype Ontology (HPO)](https://hpo.jax.org/) as provided by the `phenotype.hpoa` file, supplies detailed phenotype–disease associations.

### How It Works

1. **Phenotype Matching**
   - **Extraction & Comparison**  
     For each gene affected by the variants (annotated via VEP), compare its associated HPO annotations (from `phenotype.hpoa`) with the patient's clinical features.
   - **Lookup for Detailed Features**  
     If a match is found, perform a lookup in `genes_to_phenotypes.txt` to retrieve additional clinical features linked to that gene.
   - **Purpose**  
     This process shifts the analysis from the variant level to the gene level, aggregating multiple variant effects into a single unit (gene) and aligning this information with known disease phenotypes.

2. **Inheritance Filtering**
   - **Input Specification**  
     Receive the expected inheritance pattern as an input (e.g., autosomal dominant, autosomal recessive, etc.).
   - **Prioritization**  
     From the lookup results in `genes_to_phenotypes.txt`, prioritize those genes whose clinical features (listed in HPO annotations) indicate the specified inheritance pattern.
   - **Benefit**  
     This step leverages the inheritance information available from OMIM and HPO to further narrow down candidate genes based on whether they match the patient’s family history or clinically observed inheritance pattern.

### Summary

By integrating variant-level annotations (from tools like VEP) into a gene-level list and then applying phenotype matching using the `phenotype.hpoa` file along with supplemental gene-to-phenotype lookups, you create a refined candidate list.  
The additional filtering based on inheritance patterns ensures that only genes with clinical features and documented inheritance modes corresponding to the patient’s profile are prioritized.  
This multi-step approach results in a robust and clinically relevant method for variant prioritization in rare, Mendelian diseases. However, instead of using databases of inherited diseases exclusively,  
we can perform the lookup of the VEP's genes in other databases like TumorPortal to check for cancer variants.

### Decision (Choice of Annotation)

Excludes HPO, which is handled separately in the GUI.  
The following papers were helpful in determining the importance of several in silico prediction models:  
- [PMC9774026](https://pmc.ncbi.nlm.nih.gov/articles/PMC9774026/)  
- [Insights on variant analysis in silico tools for pathogenicity prediction](https://www.researchgate.net/publication/365863528_Insights_on_variant_analysis_in_silico_tools_for_pathogenicity_prediction)

#### Effect on cDNA/protein

1. Gene (received from VEP)
2. Transcript (retained from VEP is the transcript of the given gene)
3. Restrict transcripts on MANE
4. Type (Information from technical processing pipeline)
5. Consequence (VEP prediction)
6. LOFtee+

#### Population Genetics

7. Allele frequencies (expected occurrence in population from GnomAD)+

#### Conservedness

8. SIFT()-
9. FATHMM (hidden Markov models)+
10. GERP+

#### Sequence Structure

11. PolyPhen-2 +

#### Supervised ML models

12. REVEL +
13. MutationTaster2021 (SNPRED)
    CADD -

#### Splicing Analysis

14. SpliceAI +

#### Synthome

15. OMIM
16. ClinVar

#### Reasoning

- SIFT is chosen to calculate and predict pathogenicity with the amino acid substitution.
- Mutation Assessor: risk prediction based on 3D structure of homologous sequences.
- FATHMM and KGGSeq are the tools with the highest discriminative power, tested and proven by [Dong et al. (2015)](https://pubmed.ncbi.nlm.nih.gov/25552646/).
- Several studies show VEST3, FATHMM, and REVEL as top performers.
- FATHMM-XF provides highly accurate genome-wide SNV predictions with confidence scores for simplified interpretation and cautious classification ([source](https://www.researchgate.net/publication/365863528_Insights_on_variant_analysis_in_silico_tools_for_pathogenicity_prediction)).
- REVEL uses MutPred, VEST, and PolyPhen.
- SNPRED (2023), as per the article above, is current state of the art, using even more models not only trained on ClinVar but also on BRCA1 (not available and tested, though).
- [MutationTaster2021](https://pubmed.ncbi.nlm.nih.gov/33893808/), published and cited, is one of the newest prediction models tested on real data. It's a random forest model with Mutation Distiller, including data from HPO, OMIM, and Orphanet, and could potentially replace the HPO lookup at the end.