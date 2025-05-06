# Datasets

1. [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98969](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98969)  
   1. Columns:  
      1. Bases  
      2. Bytes  
      3. mouse\_age  
      4. Organ:   
         1. brain   
         2. spinal cord  
      5. selection\_marker  
      6. source\_name  
      7. strain  
      8. treatment  
      9. Version  
   2. DAG columns:   
      1. Experiment Condition: treatment , strain, source\_name  
      2. True Gene Expression  
      3. Sequencing Depth: Bases, Bytes, AvgSpotLen  
      4. Gene Length  
      5. Technical Noise: Platform, Center Name, selection\_marker  
2. [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74995](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74995)  
   1. Columns:  
      1. age  
      2. AvgSpotLen  
      3. Bases  
      4. brain\_region  
      5. Bytes  
      6. ReleaseDate  
      7. sex  
      8. Treatment  
   2. DAG columns:  
      1. Experiment Condition: treatment , age, brain\_region, strain, source\_name  
      2. True Gene Expression  
      3. Sequencing Depth: Bases, Bytes, AvgSpotLen  
      4. Gene Length  
      5. Technical Noise: Platform, Center Name, selection\_marker  
3. [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121654](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121654)  
   1. Columns:  
      1. AGE  
      2. AvgSpotLen  
      3. Bases  
      4. brain\_region  
      5. Bytes  
      6. ReleaseDate  
      7. sex  
      8. treatment  
   2. DAG columns:  
      1. Experiment Condition: treatment , age, brain\_region, strain, source\_name  
      2. True Gene Expression  
      3. Sequencing Depth: Bases, Bytes, AvgSpotLen  
      4. Gene Length  
      5. Technical Noise: Platform, Instrument, selection\_marker  
4. [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123025](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123025)  
   1. Columns:  
      1. AGE  
      2. Bases  
      3. Bytes  
      4. detected\_genes  
      5. gate  
      6. pooled\_library\_names  
      7. qc\_all\_3\_criteria  
      8. qc\_detected\_genes  
      9. qc\_ercc\_correlation  
      10. qc\_total\_counts  
      11. tissue  
      12. Total\_counts  
   2. DAG columns:  
      1. Experiment Condition: age, tissue, strain, gate  
      2. True Gene Expression : detected\_genes, total\_counts, qc\_detected genes  
      3. Sequencing Depth: Bases, Bytes  
      4. Gene Length  
      5. Technical Noise: AvgSpotLen, LibrarySelection, Assay type, qc\_total\_counts, qc\_ercc\_correlation  
      6. Observed Read Counts: total\_counts, qc\_total\_counts  
5. [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123030](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123030)  
   1. Columns:  
      1. AGE  
      2. Bases  
      3. BioProject  
      4. Bytes  
      5. detected\_genes  
      6. gate  
      7. genotype  
      8. qc\_all\_3\_criteria  
      9. qc\_detected\_genes  
      10. qc\_ercc\_correlation  
      11. qc\_total\_counts  
      12. source\_name  
      13. SRA Study  
      14. strain  
      15. tissue  
      16. Total\_counts  
   2. DAG Columns:  
      1. Experiment Condition: age, tissue, strain, gate, genotype, source\_name  
      2. True Gene Expression : qc\_detected\_genes  
      3. Sequencing Depth: Bases, Bytes  
      4. Gene Length  
      5. Technical Noise: qc\_total\_counts, qc\_ercc\_correlation, qc\_detected\_genes, plate\_id  
      6. Observed Read Counts: total\_counts, qc\_total\_counts  
6. [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130627](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130627)  
   1. Columns:  
      1. AGE  
      2. AvgSpotLen  
      3. Bases  
      4. BioProject  
      5. Bytes  
      6. cell\_type  
      7. create\_date  
      8. Diet  
      9. genotype  
      10. Instrument  
      11. LibraryLayout  
      12. sex  
      13. source\_name  
      14. SRA Study  
      15. Timepoint  
   2. DAG Columns:  
      1. Experiment Condition: age, tissue, sex, timepoint, diet, source\_name, cell\_type  
      2. True Gene Expression   
      3. Sequencing Depth: Bases, Bytes  
      4. Gene Length  
      5. Technical Noise: AvgSpotLen, Instrument, LibrayLayout  
         

---

| Variable Type | Variables | How to Obtain / Calculate |
| ----- | ----- | ----- |
| **From Dataset** | `age`, `brain_region`, `sex`, `treatment` | Directly from metadata |
| **Sequencing Depth** | `Bases`, `Bytes`, `AvgSpotLen` | Metadata columns |
| **Observed Read Counts** | `RNA-Seq read counts` | Provided count matrix |
| **Calculated** | `Gene Length` | Extracted from **GTF file** |
| **Calculated** | `TPM` (Transcripts Per Million) | Normalize counts |
| **Estimated** | `Technical Noise` | Use **Negative Binomial / ComBat** |
| **Outcome Variable** | `Observed Read Counts` | Direct measurement |
| **Latent Variable** | `True Gene Expression` | Bayesian inference |

---

### **Most Important Columns**

| DAG Variable | Most Commonly Observed Columns |
| ----- | ----- |
| **Experimental Condition** | `AGE`, `tissue`, `strain`, `genotype`, `sex`, `Diet`, `treatment`, `gate`, `selection_marker` |
| **True Gene Expression** | `detected_genes`, `total_counts` (from count matrices) |
| **Sequencing Depth** | `Bases`, `Bytes`, `AvgSpotLen`, `Total_reads`, `pooled_library_names` |
| **Gene Length** | **Not found in metadata** (needs external annotation) |
| **Technical Noise** | `AvgSpotLen`, `Instrument`, `LibrarySelection`, `LibraryLayout`, `qc_ercc_correlation`, `qc_detected_genes`, `qc_total_counts`, `plate_id` |
| **Observed Read Counts** | `total_counts`, `qc_total_counts` (from count matrices) |

