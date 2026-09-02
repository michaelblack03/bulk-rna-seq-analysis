# Experimental design

The dataset consisted of nine RNA-seq samples, comprising three experimental groups (A, B, C) with three samples per group.

Group labels were assigned for the purpose of the analysis and do not represent specific biological conditions.

Samples were distributed across three experimental batches, with one sample from each group represented in each batch.

|Sample|Experimental group|Batch|
|------|---|-----|
|s1|A|1|
|s2|A|2|
|s3|A|3|
|s4|B|4|
|s5|B|5|
|s6|B|6|
|s7|C|7|
|s8|C|8|
|s9|C|9|

The balanced design allowed differences between experimental groups to be assessed while also investigating potential effects associated with batch.

Differential expression was evaluated for two pairwise comparisons:
- B vs A
- C vs A

Both comparisons were assessed using DESeq2 under two log2 fold change (LFC) criteria: the standard null hypothesis of LFC = 0 and a biological-effect threshold of |LFC| ≥ 1.
