Goal: Given a sample of genes in a cell, I want to identify the top 20 most important representative non-canonical marker genes by building a model that can identify whether this cell is more OL-leaning or NB-leaning. To achieve that, this model uses all Highly Variable non-canonical, background markers (1,971 genes) as features and predicts the probability of the cell being an OL or an NB. In other words, it distinguishes a cell’s identity as OL leaning or NB leaning. As an application, I want to apply this model to an out-of-sample dataset and observe the accuracy of the model in the identification of OL or NB. Finally, I want to apply this model to a TAP dataset, which is completely excluded from the training. Here, I can compare the Marker Panel Scoring (using relative expression of OPC/COP/OL compared to the expression of NB) to the Machine Learning Classifier score.

Part 1: Conceptual Basis (Reference 1)

Whenever the brain undergoes white matter injury, or demyelination, the brain needs to produce OLs, or Oligodendrocytes, which go to the injury site and repair the damaged myelin sheath using a process called remyelination. In the brain, there are NSCs, which are neural stem cells, and they can be divided into either active NSCs or quiescent NSCs. These NSCs usually form TAPs (Transit Amplifying Progenitors), which are the ones that divide into either Neuronal Cells (which become neurons or other stuff in the brain) or OLs. Usually, TAPs become Neuronal, but the scientists in the paper that I’ve based my research on identify two growth factors/chemicals/ligands that are secreted by the microglia, which are a part of the brain’s immune cells in the V-SVZ niche, which serve as signals to the TAPS and influence them to become OLs. These two factors are IGF1 (insulin growth factor 1) and OSM (oncostatin-M). 

Part 2: Preprocessing + Understanding the Primary Dataset

Step 1  
INPUT · 10x scRNA-seq (GSE266687)  
10 samples × {barcodes.tsv, genes.tsv, matrix.mtx}  
10 GSMs, 4 conditions (Cntl/CupRap × 0wks/3wks)

The primary dataset used in this project is called GSE266687. Here, there are 10 samples of mouse V-SVZ scRNA-seq data. There are three different variables that are varied across the 10 samples: a strain, a treatment, and a time point. There are two types of Strain: CD1 and NesCre, where CD1 represents the normal, wild type of mice, whereas NesCre represents the Neuronal Stem cell targeting line. There are two types of treatment: Cntl and CupRap, where Cntl means a control treatment where everything is normal, whereas CupRap stands for Cuprizone Rapamycin, which injures the mouse and doesn’t let it naturally remyelinate (Cuprizone demyelinates the oligodendrocytes, while Rapamycin inhibits the mTOR pathway, so it inhibits natural remyelination response - Reference 1). There are two types of timepoints, 0 weeks or 3 weeks, where 0 weeks means that the mouse is given 0 weeks to proliferate/recover, whereas 3 weeks means that the mouse is given 3 weeks to proliferate/recover. In my analysis and for my purposes, 3 weeks vs 0 weeks will not be modeled explicitly because the model will not use timepoint as a feature, is trained to distinguish identity, not state. The fully formatted table with the 10 samples can be found here: [mouse V-SVZ scRNA-seq](https://docs.google.com/spreadsheets/d/1i39wCQLV8vi30uBefOcpt6HsXY3ePkHoAEIkNeEBZZs/edit?gid=0#gid=0) (link: [mouse V-SVZ scRNA-seq](https://docs.google.com/spreadsheets/d/1i39wCQLV8vi30uBefOcpt6HsXY3ePkHoAEIkNeEBZZs/edit?usp=sharing) )  
Step 2  
SCANPY · ingest (data_io.py)

* read each sample → AnnData (cells × genes, sparse counts)  
* code: adata = data_io.load_all(RAW)  
* obs metadata: strain, treatment, timepoint, replicate, sample_id  
* concatenate 10 samples on inner gene intersection → one AnnData

- There are two important files: barcodes and genes. Barcodes serve as an identity to each cell, and the genes are the gene expression values gleaned from each cell. The barcode is placed on the row labels, the genes are placed on the column labels, and a table is created with the gene expression values.   
- Everything was saved to the CSV raw_barcodes_genes_top100.csv. This can be found here: [raw_barcodes_genes_top100.csv](https://drive.google.com/file/d/1FgqykW-bC-EvUJ0CogPWKodv_uCvF8C0/view?usp=sharing) 

Step 3  
SCANPY · QC (preprocess.py)

* per-cell: n_genes_by_counts ≥ 500, pct_counts_mt ≤ 10  
* per-gene: expressed in ≥ 10 cells  
* 52,469 cells × 26,599 genes → 47,965 cells × 19,139 genes (but scanpy only keeps 2,000 HVG)

- Quality control is then performed on this data (which is all done by Scanpy - Reference 4), removing noise and ensuring that the resulting data are usable and can be applied.   
- This is saved to normalized_barcodes_genes_top100.csv. This can be found here: [normalized_barcodes_genes_top100.csv](https://drive.google.com/file/d/1ixTDt2PSAtAyBZBo_f_x72SF4ZLoVJZy/view?usp=sharing) 

Step 4 (Reference 4 - All Scanpy methods come from this reference: https://scanpy.readthedocs.io/en/stable/how-to/cell-cycle.html)
SCANPY · normalize → embed → cluster → annotate  
(preprocess.py + markers.py, orchestrated by run_pipeline.py)

* normalize_total(target_sum=1e4) + log1p   ``sc.pp.normalize_total(adata, target_sum=1e4)``, ``sc.pp.log1p(adata)``
* highly_variable_genes (2,000, Seurat flavor, batch_key=sample_id)  ``sc.pp.highly_variable_genes(
        adata, n_top_genes=n_hvg, flavor="seurat", batch_key="sample_id"
    )``
* scale → PCA(30) → neighbors(15) → UMAP ``sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=n_pcs)
    sc.pp.neighbors(adata_hvg, n_pcs=n_pcs, n_neighbors=15)
    sc.tl.umap(adata_hvg)``
* leiden(resolution=0.8)  ``sc.tl.leiden(
        adata_hvg,
        resolution=resolution,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
``
* score_genes (sc.tl.score_genes) (Reference #4) per marker panel (12 cell types) →  
  {qNSC, aNSC, TAP, Neuroblast, OPC, COP, OL,  
  Astrocyte, Microglia, Ependymal, Endothelial, Mural}  ``sc.tl.score_genes(
            adata,
            gene_list=present,
            score_name=f"score_{name}",
            ctrl_size=min(50, max(10, len(present) * 5)),
            random_state=0,
        )
``
* convert score for each panel to cell_type, then classify fate_OL_lineage

| score_qNSC | score_aNSC | score_TAP | score_Neuro | score_OPC | score_COP | score_OL | score_Astroc | score_Microg | score_Epend | score_Endotl | score_Mural | cell_type |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| -0.0245925 | 0.15729764 | 0.49285538 | 1.10410176 | -0.238148 | -0.1218039 | -0.6998098 | -0.3862443 | -0.8701482 | -0.3226306 | -0.1220567 | -0.2404025 | Neuroblast |
| -0.012597 | -0.2261938 | -0.4950369 | -0.6583718 | 0.59983178 | 0.28798455 | 3.10097177 | -0.3243312 | -0.6981846 | -0.0509202 | -0.162698 | -0.4242628 | OL |


- Since any combo of genes that are part of any of the twelve panels can be found in these cells, a concept called Leiden clustering (Reference 4) needs to be performed in order to find the relative expression of every gene to the other.  
- After this was completed, I had to look at the literature to find the canonical markers. There are 12 (technically 11, because aNSC - active NSC - folds into TAP/qNSC) types of V-SVZ cells in this niche (a few of which we don’t care about - specifically, we only care about 5). I first needed to figure out what the canonical, well-known, representative genes of each of these cells are (organized in Table 1), which was found in the actual literature. After this is determined, each cell type is classified by the canonical representative genes, and finally, a data table results (organized in Table 2):  
- This is all saved to barcode_to_cell_type_mapping.csv. This can be found here: [barcode_to_cell_type_mapping_table.xlsx](https://docs.google.com/spreadsheets/d/1t0ZHQz5KMhODr9L-O4rOQ82MiW5zbme2/edit?usp=sharing&ouid=111504262478734653360&rtpof=true&sd=true). This is a Google Drive sheet that has everything worth mentioning.  

Table 1: Canonical Marker Genes

| Cell Type | Marker Genes   |
| :---- | :---- |
| **qNSC** | Gfap, Slc1a3, Aldh1l1, Vim, Sox2, Sox9, Hes5, Aldoc, Apoe, Id3 |
| **aNSC** | Egfr, Ascl1, Mki67, Top2a, Sox2, Hes5, Nes |
| **TAP** | Egfr, Mki67, Ccnd2, Dlx1, Dlx2, Ascl1, Top2a |
| **Neuroblast** | Dcx, Sox4, Sox11, Stmn2, Stmn3, Nrxn3, Tubb3, Dlx1, Dlx2 |
| **OPC** | Pdgfra, Cspg4, Olig1, Olig2, Sox10, Cnp, Lhfpl3 |
| **COP** | Gpr17, Sox10, Olig2, Itpr2, Bmp4, Bcas1, Tcf7l2 |
| **OL** | Mbp, Mog, Mag, Plp1, Mobp, Cldn11, Mal, Trf, Cnp |
| **Astrocyte** | Aqp4, Slc1a2, Aldh1l1, S100b, Mfge8, Ntsr2, Slc1a3 |
| **Microglia** | Cx3cr1, Tmem119, P2ry12, C1qa, Csf1r, Trem2, Hexb, Aif1 |
| **Ependymal** | Foxj1, Ccdc153, Tmem212, Rsph1, Sox2 |
| **Endothelial** | Cldn5, Pecam1, Flt1, Cdh5 |
| **Mural** | Pdgfrb, Acta2, Rgs5 |

References 6-17 (This is from [markers.py](http://markers.py) - what the canonical markers are)

Table 2: Cell Type Classification by Marker Genes

| Cell type | N amount of cells |
| :---- | :---- |
| Microglia | 16270 |
| Neuroblast | 10205 |
| OL | 6400 |
| TAP | 5202 |
| Mural | 2713 |
| Astrocyte | 1784 |
| Ependymal | 1586 |
| qNSC | 1387 |
| Endothelial | 1304 |
| OPC | 862 |
| COP | 252 |

This table can be found in this Google Drive link: [barcode_to_cell_type_mapping_table.xlsx](https://docs.google.com/spreadsheets/d/1t0ZHQz5KMhODr9L-O4rOQ82MiW5zbme2/edit?usp=sharing&ouid=111504262478734653360&rtpof=true&sd=true). This table can be found in the Google Sheet tab (at the bottom) stating “Cell Type Classification by Marker Genes”. This was found by using the Pivot feature on Excel. 

Part 3: Classifier Design

XGBoost · binary classifier - Reference 2  
(train_ol_classifier.py)

* positive = OPC + COP + OL (n = 7,514), negative = Neuroblast (n = 10,205)  
* held-out sample = GSM8253799  
* features = HVGs - 29 known lineage markers  
* SHAP TreeExplainer → ranked feature importance (Reference 3)

How the Model Works (High Level)

1) First, Scanpy is used to fit the [markers.py](http://markers.py) Python file, which has the 12 labels/cells and their representative genes, into the main CSV file raw_barcodes_genes_top100.csv. The main CSV file has the barcodes/cell identity on the y-axis, and the genes on the x-axis. After the representative genes are fit onto this Python file, the resulting file barocde_to_cell_type_mapping.csv has the complete list of every cell_id matched to the cell_type. It is important to note that this process is done on the full 19k genes, not just the top 2,000 HVGs.   
2) After this is created, the bias score for every cell is created. This process is elaborated on in a little bit. The model’s results are compared to these scores.   
3) Now, the training features are chosen. These include the top 2,000 HVGs minus the 29 canonical markers, so 1,971 genes in total. These 29 canonical markers only belong to the OPC/COP/OL/NB lineages. All of the other canonical markers for the other cell types can be included as part of the features. The main CSV being used here is normalized_barcodes_genes_top100.csv, which has all of the normalized values.   
4) Now, the targets have to be created, which are saved as ol_commitment. Here, for every one of the cells, a number which is indicative of it becoming an OL (1) or NB (0) is created. It’s also marked as being positive for OL leaning and negative for NB leaning. The model now attempts to fit the features into these targets using the binary decision trees with the XGBoost classifier (Reference 2). After the trees have been created, SHAP (Reference 3) tests out the model and determines which are the most important genes that the model depends on. This is explained in more detail in a bit.  
5) The held-out sample is now confirmed using this model. This is also explained in the following section. 

Actual Features of the Model: As described in Part 3, the actual features of this model are the 1,971 genes. These are the non-canonical genes from the 2,000 gene list that surfaced. And these are the main features that are put into the model. 

Targets of the Model: y = 1 (OPC/COP/OL) or y = 0 (NB)  
During this section, I am actually training the classifier to match the 1971 background features to 1 and 0 targets. Again, this is just training the classifier.

Output of the Model: Two outputs (we are actually running the model during parts a and b)

1) A probability score. This gives the probability of a cell becoming this or that (closer to 1, closer to OL, closer to 0, closer to NB). This is run on the TAPs. This is the first deliverable of the project. This “probability score” is the same as saying a classifier score for similarity to the OL lineage class vs the neuroblast class.   
2) SHAP (Ref 3) Top 20 Genes: This is the ranked list of genes SHAP most heavily relied on. Basically, since XGBoost is unreadable to humans, and SHAP is based on game theory, SHAP played a game of “what if” and systematically hid genes from XGBoost to see how much the prediction collapsed, basically seeing how much importance XGBoost was putting onto a specific gene. This top 20 gene list is the second deliverable. 

**Data Table 3: XGBoost (Reference 2) Feature Importance: Ranked Top 20 Gene List**

| rank | gene | SHAP | gain | rough role   |
| :---: | :---- | ----: | ----: | :---- |
| **1** | **Igfbpl1** | 0.507 | 1880 | IGF-binding-protein-like; regulates IGF1 signaling — *biggest discriminator by ~3×* |
| **2** | **Pllp** | 0.378 | 583 | plasmolipin; myelin lipid raft — early OL-lineage marker |
| **3** | **Dlx6os1** | 0.370 | 1064 | lncRNA antisense to Dlx6; inhibitory-neuron / neuroblast marker — *negative predictor* |
| **4** | **Meis2** | 0.352 | 443 | striatal-neuron TF — *negative predictor* |
| 5 | Gjc3 | 0.129 | 82 |  |
| 6 | Meg3 | 0.093 | 90 |  |
| 7 | Cryab | 0.084 | 33 |  |
| 8 | Tubb2b | 0.060 | 53 |  |
| 9 | Celf4 | 0.053 | 48 |  |
| 10 | Ugt8a | 0.052 | 52 |  |
| 11 | Bcl11a | 0.046 | 41 |  |
| 12 | Tspan2 | 0.045 | 21 |  |
| 13 | Arx | 0.042 | 81 |  |
| 14 | Dock10 | 0.041 | 13 |  |
| 15 | Tmsb10 | 0.039 | 9 |  |
| 16 | Cd24a | 0.038 | 45 |  |
| 17 | Sp9 | 0.037 | 24 |  |
| 18 | Gad2 | 0.031 | 30 |  |
| 19 | Tnr | 0.028 | 29 |  |
| 20 | Grin2b | 0.026 | 36 |  |

Data Table 3 is the top 20 genes found by SHAP. (Reference 18) This is found in trigger_genes.csv. 

There are two tests for the model. These results are corroborated by using a hand-rolled bias score calculator, which uses the 29 canonical, well-known genes in the literature, and they agree 82% of the time. Another way of measuring how accurate this model is was by using a held-out sample (GSM8253799), where the cell fate commitment is known, and running the cells through the model. The accuracy is 99.85%. 

1) Model trained on cells whose fate is known.  
   1) Then, using a dataset that is excluded from the actual training, see what percentage overlaps. Basically, you run cells whose fate is already known into the model and see if the model gets it correct. This gives a 99.85% accuracy, so you know that the model is very good at identifying both extremes. The model only missed 3 cells from the OL-lineage class, and only 1 cell from the Neuroblast class, as seen in Data Table 4.  
   2) Data Table 4 - Accuracy of Model in Identification of OL or NB. This was saved to out_of_sample_comparision_test1.csv. Link: [out_of_sample_comparison_test1.xlsx](https://docs.google.com/spreadsheets/d/1wljYiTn6uA1lT9IU0YxxoFPw6pmxn-6C/edit?usp=sharing&ouid=111504262478734653360&rtpof=true&sd=true). In the second column, the correct values refer to the p_test variable, which was calculated from the XGBoost classifier, and the total values come from the barcode_to_cell_type_mapping_table.csv portion, serving as an answer key for this test.

| *fate_OL_lineage* | *ol_nb* | Count of label |
| ----- | ----- | ----- |
| 0 | 0 | 1518 |
|  | 1 | 1 |
| 0 Total |  | 1519 |
| 1 | 0 | 3 |
|  | 1 | 1210 |
| 1 Total |  | 1213 |
| **Grand Total** |  | **2732** |

   3) This dataset is a result of taking the 1213 OL-lineage cells and the 1519 neuroblast cells and running them through the model. Since we already know what they will become, we can see how many of the cells the model classified incorrectly, and for how many of the cells the model classified them correctly. The total amount of “0”s, or OLs, is 1519, and the total amount of “1”s, or NBs, is 1213. Here, the model misidentified only 1 NB (as it thought it was an OL), and it misidentified only 3 OLs (as it thought they were neuroblasts). In total, it was able to get an accuracy of 1518/1519 and 1210/1213, which in total gave a 99.85% accuracy for the model.   
2) The model works on the TAPs.  
   1) After each of the TAPs is classified (on a scale of 0-1), this result is then compared with the result given by a different scoring system, called the hand-rolled bias score. This system uses the 29 representative genes, so it doesn’t use the non-canonical data that surfaced during the XGBoost + SHAP Classification (Reference 2, 3) part. Specifically, it only uses the representative genes for the following 4 cells: OPC, COP, OL, and NB. This system takes the mean score value of a TAP becoming an OL, COP, or OPC minus the mean score value of it becoming an NB, so the value that is given is significant in the sense that if it is positive, it is more OL leaning, but if it is negative, it is more NB leaning. (example: if y-x is positive, then y has more influence than x, but if it is negative, x has more influence). These mean score values of it becoming one cell or the other are calculated using Scanpy. It calculates an expression score for each one, and that was what allowed us to compare the two systems. It can calculate these mean score values using this method:  
      1) It calculates the average expression level of each gene in each panel  
      2) Scanpy looks at the rest of the genome and selects a control group amongst the random background genes, which should have a similar expression level to each representative gene.   
      3) Then, it subtracted the representative score and the background score. If the difference is around 0, there is a high chance that the representative genes are at the same level as the background genes, meaning that the cell is probably not what you think it is. For example, if you calculate the mean expression values for all the representative genes for an OPC in Cell A, and calculate the mean expression values of the background genes in Cell A, and if you subtract these two values and get 0, then the background noise is equal to the OPC noise, meaning that the cell would most likely not be an OPC. Think of it like this: if the processes a certain cell should participate in are at the same level as the rest of the background processes, then that means the cell isn’t doing those specific activities and defining itself by them, meaning it isn’t the cell you thought it was.   
         1) However, if the difference is highly positive, then the representative genes are being expressed at a much higher rate compared to the background genes, so there is a higher chance that this cell is what you think it is. In the same example, if the OPC expression values are much greater than the background gene expression values, you know the cell is most likely an OPC.   
            1) It’s also important to note that these background genes aren’t the other 1,971 genes, it is actually the other canonical genes. It’s like a tug of war between the gene sets, seeing which specific OL lineage (either OL, OPC, or COP) or NB lineage is winning based on the sign and magnitude of the difference.   
      4) This process is then repeated for each of the 4 marker panels. This system allows you to create a bias score for each TAP, indicative of what path it will most likely take.  
   2) These two systems are compared, and they agree roughly 82.1% (410 + 16 / 519) of the time. It is important to note that there is no answer key; we are just confirming two separate ways of doing things with each other. We know that by having a score of 82%, meaning that these two independently-derived systems agree with each other 82% of the time, the derivative work we have done using the model matches the work done directly. Table #5 sums it up. Only 87 + 6 cells were disagreed upon, and the bulk of the disagreements (87) came from when the bias score system thought the TAP was heading towards an NB fate, while the model gave a probability saying the TAP was heading towards an OL fate. 

Table #5 - Agreement and Disagreement between both methods  
This is saved to out_of_sample_tap_comparision_test2. Link: [out_of_sample_tap_comparison_test2.xlsx](https://docs.google.com/spreadsheets/d/1PK6RPx_tWkVYlcZlKHx57dc3PnV4sZxQ/edit?usp=sharing&ouid=111504262478734653360&rtpof=true&sd=true) 

| *bias_is_ol* | *prediction_is_ol* | Count of Barcodes | Agree | 426 |
| ----- | ----- | ----- | :---- | ----- |
| 0 | 0 | 410 | Total | 519 |
|  | 1 | 87 | Percent | 82% |
| 0 Total |  | 497 |  |  |
| 1 | 0 | 6 |  |  |
|  | 1 | 16 |  |  |
| 1 Total |  | 22 |  |  |
| **Grand Total** |  | **519** |  |  |

In this diagram, whenever bias_is_ol and prediction_is_ol match, or give the same number, this means the Bias score and the prediction score match. Meaning, if they both give 0, they both think the TAP will become a Neuroblast, but if they both give 1, they think that the TAP will become an OL. These two systems disagree with one another when they give different numbers. The percentage was calculated by finding the total amount they agreed on (410+16) and dividing this by the total amount (519) which were tested, giving a grand total of 82%.

A summary of Marker Panel Scoring vs Machine Learning Classifier Scoring:

Cell Classification Flowchart (Table #6):

| 19,139 genes |  |
| ----- | :---- |
| Path 1: Marker Panel Scoring | Path 2: Machine Learning Classifier   |
| Step 1. Input Features 29 marker genes | Step 1. Input Features 1,971 non–marker HVGs (markers excluded) |
| Step 2. Method sc.tl.score_genes per panel  | Step 2. Method XGBoost classifier (Reference 2) |
| Step 3. Intermediate Metrics score_OPC, COP, OL, Neuroblast | Step 3. Classification Output P(OL) = 0.804 |
| Step 4. Final Value bias = +0.128 | — |
| Agreement check (bias > 0 vs P(OL) ≥ 0.5) |  |
| "both OL" quadrant |  |

In this diagram, this flowchart is run 1 time, and both of the sides agree that it is in the OL quadrant, and this flowchart has to run for all of the other 518 TAPs (as there were 519 in total). To follow the flowchart, start from the cell’s 19,139 genes and either go left or right, as the left route follows the hand-rolled bias scoring system (Path 1), while the right route follows the model system based on probability and similarity (Path 2).   
One thing that is important to understand is that the model is never trained on the TAPs; it is only trained on the OPCs, COPs, OLs, and the Neuroblasts. It can only predict the TAPs. However, we only use 519 TAPs in the actual testing run in bullet point 2 because these TAPs come from a dataset that is untouched by the model. There are 10 datasets in total, and the model is trained on all the relevant cell populations from 9 of them. The remaining dataset, which is GSM8253799, is untouched by the model during training. Even though the other TAPs are also technically untouched, the model could’ve learned batch effects and mouse-specific things for the other 9 ones, so those TAPs into the mix would be unfair and would inflate the accuracy. 