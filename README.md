##############################
# 📓 Julia Bioinformatics Portfolio
# Author: Dr Mostafa Samy
# Purpose: Showcase Julia for Bioinformatics/Omics Analysis
# Environment: Google Colab or local Julia (v1.9+)
##############################

##############
# Intro
##############
#=
# Julia Bioinformatics Portfolio
This notebook demonstrates applications of **Julia** in:
- Transcriptomics
- Sequence analysis
- Structural bioinformatics
- Systems biology & networks
- Genomics
- Clinical linkage

Each section includes:
1. Markdown explanation
2. Code snippet
3. Visualization

=#

##############
# Table of Contents
##############
# 1. Expression Analysis
# 2. Sequence Analysis
# 3. Systems Biology & Networks
# 4. Structural Bioinformatics
# 5. Genomics
# 6. Clinical Linkage
# 7. Reporting & Export
##############

using Random, Statistics, LinearAlgebra, DataFrames, Plots
Random.seed!(123)

########################################################
# 1. Expression Analysis
########################################################

# Simulate expression matrix
n_genes, n_samples = 30, 8
expr_matrix = rand(100:1000, n_genes, n_samples)
expr_log = log.(expr_matrix .+ 1)

# Heatmap
heatmap(expr_log,
        xlabel="Samples", ylabel="Genes",
        yflip=true, c=:viridis,
        title="Gene Expression Heatmap (log scale)")

# PCA (SVD method)
X = expr_log .- mean(expr_log, dims=1)
U,S,Vt = svd(X)
scores = X * Vt'
scatter(scores[:,1], scores[:,2],
        xlabel="PC1", ylabel="PC2", title="PCA of Expression Data")

########################################################
# 2. Sequence Analysis
########################################################

# GC content
seqs = ["ATGCATGC", "GGGGCCCC", "ATATATAT"]
gc_content(seq) = count(ch->ch in ['G','C'], seq)/length(seq)*100
gc_vals = [gc_content(s) for s in seqs]

# DNA → Protein (minimal codon table)
codon_table = Dict("ATG"=>"M","TAA"=>"*","TAG"=>"*","TGA"=>"*")
function translate(seq)
    protein=""
    for i in 1:3:length(seq)-2
        protein *= get(codon_table,seq[i:i+2],"X")
    end
    return protein
end
println(translate("ATGTAA"))

########################################################
# 3. Systems Biology & Networks
########################################################

# Mock GRN
genes = ["G$i" for i in 1:6]
edges = [(rand(genes), rand(genes)) for _ in 1
