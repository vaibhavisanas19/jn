import streamlit as st
import matplotlib.pyplot as plt
import networkx as nx
from Bio import AlignIO, Phylo, pairwise2
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import pandas as pd
import os

# Apply custom CSS for better styling
st.markdown("""
    <style>
        .main {
            background-color: #f0f8ff;
        }
        h1 {
            color: #004d99;
            font-size: 36px;
            font-weight: bold;
            text-align: center;
            text-shadow: 2px 2px 4px #99ccff;
        }
        h2, h3 {
            color: #003366;
            text-align: center;
        }
        .stTextArea, .stButton > button {
            border-radius: 10px;
            padding: 12px;
        }
        .stButton > button {
            background: linear-gradient(to right, #4CAF50, #2E8B57);
            color: white;
            font-weight: bold;
            border: none;
            transition: 0.3s;
            font-size: 18px;
            box-shadow: 3px 3px 6px rgba(0, 0, 0, 0.2);
        }
        .stButton > button:hover {
            background: linear-gradient(to right, #2E8B57, #1E5D3F);
            transform: scale(1.05);
        }
    </style>
""", unsafe_allow_html=True)


# Function to save user sequences to a FASTA file
def save_fasta(user_sequences):
    fasta_file = "sequences.fasta"
    with open(fasta_file, "w") as f:
        f.write(user_sequences)
    return fasta_file


# Function to analyze sequences
def analyze_sequences(file_path, ref_seq):
    alignment_results = []
    alignment = AlignIO.read(file_path, "fasta")

    for record in alignment:
        seq_str = str(record.seq)
        alignments = pairwise2.align.globalxx(ref_seq, seq_str)
        best_alignment = alignments[0]

        similarity = sum(a == b for a, b in zip(best_alignment.seqA, best_alignment.seqB))
        percent_identity = (similarity / len(best_alignment.seqA)) * 100

        snps = [(i + 1, ref, seq_str[i]) for i, ref in enumerate(ref_seq) if seq_str[i] != ref]

        alignment_results.append({
            "Sequence": record.id,
            "Percent Identity": percent_identity,
            "SNPs": str(snps),  # Convert SNPs to string to avoid PyArrow error
            "Best Alignment": (best_alignment.seqA, best_alignment.seqB)
        })

    return pd.DataFrame(alignment_results)


# Function to build a phylogenetic tree
def build_phylogenetic_tree(file_path):
    alignment = AlignIO.read(file_path, "fasta")
    calculator = DistanceCalculator("identity")
    constructor = DistanceTreeConstructor(calculator, "upgma")
    phylo_tree = constructor.build_tree(alignment)

    tree_image_file = "tree.png"
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    Phylo.draw(phylo_tree, axes=ax)
    plt.savefig(tree_image_file)
    plt.close()

    return phylo_tree, tree_image_file


# Streamlit UI
st.title("🧬 Phylogenetic Tree & Sequence Analysis")

# User input for FASTA sequences
st.subheader("📜 Enter your sequences in FASTA format:")
user_input = st.text_area("Example format:",
                          ">Seq1\nATGCGTACGTTAG\n>Seq2\nATGCGTACGTGAG\n>Seq3\nATGCGTTCGTTAG\n>Seq4\nATGCGGACGTTAG",
                          height=200)

if st.button("🚀 Generate Phylogenetic Tree & Analyze Sequences"):
    if user_input.strip() == "":
        st.error("⚠️ Please enter valid sequences in FASTA format.")
    else:
        fasta_file = save_fasta(user_input)
        st.success("✅ FASTA File Created Successfully!")

        phylogenetic_tree, tree_image = build_phylogenetic_tree(fasta_file)
        if os.path.exists(tree_image):
            st.image(tree_image, caption="🌳 Phylogenetic Tree", use_column_width=True)
        else:
            st.error("❌ Error generating tree")

        st.subheader("📊 Sequence Analysis Results")
        sequence_analysis_results = analyze_sequences(fasta_file, user_input.split("\n")[1])
        st.dataframe(sequence_analysis_results)
