import matplotlib.pyplot as plt
import streamlit as st
from analyzer import *

st.title("DNA Sequence Analyzer")
user_sequence = st.text_input("Enter a DNA sequence:")

if user_sequence:
    if validate_sequence(user_sequence):
        st.write("✅ Valid DNA sequence")

        st.write("Base counts:")
        base_counts = count_bases(user_sequence)
        for base, count in base_counts.items():
            st.write(f"{base}: {count}")

        gc_percent = gc_content(user_sequence)
        st.write(f"GC content: {gc_percent:.2f}%")

        reverse_comp = reverse_complement(user_sequence)
        st.write(f"Reverse complement: {reverse_comp}")

        rna_seq = transcribe(user_sequence)
        st.write(f"Transcribed RNA sequence: {rna_seq}")

        protein_seq = translate(user_sequence)
        st.write(f"Translated protein sequence: {protein_seq}")

        labels = list(base_counts.keys())
        values = list(base_counts.values())
        fig, ax = plt.subplots()
        ax.bar(labels, values)
        ax.set_xlabel("Base")
        ax.set_ylabel("Count")
        ax.set_title("Nucleotide Composition")
        st.pyplot(fig)

        st.write("GC vs AT Content:")
        gc = base_counts["G"] + base_counts["C"]
        at = base_counts["A"] + base_counts["T"]
        fig2, ax2 = plt.subplots()
        ax2.pie([gc, at], labels=["GC", "AT"], autopct='%1.1f%%')
        ax2.set_title("GC vs AT Content")
        st.pyplot(fig2)

        start_positions = find_start_codons(user_sequence)
        fig3, ax3 = plt.subplots()
        for index, start in enumerate(start_positions):
            orf = find_orf_from_start(user_sequence, start)
            ax3.barh(index, len(orf), left=start)
        ax3.set_xlabel("Position in sequence")
        ax3.set_ylabel("ORF")
        ax3.set_title("ORF Locations")
        st.pyplot(fig3)

        st.write("Six-Frame ORF Search:")
        six_frame_results = find_all_orfs_six_frames(user_sequence)
        for frame_label, orfs in six_frame_results.items():
            st.write(f"{frame_label}: {orfs}")

    else:
        st.write("❌ Invalid DNA sequence")





