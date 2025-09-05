# This app is created by Chanin Nantasenamat (Data Professor) https://youtube.com/dataprofessor
# Credit: Inspired by https://huggingface.co/spaces/osanseviero/esmfold

import streamlit as st
from stmol import showmol
import py3Dmol
import requests
import biotite.structure.io as bsio
import re

# Sidebar info
st.sidebar.title("🎈 ESMFold")
st.sidebar.write(
    "[*ESMFold*](https://esmatlas.com/about) is an end-to-end single sequence protein structure predictor "
    "based on the ESM-2 language model. For more information, read the "
    "[research article](https://www.biorxiv.org/content/10.1101/2022.07.20.500902v2) and the "
    "[news article](https://www.nature.com/articles/d41586-022-03539-1) published in *Nature*."
)

# -------------------
# Sequence cleaner
# -------------------
def clean_sequence(seq: str) -> str:
    seq = seq.upper()
    seq = re.sub(r"^>.*", "", seq, flags=re.MULTILINE)  # remove FASTA headers
    seq = re.sub(r"[^A-Z]", "", seq)  # remove spaces, numbers, newlines
    valid = set("ACDEFGHIKLMNPQRSTVWYBXZJ")  # standard AA + ambiguity codes
    seq = "".join([aa for aa in seq if aa in valid])
    return seq

# stmol renderer
def render_mol(pdb):
    pdbview = py3Dmol.view()
    pdbview.addModel(pdb, "pdb")
    pdbview.setStyle({"cartoon": {"color": "spectrum"}})
    pdbview.setBackgroundColor("white")
    pdbview.zoomTo()
    pdbview.zoom(2, 800)
    pdbview.spin(True)
    showmol(pdbview, height=500, width=800)

# Default sequence
DEFAULT_SEQ = (
    "MGSSHHHHHHSSGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCSLEDPAANKARKEAELAAATAEQ"
)

# Sidebar input
txt = st.sidebar.text_area("Paste sequence", DEFAULT_SEQ, height=275)

uploaded_file = st.sidebar.file_uploader("Or upload FASTA file", type=["fasta", "fa", "txt"])
if uploaded_file is not None:
    content = uploaded_file.read().decode("utf-8")
    txt = content

# Prediction function
def update(sequence=txt):
    sequence = clean_sequence(sequence)

    if not sequence:
        st.error("❌ No valid sequence after cleaning. Check your input.")
        return

    headers = {"Content-Type": "application/x-www-form-urlencoded"}
    response = requests.post(
        "https://api.esmatlas.com/foldSequence/v1/pdb/", headers=headers, data=sequence
    )

    pdb_string = response.content.decode("utf-8")

    # Check validity
    if not pdb_string.strip().startswith("HEADER"):
        st.error("❌ Failed to fetch valid PDB from ESMFold API. Response was:")
        st.code(pdb_string[:500])  # Show first 500 chars
        return

    with open("predicted.pdb", "w") as f:
        f.write(pdb_string)

    # Try loading with Biotite
    try:
        struct = bsio.load_structure("predicted.pdb", extra_fields=["b_factor"])
        b_value = round(struct.b_factor.mean(), 4)
    except Exception as e:
        st.warning(f"⚠️ Could not parse PDB with Biotite: {e}")
        b_value = "N/A"

    # Visualization
    st.subheader("Visualization of predicted protein structure")
    render_mol(pdb_string)

    # plDDT
    st.subheader("plDDT")
    st.write("plDDT is a per-residue estimate of the confidence in prediction on a scale from 0–100.")
    st.info(f"plDDT: {b_value}")

    # Download button
    st.download_button(
        label="Download PDB",
        data=pdb_string,
        file_name="predicted.pdb",
        mime="text/plain",
    )

# Sidebar button
predict = st.sidebar.button("Predict", on_click=lambda: update(txt))

if not predict:
    st.warning("👈 Paste a protein sequence or upload a file, then click **Predict**.")
