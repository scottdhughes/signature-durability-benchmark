#!/usr/bin/env python3
"""Generate synthetic cohort expression matrices for the signature-durability-benchmark.

Creates biologically realistic gene-by-sample expression matrices with
case/control phenotype labels. Each cohort models expected biology:
- Durable signatures show real effects in appropriate cohorts
- Brittle signatures show effects in at most 1 cohort
- Confounded signatures are driven by confounder gene overlap
- Stealth confounded signatures have both real + confounder signal

Usage:
    uv run python scripts/generate_cohort_data.py --seed 42
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
CONFIG_DIR = PROJECT_ROOT / "config"
DATA_DIR = PROJECT_ROOT / "data"
CURATION_DIR = DATA_DIR / "curation"
FREEZE_DIR = DATA_DIR / "freeze"
MATRIX_DIR = FREEZE_DIR / "cohort_matrices"
PHENO_DIR = FREEZE_DIR / "cohort_phenotypes"

# ---------------------------------------------------------------------------
# Real human gene symbols for background universe
# These are common, well-annotated HGNC symbols spanning diverse functions.
# We use a curated list of ~8200 real genes to avoid placeholder names.
# ---------------------------------------------------------------------------

# fmt: off
BACKGROUND_GENES = [
    # Transcription factors (200)
    "TP63", "TP73", "MYC", "MYCN", "MYCL", "MAX", "MAD1L1", "MXD1",
    "MNT", "MLX", "MLXIP", "MLXIPL", "HEY1", "HEY2", "HEYL", "HES1",
    "HES5", "HES6", "HES7", "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4",
    "JAG1", "JAG2", "DLL1", "DLL3", "DLL4", "RBPJ", "MAML1", "MAML2",
    "MAML3", "SOX2", "SOX4", "SOX9", "SOX11", "SOX17", "POU5F1", "NANOG",
    "KLF4", "KLF2", "KLF5", "KLF9", "KLF10", "KLF11", "KLF13", "KLF15",
    "SP1", "SP3", "SP4", "EGR2", "EGR3", "EGR4", "ETS1", "ETS2",
    "ERG", "FLI1", "ELF1", "ELF3", "ELF4", "ELF5", "ETV1", "ETV4",
    "ETV5", "ETV6", "GATA1", "GATA2", "GATA3", "GATA4", "GATA5", "GATA6",
    "FOXA1", "FOXA2", "FOXA3", "FOXB1", "FOXC1", "FOXC2", "FOXD1", "FOXD3",
    "FOXE1", "FOXF1", "FOXF2", "FOXG1", "FOXH1", "FOXI1", "FOXJ1", "FOXK1",
    "FOXK2", "FOXL1", "FOXL2", "FOXM1", "FOXN1", "FOXN3", "FOXO1", "FOXO3",
    "FOXO4", "FOXQ1", "FOXR1", "FOXS1", "PAX2", "PAX3", "PAX5", "PAX6",
    "PAX7", "PAX8", "PAX9", "TBX1", "TBX2", "TBX3", "TBX4", "TBX5",
    "TBX6", "TBX15", "TBX18", "TBX19", "TBX20", "TBX21", "EOMES",
    "BRACHYURY", "HAND1", "HAND2", "TWIST2", "TCF3", "TCF4", "TCF7",
    "TCF7L1", "TCF7L2", "LEF1", "ID1", "ID2", "ID3", "ID4", "SMAD1",
    "SMAD2", "SMAD3", "SMAD4", "SMAD5", "SMAD6", "SMAD7", "SMAD9",
    "RUNX1", "RUNX2", "RUNX3", "CBFB", "CEBPA", "CEBPB", "CEBPD",
    "CEBPE", "CEBPG", "NFE2", "NFE2L1", "NFE2L2", "NFE2L3", "BACH1",
    "BACH2", "MAFB", "MAFG", "MAFK", "MAFF", "NR1H2", "NR1H3", "NR1H4",
    "NR2F1", "NR2F2", "NR2F6", "NR3C1", "NR3C2", "NR4A1", "NR4A2",
    "NR4A3", "NR5A1", "NR5A2", "PPARA", "PPARD", "PPARG", "RXRA",
    "RXRB", "RXRG", "RARA", "RARB", "RARG", "VDR", "ESR1", "ESR2",
    "AR", "PGR", "THRA", "THRB", "NR1I2", "NR1I3", "HNF1A", "HNF1B",
    "HNF4A", "HNF4G", "TFAP2A", "TFAP2B", "TFAP2C", "TFAP2D", "TFAP2E",
    "MITF", "TFE3", "TFEB", "TFEC", "USF1", "USF2", "CLOCK", "BMAL1",
    "PER1", "PER2", "PER3", "CRY1", "CRY2",
    # Kinases (200)
    "ABL1", "ABL2", "ACK1", "AKT1", "AKT2", "AKT3", "ALK", "AXL",
    "BRAF", "BTK", "CAMK1", "CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G",
    "CAMK4", "CASK", "CDC42BPA", "CDC42BPB", "CDK3", "CDK4", "CDK5",
    "CDK6", "CDK7", "CDK8", "CDK9", "CDK10", "CDK11A", "CDK11B", "CDK12",
    "CDK13", "CDK14", "CDK16", "CDK17", "CDK18", "CDK19", "CDK20", "CHEK1",
    "CHEK2", "CLK1", "CLK2", "CLK3", "CLK4", "CSK", "CSF1R", "DAPK1",
    "DAPK2", "DAPK3", "DDR1", "DDR2", "DYRK1A", "DYRK1B", "DYRK2",
    "DYRK3", "DYRK4", "EGFR", "EPHA1", "EPHA2", "EPHA3", "EPHA4",
    "EPHA5", "EPHA6", "EPHA7", "EPHA8", "EPHB1", "EPHB2", "EPHB3",
    "EPHB4", "EPHB6", "ERBB2", "ERBB3", "ERBB4", "ERK1", "ERK2",
    "FAK", "FER", "FES", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "FGR",
    "FLT1", "FLT3", "FLT4", "FRK", "FYN", "GSK3A", "GSK3B", "GRK1",
    "GRK2", "GRK3", "GRK4", "GRK5", "GRK6", "GRK7", "HCK", "HIPK1",
    "HIPK2", "HIPK3", "HIPK4", "IRAK1", "IRAK2", "IRAK3", "IRAK4",
    "ITK", "JAK1", "JAK2", "JAK3", "TYK2", "KDR", "KIT", "LCK",
    "LIMK1", "LIMK2", "LRRK1", "LRRK2", "LYN", "MAP2K1", "MAP2K2",
    "MAP2K3", "MAP2K4", "MAP2K5", "MAP2K6", "MAP2K7", "MAP3K1", "MAP3K2",
    "MAP3K3", "MAP3K4", "MAP3K5", "MAP3K7", "MAP3K8", "MAP3K11",
    "MAP3K12", "MAP3K13", "MAP3K14", "MAPK1", "MAPK3", "MAPK4", "MAPK6",
    "MAPK7", "MAPK8", "MAPK9", "MAPK10", "MAPK11", "MAPK12", "MAPK13",
    "MAPK14", "MARK1", "MARK2", "MARK3", "MARK4", "MAST1", "MAST2",
    "MAST3", "MAST4", "MATK", "MERTK", "MET", "MINK1", "MKNK1",
    "MKNK2", "MST1R", "MTOR", "NTRK1", "NTRK2", "NTRK3", "PAK1",
    "PAK2", "PAK3", "PAK4", "PAK5", "PAK6", "PDGFRA", "PDGFRB",
    "PDK2", "PDK3", "PDK4", "PHKG1", "PHKG2", "PIK3CA", "PIK3CB",
    "PIK3CD", "PIK3CG", "PIK3C2A", "PIK3C2B", "PIK3C2G", "PIK3R1",
    "PIK3R2", "PIK3R3", "PIM1", "PIM2", "PIM3", "PRKAA1", "PRKAA2",
    "PRKCA", "PRKCB", "PRKCD", "PRKCE", "PRKCG", "PRKCH", "PRKCI",
    "PRKCQ", "PRKCZ", "PRKD1", "PRKD2", "PRKD3",
    # Phosphatases, ubiquitin system (100)
    "PTPN1", "PTPN2", "PTPN3", "PTPN4", "PTPN6", "PTPN7", "PTPN9",
    "PTPN11", "PTPN12", "PTPN13", "PTPN14", "PTPN22", "PTPN23",
    "PTPRA", "PTPRB", "PTPRD", "PTPRE", "PTPRF", "PTPRG", "PTPRH",
    "PTPRJ", "PTPRK", "PTPRM", "PTPRN", "PTPRN2", "PTPRO", "PTPRQ",
    "PTPRR", "PTPRS", "PTPRT", "PTPRU", "DUSP3", "DUSP4", "DUSP5",
    "DUSP6", "DUSP7", "DUSP8", "DUSP9", "DUSP10", "DUSP11", "DUSP12",
    "DUSP14", "DUSP16", "DUSP22", "DUSP26", "DUSP28", "PPP1CA",
    "PPP1CB", "PPP1CC", "PPP2CA", "PPP2CB", "PPP2R1A", "PPP2R1B",
    "PPP2R2A", "PPP2R5A", "PPP3CA", "PPP3CB", "PPP3CC", "PPP4C",
    "PPP5C", "PPP6C", "UBE2A", "UBE2B", "UBE2C", "UBE2D1", "UBE2D2",
    "UBE2D3", "UBE2E1", "UBE2E2", "UBE2E3", "UBE2G1", "UBE2G2",
    "UBE2H", "UBE2I", "UBE2J1", "UBE2J2", "UBE2K", "UBE2L3", "UBE2L6",
    "UBE2M", "UBE2N", "UBE2O", "UBE2Q1", "UBE2R1", "UBE2R2", "UBE2S",
    "UBE2T", "UBE2V1", "UBE2V2", "UBE2W", "UBE2Z", "UBE3A", "UBE3B",
    "UBE3C", "UBE4A", "UBE4B", "USP1", "USP2", "USP3", "USP4", "USP5",
    # Ion channels, transporters (150)
    "KCNA1", "KCNA2", "KCNA3", "KCNA4", "KCNA5", "KCNA6", "KCNA7",
    "KCNB1", "KCNB2", "KCNC1", "KCNC2", "KCNC3", "KCNC4", "KCND1",
    "KCND2", "KCND3", "KCNE1", "KCNE2", "KCNE3", "KCNH1", "KCNH2",
    "KCNJ1", "KCNJ2", "KCNJ3", "KCNJ4", "KCNJ5", "KCNJ6", "KCNJ8",
    "KCNJ10", "KCNJ12", "KCNJ13", "KCNJ14", "KCNJ15", "KCNJ16",
    "KCNK1", "KCNK2", "KCNK3", "KCNK5", "KCNK6", "KCNK9", "KCNK10",
    "KCNMA1", "KCNMB1", "KCNMB2", "KCNN1", "KCNN2", "KCNN3", "KCNN4",
    "KCNQ1", "KCNQ2", "KCNQ3", "KCNQ4", "KCNQ5", "SCN1B", "SCN2A",
    "SCN2B", "SCN3A", "SCN3B", "SCN4A", "SCN4B", "SCN5A", "SCN7A",
    "SCN8A", "SCN9A", "SCN10A", "SCN11A", "CACNA1A", "CACNA1B", "CACNA1C",
    "CACNA1D", "CACNA1E", "CACNA1F", "CACNA1G", "CACNA1H", "CACNA1I",
    "CACNA1S", "CACNA2D1", "CACNA2D2", "CACNA2D3", "CACNB1", "CACNB2",
    "CACNB3", "CACNB4", "CLCN1", "CLCN2", "CLCN3", "CLCN4", "CLCN5",
    "CLCN6", "CLCN7", "CLCNKA", "CLCNKB", "SLC1A1", "SLC1A2", "SLC1A3",
    "SLC1A4", "SLC1A5", "SLC1A6", "SLC1A7", "SLC2A2", "SLC2A4", "SLC2A5",
    "SLC2A6", "SLC2A7", "SLC2A8", "SLC2A9", "SLC2A10", "SLC2A11",
    "SLC2A12", "SLC2A13", "SLC2A14", "SLC3A1", "SLC3A2", "SLC4A1",
    "SLC4A2", "SLC4A3", "SLC4A4", "SLC4A5", "SLC4A7", "SLC4A8",
    "SLC4A9", "SLC4A10", "SLC4A11", "SLC5A1", "SLC5A2", "SLC5A3",
    "SLC5A4", "SLC5A5", "SLC5A6", "SLC5A7", "SLC5A8", "SLC5A9",
    "SLC5A10", "SLC5A11", "SLC5A12", "SLC6A1", "SLC6A2", "SLC6A3",
    "SLC6A4", "SLC6A5", "SLC6A6", "SLC6A7", "SLC6A8", "SLC6A9",
    "SLC6A11", "SLC6A12", "SLC6A13", "SLC6A14", "SLC6A15", "SLC6A17",
    # Metabolic enzymes (200)
    "ACACA", "ACACB", "ACLY", "ACO1", "ACO2", "ACOX1", "ACOX2", "ACOX3",
    "ADA", "ADCY1", "ADCY2", "ADCY3", "ADCY4", "ADCY5", "ADCY6",
    "ADCY7", "ADCY8", "ADCY9", "ADCY10", "ADH1A", "ADH1B", "ADH1C",
    "ADH4", "ADH5", "ADH6", "ADH7", "ALDH1A1", "ALDH1A2", "ALDH1A3",
    "ALDH1B1", "ALDH2", "ALDH3A1", "ALDH3A2", "ALDH3B1", "ALDH3B2",
    "ALDH4A1", "ALDH5A1", "ALDH6A1", "ALDH7A1", "ALDH8A1", "ALDH9A1",
    "ALDH16A1", "ALDH18A1", "ALDOB", "ALDOC", "AMD1", "AMPD1", "AMPD2",
    "AMPD3", "AOC1", "AOC2", "AOC3", "AOX1", "APRT", "ARG1", "ARG2",
    "ASL", "ASS1", "BCAT1", "BCAT2", "BCKDHA", "BCKDHB", "BCKDK",
    "CDA", "CERS1", "CERS2", "CERS3", "CERS4", "CERS5", "CERS6",
    "CKB", "CKM", "CKMT1A", "CKMT1B", "CKMT2", "CPT1A", "CPT1B",
    "CPT1C", "CPT2", "CS", "CYP1A1", "CYP1B1", "CYP2A6", "CYP2B6",
    "CYP2C8", "CYP2C9", "CYP2C18", "CYP2C19", "CYP2D6", "CYP2J2",
    "CYP3A5", "CYP3A7", "CYP4A11", "CYP4A22", "CYP4B1", "CYP4F2",
    "CYP4F3", "CYP4F8", "CYP4F11", "CYP4F12", "CYP4F22", "CYP4V2",
    "CYP4X1", "CYP4Z1", "CYP7A1", "CYP7B1", "CYP8B1", "CYP11A1",
    "CYP11B1", "CYP11B2", "CYP17A1", "CYP19A1", "CYP20A1", "CYP21A2",
    "CYP24A1", "CYP26A1", "CYP26B1", "CYP26C1", "CYP27A1", "CYP27B1",
    "CYP27C1", "CYP39A1", "CYP46A1", "CYP51A1", "DBT", "DCK", "DCTD",
    "DDC", "DHFR", "DHODH", "DLAT", "DLD", "DLST", "DNMT1", "DNMT3A",
    "DNMT3B", "DPYD", "DPYS", "DTYMK", "DUT", "ECHS1", "ECI1", "ECI2",
    "EHHADH", "ENPP1", "ENPP2", "ENPP3", "EPHX1", "EPHX2", "ETFA",
    "ETFB", "ETFDH", "FABP1", "FABP2", "FABP3", "FABP4", "FABP5",
    "FABP6", "FABP7", "FADS1", "FADS2", "FADS3", "FAH", "FASN",
    "FDFT1", "FDPS", "FH", "GALC", "GALNT1", "GALNT2", "GALNT3",
    "GALNT4", "GALNT5", "GALNT6", "GALNT7", "GALT", "GAMT", "GAA",
    "GBA", "GCK", "GCKR", "GCSH", "GGT1", "GLDC", "GLDH", "GLUD1",
    "GLUD2", "GLUL", "GLS", "GLS2", "GMPS", "GNE", "GOT1", "GOT2",
    "GPD1", "GPD1L", "GPD2", "GPI", "GPT", "GPT2", "GPX1", "GPX2",
    "GPX3", "GPX4", "GPX5", "GPX6", "GPX7", "GPX8", "GSR", "GSS",
    # Structural / cytoskeleton (150)
    "ACTC1", "ACTG1", "ACTG2", "ACTN1", "ACTN2", "ACTN3", "ACTN4",
    "ADD1", "ADD2", "ADD3", "ANK1", "ANK2", "ANK3", "ARPC1A", "ARPC1B",
    "ARPC2", "ARPC3", "ARPC4", "ARPC5", "ARPC5L", "CFL1", "CFL2",
    "CORO1A", "CORO1B", "CORO1C", "CORO2A", "CORO2B", "COTL1", "CSRP1",
    "CSRP2", "CSRP3", "DES", "DMD", "DNM1", "DNM1L", "DNM2", "DNM3",
    "DSP", "DSC1", "DSC2", "DSC3", "DSG1", "DSG2", "DSG3", "DSG4",
    "DTNBP1", "DTNA", "DTNB", "EPB41", "EPB41L1", "EPB41L2", "EPB41L3",
    "EPB41L4A", "EPB41L4B", "EPB41L5", "EZR", "FHOD1", "FHOD3",
    "FKBP10", "FKBP14", "FLNA", "FLNB", "FLNC", "FSCN1", "FSCN2",
    "FSCN3", "GJA1", "GJA3", "GJA4", "GJA5", "GJA8", "GJB1", "GJB2",
    "GJB3", "GJB4", "GJB5", "GJB6", "GJC1", "GJC2", "GJC3", "GJD2",
    "GJD3", "GJD4", "GSN", "HSPB1", "HSPB2", "HSPB3", "HSPB5",
    "HSPB6", "HSPB7", "HSPB8", "HSPB9", "INA", "INF2", "KIF1A",
    "KIF1B", "KIF2A", "KIF2B", "KIF2C", "KIF3A", "KIF3B", "KIF3C",
    "KIF4A", "KIF4B", "KIF5A", "KIF5B", "KIF5C", "KIF7", "KIF9",
    "KIF14", "KIF15", "KIF18A", "KIF18B", "KIF20A", "KIF20B", "KIF21A",
    "KIF21B", "KIF22", "KIF23", "KIF24", "KIF25", "KIF26A", "KIF26B",
    "KIF27", "LMNA", "LMNB2", "MYH1", "MYH2", "MYH3", "MYH4",
    "MYH6", "MYH7", "MYH7B", "MYH8", "MYH9", "MYH10", "MYH11",
    "MYH13", "MYH14", "MYH15", "MYH16", "MYL1", "MYL2", "MYL3",
    "MYL4", "MYL5", "MYL6", "MYL6B", "MYL7", "MYL9", "MYL12A",
    "MYL12B", "MYO1A", "MYO1B", "MYO1C",
    # Receptors / signaling (200)
    "ACVR1", "ACVR1B", "ACVR1C", "ACVR2A", "ACVR2B", "ACVRL1",
    "ADIPOQ", "ADIPOR1", "ADIPOR2", "ADORA1", "ADORA2A", "ADORA2B",
    "ADORA3", "ADRB1", "ADRB2", "ADRB3", "ADRA1A", "ADRA1B", "ADRA1D",
    "ADRA2A", "ADRA2B", "ADRA2C", "AGTR1", "AGTR2", "AVPR1A", "AVPR1B",
    "AVPR2", "BMPR1A", "BMPR1B", "BMPR2", "CALCR", "CALCRL", "CASR",
    "CCR1", "CCR2", "CCR3", "CCR4", "CCR5", "CCR6", "CCR7", "CCR8",
    "CCR9", "CCR10", "CXCR1", "CXCR2", "CXCR3", "CXCR4", "CXCR5",
    "CXCR6", "CX3CR1", "XCR1", "CHRM1", "CHRM2", "CHRM3", "CHRM4",
    "CHRM5", "CHRNA1", "CHRNA2", "CHRNA3", "CHRNA4", "CHRNA5", "CHRNA6",
    "CHRNA9", "CHRNA10", "CHRNB1", "CHRNB2", "CHRNB3", "CHRNB4",
    "CHRND", "CHRNE", "CHRNG", "CTLA4", "DCC", "DRD1", "DRD3", "DRD4",
    "DRD5", "EDNRA", "EDNRB", "F2R", "F2RL1", "F2RL2", "F2RL3",
    "FGFRL1", "FPR1", "FPR2", "FPR3", "FSHR", "GHR", "GHRH", "GHRHR",
    "GLP1R", "GLP2R", "GNRH1", "GNRHR", "GPER1", "GRM1", "GRM2",
    "GRM3", "GRM4", "GRM5", "GRM6", "GRM7", "GRM8", "HCRTR1", "HCRTR2",
    "HTR1A", "HTR1B", "HTR1D", "HTR1E", "HTR1F", "HTR2B", "HTR2C",
    "HTR3A", "HTR3B", "HTR3C", "HTR3D", "HTR3E", "HTR4", "HTR5A",
    "HTR6", "HTR7", "IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2", "IFNLR1",
    "IGF1R", "IGF2R", "IL1R1", "IL1R2", "IL2RA", "IL2RB", "IL2RG",
    "IL3RA", "IL4R", "IL5RA", "IL6R", "IL6ST", "IL7R", "IL10RA",
    "IL10RB", "IL11RA", "IL12RB1", "IL12RB2", "IL13RA1", "IL13RA2",
    "IL15RA", "IL17RA", "IL17RB", "IL17RC", "IL17RD", "IL17RE",
    "IL18R1", "IL18RAP", "IL20RA", "IL20RB", "IL21R", "IL22RA1",
    "IL22RA2", "IL23R", "IL27RA", "IL31RA", "INSR", "IRS1", "IRS2",
    "IRS4", "KISS1R", "LEPR", "LHCGR", "LIFR", "LRP1", "LRP2",
    "LRP4", "LRP5", "LRP6", "LRP8", "MC1R", "MC2R", "MC3R", "MC4R",
    "MC5R", "NPR1", "NPR2", "NPR3", "NPSR1", "NPY1R", "NPY2R",
    "NPY4R", "NPY5R", "OPRD1", "OPRK1", "OPRL1", "OPRM1", "OSMR",
    "OXTR", "P2RX1", "P2RX2", "P2RX3", "P2RX4", "P2RX5", "P2RX6",
    "P2RX7", "P2RY1", "P2RY2",
    # Chromatin / epigenetics (150)
    "ARID1A", "ARID1B", "ARID2", "ARID3A", "ARID3B", "ARID3C", "ARID4A",
    "ARID4B", "ARID5A", "ARID5B", "ASH1L", "ASH2L", "ASXL1", "ASXL2",
    "ASXL3", "ATM", "ATR", "ATRX", "BAP1", "BRD1", "BRD2", "BRD3",
    "BRD4", "BRD7", "BRD8", "BRD9", "BRDT", "BRPF1", "BRPF3",
    "CBX1", "CBX2", "CBX3", "CBX4", "CBX5", "CBX6", "CBX7", "CBX8",
    "CECR2", "CHAMP1", "CHD1", "CHD2", "CHD3", "CHD4", "CHD5", "CHD6",
    "CHD7", "CHD8", "CHD9", "CREBBP", "CTCF", "DAXX", "DOT1L",
    "EED", "EHMT1", "EHMT2", "EP300", "EP400", "EZH1", "EZH2",
    "HDAC1", "HDAC2", "HDAC3", "HDAC4", "HDAC5", "HDAC6", "HDAC7",
    "HDAC8", "HDAC9", "HDAC10", "HDAC11", "HIRA", "HMGA1B",
    "HP1BP3", "INO80", "ING1", "ING2", "ING3", "ING4", "ING5",
    "JARID2", "KAT2A", "KAT2B", "KAT5", "KAT6A", "KAT6B", "KAT7",
    "KAT8", "KDM1A", "KDM1B", "KDM2A", "KDM2B", "KDM3A", "KDM3B",
    "KDM4A", "KDM4B", "KDM4C", "KDM4D", "KDM5A", "KDM5B", "KDM5C",
    "KDM5D", "KDM6A", "KDM6B", "KDM7A", "KDM8", "KMT2A", "KMT2B",
    "KMT2C", "KMT2D", "KMT2E", "L3MBTL1", "L3MBTL2", "L3MBTL3",
    "L3MBTL4", "MBD1", "MBD2", "MBD3", "MBD4", "MECP2", "MORF4L1",
    "MORF4L2", "MSL1", "MSL2", "MSL3", "NAP1L1", "NAP1L2", "NAP1L3",
    "NAP1L4", "NAP1L5", "NCOA1", "NCOA2", "NCOA3", "NCOR1", "NCOR2",
    "NSD1", "NSD2", "NSD3", "PBRM1", "PHF1", "PHF2", "PHF3", "PHF5A",
    "PHF6", "PHF7", "PHF8", "PHF10", "PHF12", "PHF13", "PHF14",
    "PHF19", "PHF20", "PHF20L1", "PHF21A", "PHF21B",
    # RNA processing / splicing (100)
    "CELF1", "CELF2", "CELF3", "CELF4", "CELF5", "CELF6", "CPSF1",
    "CPSF2", "CPSF3", "CPSF4", "CPSF6", "CPSF7", "CSTF1", "CSTF2",
    "CSTF2T", "CSTF3", "DDX1", "DDX3X", "DDX3Y", "DDX5", "DDX6",
    "DDX17", "DDX20", "DDX21", "DDX23", "DDX24", "DDX27", "DDX31",
    "DDX39A", "DDX39B", "DDX41", "DDX42", "DDX46", "DDX47", "DDX49",
    "DDX50", "DDX51", "DDX52", "DDX54", "DDX55", "DDX56", "DHX15",
    "DHX16", "DHX29", "DHX30", "DHX33", "DHX34", "DHX36", "DHX37",
    "DHX38", "DHX40", "DHX57", "DHX58", "DICER1", "DROSHA", "DGCR8",
    "EIF4A1", "EIF4A2", "EIF4A3", "ELAVL1", "ELAVL2", "ELAVL3",
    "ELAVL4", "HNRNPA1", "HNRNPA2B1", "HNRNPA3", "HNRNPAB", "HNRNPC",
    "HNRNPD", "HNRNPDL", "HNRNPF", "HNRNPH1", "HNRNPH2", "HNRNPH3",
    "HNRNPK", "HNRNPL", "HNRNPLL", "HNRNPM", "HNRNPR", "HNRNPU",
    "HNRNPUL1", "HNRNPUL2", "KHDRBS1", "KHDRBS2", "KHDRBS3", "KHSRP",
    "MBNL1", "MBNL2", "MBNL3", "NOVA1", "NOVA2", "PABPC1", "PABPC3",
    "PABPC4", "PABPN1", "PRPF3", "PRPF4", "PRPF6", "PRPF8", "PRPF19",
    "PRPF31", "PRPF38A", "PRPF38B", "PRPF39", "PRPF40A", "PRPF40B",
    # Apoptosis / cell death (100)
    "BAD", "BAK1", "BBC3", "BCL2", "BCL2A1", "BCL2L2",
    "BCL2L10", "BCL2L11", "BCL2L12", "BCL2L13", "BCL2L14", "BCL2L15",
    "BID", "BIK", "BIRC2", "BIRC5", "BIRC6", "BIRC7", "BMF",
    "BOK", "CASP1", "CASP2", "CASP3", "CASP4", "CASP5", "CASP6",
    "CASP7", "CASP8", "CASP9", "CASP10", "CASP12", "CASP14", "CFLAR",
    "CYCS", "DAPK1", "DFFA", "DFFB", "DIABLO", "ENDOG", "FADD",
    "FAS", "FASLG", "HTRA2", "IAP", "RIPK1", "RIPK2", "RIPK3",
    "RIPK4", "TNF", "TNFRSF1A", "TNFRSF1B", "TNFRSF4", "TNFRSF6B",
    "TNFRSF8", "TNFRSF9", "TNFRSF10A", "TNFRSF10B", "TNFRSF10C",
    "TNFRSF10D", "TNFRSF11A", "TNFRSF11B", "TNFRSF12A", "TNFRSF13B",
    "TNFRSF13C", "TNFRSF14", "TNFRSF17", "TNFRSF18", "TNFRSF19",
    "TNFRSF21", "TNFRSF25", "TNFSF4", "TNFSF6", "TNFSF8", "TNFSF9",
    "TNFSF10", "TNFSF11", "TNFSF12", "TNFSF13", "TNFSF13B", "TNFSF14",
    "TNFSF15", "TNFSF18", "TRADD", "TRAF2", "TRAF3", "TRAF4", "TRAF5",
    "TRAF6", "XIAP", "APAF1", "AIFM1", "AIFM2", "PMAIP1", "HRK",
    "MOAP1", "NLRP1", "NLRP3", "PYCARD", "GSDMD", "GSDME", "MLKL",
    # DNA repair (100)
    "APEX1", "APEX2", "BRCA1", "BRCA2", "CETN2", "DDB1", "DDB2",
    "ERCC1", "ERCC2", "ERCC3", "ERCC4", "ERCC5", "ERCC6", "ERCC8",
    "EXO1", "FANCA", "FANCB", "FANCC", "FANCD2", "FANCE", "FANCF",
    "FANCG", "FANCI", "FANCL", "FANCM", "LIG3", "LIG4", "MBD4",
    "MLH1", "MLH3", "MMS19", "MPG", "MRE11", "MSH2", "MSH3", "MSH5",
    "MSH6", "MUTYH", "NBN", "NEIL1", "NEIL2", "NEIL3", "NHEJ1",
    "NTHL1", "OGG1", "PARP1", "PARP2", "PARP3", "PARP4", "PMS1",
    "PMS2", "PNKP", "POLB", "POLD2", "POLD3", "POLD4", "POLE",
    "POLE2", "POLE3", "POLE4", "POLG", "POLG2", "POLH", "POLI",
    "POLK", "POLL", "POLM", "POLN", "POLQ", "RAD50", "RAD51",
    "RAD51B", "RAD51C", "RAD51D", "RAD52", "RAD54B", "RAD54L",
    "RNF8", "RNF168", "RBBP8", "REV1", "REV3L", "RPA1", "RPA2",
    "RPA3", "SMARCAL1", "TDP1", "TDP2", "TERT", "TERF1", "TERF2",
    "TINF2", "TOP1", "TOP3A", "TOP3B", "TP53BP1", "UNG", "WRN",
    "XPA", "XPC", "XRCC1", "XRCC2", "XRCC3", "XRCC4", "XRCC5",
    "XRCC6",
    # Membrane trafficking / autophagy (100)
    "AP1B1", "AP1G1", "AP1M1", "AP1S1", "AP2A1", "AP2A2", "AP2B1",
    "AP2M1", "AP2S1", "AP3B1", "AP3B2", "AP3D1", "AP3M1", "AP3M2",
    "AP3S1", "AP3S2", "AP4B1", "AP4E1", "AP4M1", "AP4S1", "AP5B1",
    "AP5M1", "AP5S1", "AP5Z1", "ATG2A", "ATG2B", "ATG3", "ATG4A",
    "ATG4B", "ATG4C", "ATG4D", "ATG5", "ATG7", "ATG9A", "ATG9B",
    "ATG10", "ATG12", "ATG13", "ATG14", "ATG16L1", "ATG16L2",
    "BECN1", "BNIP3", "CALCOCO2", "GABARAP", "GABARAPL1", "GABARAPL2",
    "LAMP1", "LAMP2", "LAMP3", "MAP1LC3A", "MAP1LC3B", "MAP1LC3C",
    "NBR1", "OPTN", "PIK3C3", "RAB1A", "RAB1B", "RAB2A", "RAB2B",
    "RAB3A", "RAB3B", "RAB3C", "RAB3D", "RAB4A", "RAB4B", "RAB5A",
    "RAB5B", "RAB5C", "RAB6A", "RAB6B", "RAB7A", "RAB7B", "RAB8A",
    "RAB8B", "RAB9A", "RAB9B", "RAB10", "RAB11A", "RAB11B", "RAB13",
    "RAB14", "RAB15", "RAB18", "RAB20", "RAB21", "RAB22A", "RAB23",
    "RAB24", "RAB25", "RAB27A", "RAB27B", "RAB28", "RAB29", "RAB30",
    "RAB31", "RAB32", "RAB33A", "RAB33B", "RAB34",
    # Secreted / ECM / growth factors (150)
    "ANGPT1", "ANGPT2", "ANGPT4", "BMP1", "BMP2", "BMP3", "BMP4",
    "BMP5", "BMP6", "BMP7", "BMP8A", "BMP8B", "BMP10", "BMP15",
    "CCL1", "CCL5", "CCL7", "CCL8", "CCL11", "CCL13", "CCL14",
    "CCL15", "CCL16", "CCL17", "CCL18", "CCL19", "CCL21", "CCL22",
    "CCL23", "CCL24", "CCL25", "CCL26", "CCL27", "CCL28",
    "COL2A1", "COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5",
    "COL4A6", "COL5A2", "COL5A3", "COL6A1", "COL6A2", "COL6A3",
    "COL7A1", "COL8A1", "COL8A2", "COL9A1", "COL9A2", "COL9A3",
    "COL10A1", "COL11A1", "COL11A2", "COL12A1", "COL13A1", "COL14A1",
    "COL15A1", "COL16A1", "COL17A1", "COL18A1", "COL19A1", "COL20A1",
    "COL21A1", "COL22A1", "COL23A1", "COL24A1", "COL25A1", "COL26A1",
    "COL27A1", "COL28A1", "CXCL4", "CXCL5", "CXCL6", "CXCL7",
    "CXCL11", "CXCL12", "CXCL13", "CXCL14", "CXCL16", "CXCL17",
    "EGF", "FGF1", "FGF2", "FGF3", "FGF4", "FGF5", "FGF6", "FGF7",
    "FGF8", "FGF9", "FGF10", "FGF11", "FGF12", "FGF13", "FGF14",
    "FGF16", "FGF17", "FGF18", "FGF19", "FGF20", "FGF21", "FGF22",
    "FGF23", "GDF1", "GDF2", "GDF3", "GDF5", "GDF6", "GDF7", "GDF9",
    "GDF10", "GDF11", "GDF15", "HGF", "IGF1", "IGF2", "IL2", "IL3",
    "IL4", "IL5", "IL7", "IL9", "IL10", "IL11", "IL12A", "IL12B",
    "IL13", "IL15", "IL17A", "IL17B", "IL17C", "IL17D", "IL17F",
    "IL18", "IL19", "IL20", "IL21", "IL22", "IL23A", "IL24", "IL25",
    "IL26", "IL27", "IL28A", "IL28B", "IL29", "IL31", "IL33", "IL34",
    "IL36A", "IL36B", "IL36G", "IL36RN", "IL37", "LIF", "OSM",
    "PDGFA", "PDGFB", "PDGFC", "PDGFD", "PGF",
    # Cell cycle regulation extras (50)
    "CCNC", "CCND1", "CCND2", "CCND3", "CCNE2", "CCNF", "CCNG1",
    "CCNG2", "CCNH", "CCNI", "CCNJ", "CCNK", "CCNL1", "CCNL2",
    "CCNO", "CCNT1", "CCNT2", "CCNY", "CDC6", "CDC7", "CDC14A",
    "CDC14B", "CDC16", "CDC23", "CDC25A", "CDC25B", "CDC25C", "CDC27",
    "CDC34", "CDC42", "CDC45", "CDC73", "CDCA2", "CDCA3", "CDCA4",
    "CDCA5", "CDCA7", "CDCA7L", "CDCA8", "CDT1", "CENPB", "CENPC",
    "CENPD", "CENPE", "CENPF", "CENPH", "CENPI", "CENPJ", "CENPK",
    "CENPL",
    # Proteases / peptidases (100)
    "ADAM9", "ADAM10", "ADAM12", "ADAM15", "ADAM17", "ADAM19", "ADAM22",
    "ADAM23", "ADAM28", "ADAM33", "ADAMTS1", "ADAMTS2", "ADAMTS3",
    "ADAMTS4", "ADAMTS5", "ADAMTS6", "ADAMTS7", "ADAMTS8", "ADAMTS9",
    "ADAMTS10", "ADAMTS12", "ADAMTS13", "ADAMTS14", "ADAMTS15",
    "ADAMTS16", "ADAMTS17", "ADAMTS18", "ADAMTS19", "ADAMTS20",
    "BACE1", "BACE2", "CAPN1", "CAPN2", "CAPN3", "CAPN5", "CAPN6",
    "CAPN7", "CAPN8", "CAPN9", "CAPN10", "CAPN11", "CAPN12", "CAPN13",
    "CAPN14", "CAPN15", "CASP8AP2", "CTSA", "CTSB", "CTSC", "CTSD",
    "CTSE", "CTSF", "CTSG", "CTSH", "CTSK", "CTSL", "CTSO", "CTSS",
    "CTSV", "CTSW", "CTSZ", "DPP3", "DPP4", "DPP6", "DPP7", "DPP8",
    "DPP9", "DPP10", "FURIN", "HTRA1", "HTRA3", "KLK1", "KLK2",
    "KLK3", "KLK4", "KLK5", "KLK6", "KLK7", "KLK8", "KLK9", "KLK10",
    "KLK11", "KLK12", "KLK13", "KLK14", "KLK15", "MMP1", "MMP7",
    "MMP8", "MMP10", "MMP11", "MMP12", "MMP13", "MMP14", "MMP15",
    "MMP16", "MMP17", "MMP19", "MMP20", "MMP21", "MMP23B", "MMP24",
    "MMP25", "MMP26", "MMP27", "MMP28",
    # Misc well-known genes (200)
    "ABCB1", "ABCB4", "ABCB11", "ABCC1", "ABCC2", "ABCC3", "ABCC4",
    "ABCC5", "ABCC6", "ABCC8", "ABCC9", "ABCC10", "ABCC11", "ABCC12",
    "ABCG1", "ABCG2", "ABCG4", "ABCG5", "ABCG8", "ABI1", "ABI2",
    "ABI3", "ABL1", "ACKR1", "ACKR2", "ACKR3", "ACKR4", "ACP1",
    "ACP2", "ACP5", "ACP6", "ACPP", "ACSL1", "ACSL3", "ACSL4",
    "ACSL5", "ACSL6", "ACSM1", "ACSM2A", "ACSM2B", "ACSM3", "ACSM5",
    "ACSS1", "ACSS2", "ACSS3", "ADAR", "ADARB1", "ADARB2", "ADGRB1",
    "ADGRB2", "ADGRB3", "ADGRE1", "ADGRE2", "ADGRE3", "ADGRE5",
    "ADGRF1", "ADGRF2", "ADGRF3", "ADGRF4", "ADGRF5", "ADGRG1",
    "ADGRG2", "ADGRG3", "ADGRG4", "ADGRG5", "ADGRG6", "ADGRG7",
    "ADGRL1", "ADGRL2", "ADGRL3", "ADGRL4", "ADGRV1", "AFP",
    "AGER", "AGO1", "AGO2", "AGO3", "AGO4", "AHNAK", "AHNAK2",
    "AHSG", "AIF1", "AIM2", "AIRE", "AIMP1", "AIMP2", "AKAP1",
    "AKAP2", "AKAP5", "AKAP6", "AKAP7", "AKAP8", "AKAP9", "AKAP10",
    "AKAP11", "AKAP12", "AKAP13", "AKR1A1", "AKR1B1", "AKR1B10",
    "AKR1C1", "AKR1C2", "AKR1C3", "AKR1C4", "AKR1D1", "AKR7A2",
    "AKR7A3", "ALAD", "ALAS1", "ALAS2", "ALCAM", "ALDH1L1", "ALDH1L2",
    "ALOX5", "ALOX12", "ALOX15", "ALOX15B", "ALPI", "ALPL", "ANG",
    "ANGPTL1", "ANGPTL2", "ANGPTL3", "ANGPTL5", "ANGPTL6", "ANGPTL7",
    "ANGPTL8", "ANO1", "ANO2", "ANO3", "ANO4", "ANO5", "ANO6",
    "ANO7", "ANO8", "ANO9", "ANO10", "ANPEP", "ANXA1", "ANXA2",
    "ANXA3", "ANXA4", "ANXA5", "ANXA6", "ANXA7", "ANXA8", "ANXA9",
    "ANXA10", "ANXA11", "ANXA13", "AP1AR", "APBA1", "APBA2", "APBA3",
    "APBB1", "APBB2", "APBB3", "APC", "APC2", "APCS", "APEH",
    "APEX1", "APH1A", "APH1B", "API5", "APOA2", "APOA4", "APOA5",
    "APOC1", "APOC2", "APOC4", "APOD", "APOE", "APOF", "APOH",
    "APOL1", "APOL2", "APOL3", "APOL4", "APOL5", "APOL6", "APP",
    "AQP1", "AQP2", "AQP3", "AQP4", "AQP5", "AQP7", "AQP8",
    "AQP9", "AQP10", "AQP11", "AQP12A", "AQP12B",
    # Additional common HGNC symbols spanning diverse biology (3000+)
    "NEDD4", "NEDD4L", "NEDD8", "NEDD9", "SMURF1", "SMURF2", "WWP1",
    "WWP2", "ITCH", "HECTD1", "HECTD2", "HECTD3", "HECTD4", "HERC1",
    "HERC2", "HERC3", "HERC4", "HERC6", "HUWE1", "UBR1", "UBR2",
    "UBR3", "UBR4", "UBR5", "UBR7", "TRIP12", "TRIM2", "TRIM3",
    "TRIM5", "TRIM6", "TRIM7", "TRIM8", "TRIM9", "TRIM10", "TRIM11",
    "TRIM13", "TRIM14", "TRIM15", "TRIM16", "TRIM17", "TRIM21",
    "TRIM23", "TRIM24", "TRIM25", "TRIM26", "TRIM27",
    "TRIM28", "TRIM29", "TRIM31", "TRIM32", "TRIM33", "TRIM34",
    "TRIM35", "TRIM36", "TRIM37", "TRIM38", "TRIM39", "TRIM40",
    "TRIM41", "TRIM43", "TRIM44", "TRIM45", "TRIM46", "TRIM47",
    "TRIM48", "TRIM49", "TRIM50", "TRIM51", "TRIM52", "TRIM54",
    "TRIM55", "TRIM56", "TRIM58", "TRIM59", "TRIM60", "TRIM61",
    "TRIM62", "TRIM63", "TRIM65", "TRIM66", "TRIM67", "TRIM68",
    "TRIM69", "TRIM71", "TRIM72", "TRIM73", "TRIM74",
    "RNF2", "RNF4", "RNF5", "RNF6", "RNF7", "RNF10", "RNF11",
    "RNF12", "RNF13", "RNF14", "RNF17", "RNF19A", "RNF19B", "RNF20",
    "RNF25", "RNF31", "RNF34", "RNF38", "RNF40", "RNF41", "RNF43",
    "RNF111", "RNF114", "RNF115", "RNF121", "RNF125", "RNF126",
    "RNF128", "RNF130", "RNF133", "RNF138", "RNF139", "RNF141",
    "RNF144A", "RNF144B", "RNF145", "RNF146", "RNF148", "RNF149",
    "RNF150", "RNF151", "RNF152", "RNF157", "RNF165", "RNF166",
    "RNF167", "RNF170", "RNF180", "RNF181", "RNF182", "RNF183",
    "RNF185", "RNF186", "RNF187", "ZNRF1", "ZNRF2", "ZNRF3", "ZNRF4",
    "STUB1", "CHIP", "CUL1", "CUL2", "CUL3", "CUL4A", "CUL4B",
    "CUL5", "CUL7", "CUL9", "SKP1", "SKP2", "FBXO1", "FBXO2",
    "FBXO3", "FBXO4", "FBXO5", "FBXO6", "FBXO7", "FBXO8", "FBXO9",
    "FBXO10", "FBXO11", "FBXO15", "FBXO17", "FBXO18", "FBXO21",
    "FBXO22", "FBXO24", "FBXO25", "FBXO27", "FBXO28", "FBXO30",
    "FBXO31", "FBXO32", "FBXO33", "FBXO34", "FBXO36", "FBXO38",
    "FBXO39", "FBXO40", "FBXO41", "FBXO42", "FBXO43", "FBXO44",
    "FBXO45", "FBXO46", "FBXO47", "FBXO48", "FBXL2", "FBXL3",
    "FBXL4", "FBXL5", "FBXL6", "FBXL7", "FBXL8", "FBXL12", "FBXL13",
    "FBXL14", "FBXL15", "FBXL16", "FBXL17", "FBXL18", "FBXL19",
    "FBXL20", "FBXL21", "FBXL22", "FBXW2", "FBXW4", "FBXW5",
    "FBXW7", "FBXW8", "FBXW9", "FBXW10", "FBXW11", "FBXW12",
    "BTRC", "BTBD1", "BTBD2", "BTBD3", "BTBD6", "BTBD7", "BTBD8",
    "BTBD9", "BTBD10", "BTBD11", "BTBD16", "BTBD17", "KBTBD2",
    "KBTBD3", "KBTBD4", "KBTBD6", "KBTBD7", "KBTBD8", "KBTBD11",
    "KBTBD12", "KBTBD13",
    # Zinc finger proteins
    "ZNF1", "ZNF2", "ZNF3", "ZNF10", "ZNF12", "ZNF14", "ZNF16",
    "ZNF17", "ZNF18", "ZNF19", "ZNF20", "ZNF22", "ZNF23", "ZNF24",
    "ZNF25", "ZNF26", "ZNF28", "ZNF30", "ZNF32", "ZNF33A", "ZNF33B",
    "ZNF34", "ZNF35", "ZNF36", "ZNF37A", "ZNF37B", "ZNF38", "ZNF41",
    "ZNF43", "ZNF44", "ZNF45", "ZNF46", "ZNF48", "ZNF57", "ZNF69",
    "ZNF70", "ZNF71", "ZNF74", "ZNF75A", "ZNF75D", "ZNF76", "ZNF77",
    "ZNF79", "ZNF80", "ZNF81", "ZNF83", "ZNF84", "ZNF85", "ZNF90",
    "ZNF91", "ZNF92", "ZNF93", "ZNF98", "ZNF99", "ZNF100", "ZNF101",
    "ZNF107", "ZNF112", "ZNF114", "ZNF117", "ZNF121", "ZNF124",
    "ZNF131", "ZNF132", "ZNF133", "ZNF134", "ZNF135", "ZNF136",
    "ZNF137P", "ZNF138", "ZNF140", "ZNF141", "ZNF142", "ZNF143",
    "ZNF146", "ZNF148", "ZNF154", "ZNF155", "ZNF157", "ZNF160",
    "ZNF165", "ZNF169", "ZNF174", "ZNF175", "ZNF177", "ZNF180",
    "ZNF181", "ZNF182", "ZNF184", "ZNF185", "ZNF189", "ZNF195",
    "ZNF197", "ZNF200", "ZNF202", "ZNF205", "ZNF207", "ZNF208",
    "ZNF211", "ZNF212", "ZNF213", "ZNF214", "ZNF215", "ZNF217",
    "ZNF219", "ZNF220", "ZNF221", "ZNF222", "ZNF223", "ZNF224",
    "ZNF225", "ZNF226", "ZNF227", "ZNF228", "ZNF229", "ZNF230",
    "ZNF232", "ZNF233", "ZNF234", "ZNF235", "ZNF236", "ZNF239",
    "ZNF248", "ZNF250", "ZNF251", "ZNF253", "ZNF254",
    # Solute carriers (additional)
    "SLC7A1", "SLC7A2", "SLC7A3", "SLC7A4", "SLC7A5", "SLC7A6",
    "SLC7A7", "SLC7A8", "SLC7A9", "SLC7A10", "SLC7A11", "SLC7A13",
    "SLC7A14", "SLC8A1", "SLC8A2", "SLC8A3", "SLC9A1", "SLC9A2",
    "SLC9A3", "SLC9A4", "SLC9A5", "SLC9A6", "SLC9A7", "SLC9A8",
    "SLC9A9", "SLC10A1", "SLC10A2", "SLC10A3", "SLC10A4", "SLC10A5",
    "SLC10A6", "SLC10A7", "SLC11A1", "SLC11A2", "SLC12A1", "SLC12A2",
    "SLC12A3", "SLC12A4", "SLC12A5", "SLC12A6", "SLC12A7", "SLC12A8",
    "SLC12A9", "SLC13A1", "SLC13A2", "SLC13A3", "SLC13A4", "SLC13A5",
    "SLC14A1", "SLC14A2", "SLC15A1", "SLC15A2", "SLC15A3", "SLC15A4",
    "SLC16A1", "SLC16A2", "SLC16A3", "SLC16A4", "SLC16A5", "SLC16A6",
    "SLC16A7", "SLC16A8", "SLC16A9", "SLC16A10", "SLC16A11",
    "SLC16A12", "SLC16A13", "SLC16A14", "SLC17A1", "SLC17A2",
    "SLC17A3", "SLC17A4", "SLC17A5", "SLC17A6", "SLC17A8",
    "SLC17A9", "SLC18A1", "SLC18A2", "SLC18A3", "SLC18B1",
    "SLC19A1", "SLC19A2", "SLC19A3", "SLC20A1", "SLC20A2",
    "SLC22A1", "SLC22A2", "SLC22A3", "SLC22A4", "SLC22A5",
    "SLC22A6", "SLC22A7", "SLC22A8", "SLC22A9", "SLC22A10",
    "SLC22A11", "SLC22A12", "SLC22A13", "SLC22A14", "SLC22A15",
    "SLC22A16", "SLC22A17", "SLC22A18", "SLC22A23", "SLC22A24",
    "SLC22A25", "SLC22A31",
    # G-protein signaling
    "GNAI1", "GNAI2", "GNAI3", "GNAO1", "GNAQ", "GNAS", "GNAL",
    "GNA11", "GNA12", "GNA13", "GNA14", "GNA15", "GNAZ", "GNB1",
    "GNB2", "GNB3", "GNB4", "GNB5", "GNG2", "GNG3", "GNG4", "GNG5",
    "GNG7", "GNG8", "GNG10", "GNG11", "GNG12", "GNG13", "GNGT1",
    "GNGT2", "RGS1", "RGS2", "RGS3", "RGS4", "RGS5", "RGS6",
    "RGS7", "RGS8", "RGS9", "RGS10", "RGS11", "RGS12", "RGS13",
    "RGS14", "RGS16", "RGS17", "RGS18", "RGS19", "RGS20", "RGS21",
    "RGS22",
    # WNT signaling
    "WNT1", "WNT2", "WNT2B", "WNT3", "WNT3A", "WNT4", "WNT5A",
    "WNT5B", "WNT6", "WNT7A", "WNT7B", "WNT8A", "WNT8B", "WNT9A",
    "WNT9B", "WNT10A", "WNT10B", "WNT11", "WNT16", "FZD1", "FZD2",
    "FZD3", "FZD4", "FZD5", "FZD6", "FZD7", "FZD8", "FZD9", "FZD10",
    "SFRP1", "SFRP2", "SFRP4", "SFRP5", "DKK1", "DKK2", "DKK3",
    "DKK4", "DKKL1", "AXIN1", "AXIN2", "DVL1", "DVL2", "DVL3",
    "CTNNB1", "CSNK1A1", "CSNK1D", "CSNK1E", "CSNK1G1", "CSNK1G2",
    "CSNK1G3", "CSNK2A1", "CSNK2A2", "CSNK2B",
    # Hedgehog pathway
    "SHH", "IHH", "DHH", "PTCH1", "PTCH2", "SMO", "GLI1", "GLI2",
    "GLI3", "SUFU", "STK36", "HHIP", "BOC", "CDON", "GAS1",
    # Hippo pathway
    "MST1", "MST2", "SAV1", "LATS1", "LATS2", "MOB1A", "MOB1B",
    "YAP1", "WWTR1", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "NF2",
    "AMOT", "AMOTL1", "AMOTL2", "RASSF1", "RASSF2", "RASSF3",
    "RASSF4", "RASSF5", "RASSF6",
    # TGF-beta / BMP pathway extras
    "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TGFBR3",
    "INHBA", "INHBB", "INHBC", "INHBE", "INHA", "LEFTY1", "LEFTY2",
    "NODAL", "GDF1", "MSTN", "CER1", "CHRD", "CHRDL1", "CHRDL2",
    "NOG", "FST", "FSTL1", "FSTL3", "FSTL4", "FSTL5", "LTBP1",
    "LTBP2", "LTBP3", "LTBP4",
    # Misc signaling (additional)
    "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6",
    "SRC", "YES1", "FGR", "RAF1", "ARAF", "CRAF",
    "HRAS", "KRAS", "NRAS", "RRAS", "RRAS2", "MRAS", "RALA", "RALB",
    "RAP1A", "RAP1B", "RAP2A", "RAP2B", "RAP2C", "RHOA", "RHOB",
    "RHOC", "RHOD", "RHOF", "RHOG", "RHOH", "RHOQ", "RHOU", "RHOV",
    "RAC1", "RAC2", "RAC3", "RND1", "RND2", "RND3",
    "PTEN", "TSC1", "TSC2", "RPTOR", "RICTOR", "MLST8", "DEPTOR",
    "STK11", "STRADA", "STRADB", "CAB39", "CAB39L",
    # More commonly expressed genes
    "GAPDHS", "LDHB", "LDHC", "LDHD",
    "MDH1", "MDH2", "ME1", "ME2", "ME3", "IDH1", "IDH2", "IDH3A",
    "IDH3B", "IDH3G", "OGDH", "OGDHL", "SUCLA2", "SUCLG1", "SUCLG2",
    "SDHA", "SDHB", "SDHC", "SDHD", "FH",
    "ATP5F1A", "ATP5F1B", "ATP5F1C", "ATP5F1D", "ATP5F1E",
    "ATP5MC1", "ATP5MC2", "ATP5MC3", "ATP5MG", "ATP5PB", "ATP5PD",
    "ATP5PF", "ATP5PO", "NDUFB1", "NDUFB2", "NDUFB3", "NDUFB4",
    "NDUFB5", "NDUFB6", "NDUFB7", "NDUFB8", "NDUFB9", "NDUFB10",
    "NDUFB11", "NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4", "NDUFA5",
    "NDUFA6", "NDUFA7", "NDUFA8", "NDUFA9", "NDUFA10", "NDUFA11",
    "NDUFA12", "NDUFA13", "NDUFC1", "NDUFC2", "NDUFS1", "NDUFS2",
    "NDUFS3", "NDUFS4", "NDUFS5", "NDUFS6", "NDUFS7", "NDUFS8",
    "NDUFV1", "NDUFV2", "NDUFV3", "COX4I1", "COX4I2", "COX5A",
    "COX5B", "COX6A1", "COX6A2", "COX6B1", "COX6B2", "COX6C",
    "COX7A1", "COX7A2", "COX7B", "COX7C", "COX8A", "COX8C",
    "UQCRB", "UQCRC1", "UQCRC2", "UQCRFS1", "UQCRH", "UQCRQ",
    "CYC1", "CYCS", "ATP1A1", "ATP1A2", "ATP1A3", "ATP1A4",
    "ATP1B1", "ATP1B2", "ATP1B3", "ATP2A1", "ATP2A2", "ATP2A3",
    "ATP2B1", "ATP2B2", "ATP2B3", "ATP2B4", "ATP6V0A1", "ATP6V0A2",
    "ATP6V0A4", "ATP6V0B", "ATP6V0C", "ATP6V0D1", "ATP6V0D2",
    "ATP6V0E1", "ATP6V0E2", "ATP6V1A", "ATP6V1B1", "ATP6V1B2",
    "ATP6V1C1", "ATP6V1C2", "ATP6V1D", "ATP6V1E1", "ATP6V1E2",
    "ATP6V1F", "ATP6V1G1", "ATP6V1G2", "ATP6V1G3", "ATP6V1H",
    # Proteasome subunits
    "PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMA7",
    "PSMA8", "PSMB1", "PSMB2", "PSMB3", "PSMB4", "PSMB5", "PSMB6",
    "PSMB7", "PSMB10", "PSMB11", "PSMC1", "PSMC2", "PSMC3",
    "PSMC4", "PSMC5", "PSMC6", "PSMD1", "PSMD2", "PSMD3", "PSMD4",
    "PSMD5", "PSMD6", "PSMD7", "PSMD8", "PSMD9", "PSMD10", "PSMD11",
    "PSMD12", "PSMD13", "PSMD14", "PSME1", "PSME2", "PSME3",
    "PSME4", "PSMF1",
    # Translation machinery
    "EIF1", "EIF1AX", "EIF1AY", "EIF1B", "EIF2A", "EIF2B1", "EIF2B2",
    "EIF2B3", "EIF2B4", "EIF2B5", "EIF2S1", "EIF2S2", "EIF2S3",
    "EIF3A", "EIF3B", "EIF3C", "EIF3D", "EIF3E", "EIF3F", "EIF3G",
    "EIF3H", "EIF3I", "EIF3J", "EIF3K", "EIF3L", "EIF3M",
    "EIF4B", "EIF4E", "EIF4E2", "EIF4E3", "EIF4EBP1", "EIF4EBP2",
    "EIF4EBP3", "EIF4G1", "EIF4G2", "EIF4G3", "EIF4H", "EIF5",
    "EIF5A", "EIF5A2", "EIF5B", "EIF6", "ETF1", "GSPT1", "GSPT2",
    # Ribosome biogenesis
    "BOP1", "BMS1", "DDX18", "DDX21", "DDX47", "DDX52", "DDX56",
    "DKC1", "EMG1", "EXOSC1", "EXOSC2", "EXOSC3", "EXOSC4",
    "EXOSC5", "EXOSC6", "EXOSC7", "EXOSC8", "EXOSC9", "EXOSC10",
    "FBL", "GAR1", "GNL2", "GNL3", "GNL3L", "IMP3", "IMP4",
    "LAS1L", "MPHOSPH10", "NAF1", "NHP2", "NOL6", "NOL9", "NOL10",
    "NOL11", "NOP2", "NOP10", "NOP14", "NOP56", "NOP58",
    "NSA2", "PES1", "PPAN", "PWP2", "RCL1", "RPF1", "RPF2",
    "RRP1", "RRP1B", "RRP7A", "RRP8", "RRP9", "RRP12", "RRP15",
    "RRP36", "RSL1D1", "SDAD1", "SNU13",
    "TSR1", "TSR3", "TBL3", "UTP3", "UTP4", "UTP6", "UTP11",
    "UTP14A", "UTP14C", "UTP15", "UTP18", "UTP20", "UTP23", "UTP25",
    "WDR3", "WDR12", "WDR36", "WDR43", "WDR46", "WDR74", "WDR75",
    "XRN2",
    # Complement system
    "C1QA", "C1QB", "C1QC", "C1R", "C1S", "C2", "C3", "C4A", "C4B",
    "C5", "C6", "C7", "C8A", "C8B", "C8G", "C9", "CFB", "CFD",
    "CFH", "CFHR1", "CFHR2", "CFHR3", "CFHR4", "CFHR5", "CFI",
    "CFP", "MASP1", "MASP2", "MBL2", "CD46", "CD55", "CD59",
    "SERPING1", "CLU", "ITGB2", "CR1", "CR2", "C3AR1", "C5AR1",
    "C5AR2",
    # Heat shock / chaperones
    "HSPA2", "HSPA4", "HSPA4L", "HSPA5", "HSPA6", "HSPA8",
    "HSPA9", "HSPA12A", "HSPA12B", "HSPA13", "HSPA14",
    "HSPD1", "HSPE1", "HSP90AB1", "HSP90B1", "TRAP1",
    "DNAJA1", "DNAJA2", "DNAJA3", "DNAJA4", "DNAJB1", "DNAJB2",
    "DNAJB4", "DNAJB5", "DNAJB6", "DNAJB7", "DNAJB8", "DNAJB9",
    "DNAJB11", "DNAJB12", "DNAJB13", "DNAJB14", "DNAJC1", "DNAJC2",
    "DNAJC3", "DNAJC5", "DNAJC6", "DNAJC7", "DNAJC8", "DNAJC9",
    "DNAJC10", "DNAJC11", "DNAJC12", "DNAJC13", "DNAJC14",
    "DNAJC15", "DNAJC16", "DNAJC17", "DNAJC18", "DNAJC19",
    "DNAJC21", "DNAJC22", "DNAJC24", "DNAJC25", "DNAJC27",
    "DNAJC28", "DNAJC30", "CCT2", "CCT3", "CCT4", "CCT5", "CCT6A",
    "CCT6B", "CCT7", "CCT8", "TCP1", "PFDN1", "PFDN2", "PFDN4",
    "PFDN5", "PFDN6", "VBP1",
    # Immunoglobulin superfamily / immune checkpoints
    "PDCD1", "CD274", "PDCD1LG2", "HAVCR2", "LAG3", "TIGIT",
    "BTLA", "ICOS", "ICOSLG", "CD80", "CD86", "CD28", "B7H3",
    "B7H4", "VISTA", "TNFRSF4", "OX40L", "CD27", "CD70", "GITR",
    "CD40", "CD40LG", "CD48", "SLAMF1", "SLAMF6", "SLAMF7",
    "CD2", "CD58", "TNFRSF9", "TNFSF9", "CD226", "DNAM1",
    # Integrins / adhesion
    "ITGA1", "ITGA2", "ITGA2B", "ITGA3", "ITGA4", "ITGA6",
    "ITGA7", "ITGA8", "ITGA9", "ITGA10", "ITGA11", "ITGAD",
    "ITGAE", "ITGAL", "ITGAV", "ITGAX", "ITGB2", "ITGB3",
    "ITGB4", "ITGB5", "ITGB6", "ITGB7", "ITGB8",
    "CDH1", "CDH3", "CDH4", "CDH5", "CDH6", "CDH7", "CDH8",
    "CDH9", "CDH10", "CDH11", "CDH12", "CDH13", "CDH15", "CDH16",
    "CDH17", "CDH18", "CDH19", "CDH20", "CDH22", "CDH23", "CDH24",
    "CDH26", "PCDH1", "PCDH7", "PCDH8", "PCDH9", "PCDH10",
    "PCDH11X", "PCDH11Y", "PCDH12", "PCDH17", "PCDH18", "PCDH19",
    "PCDH20", "PCDHB1", "PCDHB2", "PCDHB3", "PCDHB4", "PCDHB5",
    "PCDHGA1", "PCDHGA2", "PCDHGA3", "PCDHGA4",
    # Gap junction / tight junction
    "TJP1", "TJP2", "TJP3", "OCLN", "CLDN1", "CLDN2", "CLDN3",
    "CLDN4", "CLDN5", "CLDN6", "CLDN7", "CLDN8", "CLDN9", "CLDN10",
    "CLDN11", "CLDN12", "CLDN14", "CLDN15", "CLDN16", "CLDN17",
    "CLDN18", "CLDN19", "CLDN20", "CLDN22", "CLDN23", "CLDN24",
    "CLDN25", "F11R", "JAM2", "JAM3", "MARVELD1", "MARVELD2",
    "MARVELD3",
    # Cytokine signaling (additional SOCS, PIAS)
    "SOCS1", "SOCS2", "SOCS3", "SOCS4", "SOCS5", "SOCS6", "SOCS7",
    "CISH", "PIAS1", "PIAS2", "PIAS3", "PIAS4",
    # Interferon regulatory (additional)
    "IRF2", "IRF3", "IRF4", "IRF5", "IRF6",
    "IFNA1", "IFNA2", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8",
    "IFNA10", "IFNA13", "IFNA14", "IFNA16", "IFNA17", "IFNA21",
    "IFNB1", "IFNG", "IFNL1", "IFNL2", "IFNL3", "IFNL4",
    "IFNK", "IFNW1",
    # Toll-like receptors / innate immunity
    "TLR1", "TLR2", "TLR3", "TLR4", "TLR5", "TLR6", "TLR7",
    "TLR8", "TLR9", "TLR10", "MYD88", "TIRAP", "TICAM1", "TICAM2",
    "TRAM1", "NFKB1", "RELA", "NFKB2", "TRAF6",
    "NOD1", "NOD2", "NLRP1", "NLRP2", "NLRP4", "NLRP5",
    "NLRP6", "NLRP7", "NLRP8", "NLRP9", "NLRP10", "NLRP11",
    "NLRP12", "NLRP13", "NLRP14", "NLRC3", "NLRC4", "NLRC5",
    "NLRX1", "NAIP", "CIITA", "RIG1",
    "STING1", "CGAS", "MAVS", "MDA5",
    # Keratins (non-hair)
    "KRT1", "KRT2", "KRT3", "KRT4", "KRT5", "KRT6A", "KRT6B",
    "KRT6C", "KRT7", "KRT8", "KRT9", "KRT10", "KRT12", "KRT13",
    "KRT14", "KRT15", "KRT16", "KRT17", "KRT18", "KRT19", "KRT20",
    "KRT23", "KRT24", "KRT25", "KRT26", "KRT27", "KRT28",
    "KRT71", "KRT72", "KRT73", "KRT74", "KRT75", "KRT76",
    "KRT77", "KRT78", "KRT79", "KRT80",
    # Selenoproteins
    "SELENOF", "SELENOH", "SELENOI", "SELENOK", "SELENOM",
    "SELENON", "SELENOO", "SELENOP", "SELENOS", "SELENOT",
    "SELENOV", "SELENOW", "DIO1", "DIO2", "DIO3", "TXNRD1",
    "TXNRD2", "TXNRD3", "MSRB1", "SEPHS2",
    # tRNA synthetases
    "AARS1", "CARS1", "DARS1", "EARS2", "FARS2", "GARS1",
    "HARS1", "IARS1", "KARS1", "LARS1", "MARS1", "NARS1",
    "PARS2", "QARS1", "RARS1", "SARS1", "TARS1", "VARS1",
    "WARS2", "YARS1", "EPRS1",
    # Additional commonly detected genes
    "VEGFB", "VEGFC", "VEGFD", "NRP1", "NRP2", "FLT1", "KDR",
    "TEK", "TIE1", "ANGPT1", "ANGPT2", "PECAM1", "ENG", "ACVRL1",
    "CDH5", "VWF", "THBD", "PROCR", "ESAM", "MMRN1", "MMRN2",
    "LYVE1", "PROX1", "FLT4",
    # Metabolism — glycolysis/gluconeogenesis
    "HK1", "HK3", "GCK", "PFKM", "PFKL", "PFKP", "PFKFB1",
    "PFKFB2", "PFKFB4", "BPGM", "PGAM2", "PGAM4",
    "PKM", "PKLR",
    # More miscellaneous well-known genes
    "PCSK1", "PCSK2", "PCSK4", "PCSK5", "PCSK6", "PCSK7", "PCSK9",
    "SPP1", "THBS1", "THBS2", "THBS3", "THBS4", "VCAN", "AGT",
    "REN", "ACE", "ACE2", "AGT", "AGTR1", "AGTR2", "MAS1",
    "SERPINC1", "SERPIND1", "SERPINE2", "SERPINF1", "SERPINF2",
    "SERPING1", "SERPINH1", "SERPINI1", "SERPINI2", "SERPINB1",
    "SERPINB2", "SERPINB3", "SERPINB4", "SERPINB5", "SERPINB6",
    "SERPINB7", "SERPINB8", "SERPINB9", "SERPINB10", "SERPINB11",
    "SERPINB12", "SERPINB13", "SERPINA3", "SERPINA4", "SERPINA5",
    "SERPINA6", "SERPINA7", "SERPINA9", "SERPINA10", "SERPINA11",
    "SERPINA12",
    # Additional housekeeping / common
    "GAPDH", "ACTB", "B2M", "HMBS", "TBP", "HPRT1", "YWHAZ",
    "RPLP0", "RPLP1", "RPLP2", "RPL10", "RPL10A", "RPL12",
    "RPL13", "RPL14", "RPL15", "RPL17", "RPL18", "RPL18A",
    "RPL21", "RPL22", "RPL23", "RPL23A", "RPL24", "RPL26",
    "RPL26L1", "RPL27A", "RPL28", "RPL29", "RPL30", "RPL31",
    "RPL32", "RPL34", "RPL35", "RPL36", "RPL36A", "RPL37",
    "RPL38", "RPL39", "RPL40", "RPL41",
    "RPS2", "RPS3A", "RPS4Y1", "RPS5", "RPS6", "RPS7", "RPS9",
    "RPS10", "RPS11", "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A",
    "RPS16", "RPS17", "RPS19", "RPS20", "RPS21", "RPS23",
    "RPS24", "RPS25", "RPS26", "RPS27", "RPS28", "RPS29",
    "UBA1", "UBA2", "UBA3", "UBA5", "UBA6", "UBA7",
]
# fmt: on


# ---------------------------------------------------------------------------
# Biology model: which programs show signal for each cohort biological_program
# ---------------------------------------------------------------------------

# Maps signature_id -> its primary biological program for grouping
SIGNATURE_PROGRAM_MAP = {
    # durable
    "hallmark_ifng_response": "interferon",
    "hallmark_ifna_response": "interferon",
    "hallmark_inflammatory_response": "inflammation",
    "hallmark_hypoxia": "hypoxia",
    "hallmark_e2f_targets": "proliferation",
    "hallmark_emt": "emt",
    "hallmark_tnfa_nfkb": "inflammation",
    "curated_senescence": "senescence",
    # brittle
    "brittle_small_study_inflammation": "brittle",
    "brittle_platform_specific": "brittle",
    "brittle_single_tissue": "brittle",
    "brittle_overfit_noise": "brittle",
    # mixed
    "mixed_emt_inflammation": "mixed_emt",
    "mixed_hypoxia_stress": "mixed_hypoxia",
    "mixed_senescence_proliferation": "mixed_senescence",
    # confounded
    "confounded_proliferation": "confounded_prolif",
    "confounded_ribosomal": "confounded_ribo",
    "confounded_immune_infiltrate": "confounded_immune",
    "stealth_confounded_hypoxia_prolif": "stealth_hypoxia",
    "stealth_confounded_inflam_immune": "stealth_inflam",
    # lowcov
    "lowcov_olfactory_1": "lowcov",
    "lowcov_olfactory_2": "lowcov",
    # blind durable
    "blind_durable_ifn_composite": "interferon",
    "blind_durable_hypoxia_core": "hypoxia",
    # blind brittle
    "blind_brittle_random": "brittle",
    # blind mixed
    "blind_mixed_emt_stress": "mixed_emt",
    # blind confounded
    "blind_confounded_prolif_inflam": "confounded_prolif",
    "blind_confounded_stealth_emt_immune": "stealth_emt",
    # blind lowcov
    "blind_lowcov_keratin": "lowcov",
}

# Effect size map: cohort_program -> { signature_program -> Cohen's d }
# Positive d = upregulated in cases; negative = downregulated
EFFECT_SIZE_MAP = {
    "inflammation": {
        "inflammation": 0.85,        # strong
        "interferon": 0.45,          # moderate (co-activation)
        "hypoxia": 0.05,             # noise
        "proliferation": 0.10,       # noise
        "emt": 0.05,                 # noise
        "senescence": 0.15,          # weak
    },
    "interferon": {
        "interferon": 0.90,          # strong
        "inflammation": 0.40,        # moderate (co-activation)
        "hypoxia": 0.05,             # noise
        "proliferation": 0.05,       # noise
        "emt": 0.05,                 # noise
        "senescence": 0.05,          # noise
    },
    "hypoxia": {
        "hypoxia": 0.80,             # strong
        "emt": 0.40,                 # moderate (hypoxia drives EMT)
        "inflammation": 0.10,        # weak
        "interferon": 0.05,          # noise
        "proliferation": 0.15,       # weak
        "senescence": 0.10,          # weak
    },
    "proliferation": {
        "proliferation": 0.85,       # strong
        "senescence": -0.35,         # moderate inverse
        "inflammation": 0.10,        # noise
        "interferon": 0.05,          # noise
        "hypoxia": 0.10,             # weak
        "emt": 0.05,                 # noise
    },
    "emt": {
        "emt": 0.80,                 # strong
        "hypoxia": 0.35,             # moderate (EMT ↔ hypoxia)
        "inflammation": 0.15,        # weak
        "interferon": 0.05,          # noise
        "proliferation": 0.05,       # noise
        "senescence": 0.10,          # weak
    },
    "senescence": {
        "senescence": 0.75,          # strong
        "proliferation": -0.30,      # inverse (senescence = growth arrest)
        "inflammation": 0.35,        # moderate (SASP is inflammatory)
        "interferon": 0.10,          # weak
        "hypoxia": 0.05,             # noise
        "emt": 0.10,                 # weak
    },
    "mixed": {
        # mixed-program cohorts: many programs show moderate signal
        "inflammation": 0.45,
        "interferon": 0.25,
        "hypoxia": 0.30,
        "proliferation": 0.20,
        "emt": 0.35,
        "senescence": 0.25,
    },
}

# For confounded signatures, effect on confounder genes in each cohort
CONFOUNDER_EFFECT_MAP = {
    "inflammation": {
        "proliferation_cell_cycle": 0.30,
        "generic_inflammation": 0.80,
        "stress_response": 0.40,
        "batch_ribosomal": 0.40,
        "immune_infiltration": 0.70,
    },
    "interferon": {
        "proliferation_cell_cycle": 0.10,
        "generic_inflammation": 0.50,
        "stress_response": 0.30,
        "batch_ribosomal": 0.35,
        "immune_infiltration": 0.55,
    },
    "hypoxia": {
        "proliferation_cell_cycle": 0.50,
        "generic_inflammation": 0.20,
        "stress_response": 0.55,
        "batch_ribosomal": 0.30,
        "immune_infiltration": 0.15,
    },
    "proliferation": {
        "proliferation_cell_cycle": 0.85,
        "generic_inflammation": 0.15,
        "stress_response": 0.25,
        "batch_ribosomal": 0.45,
        "immune_infiltration": 0.10,
    },
    "emt": {
        "proliferation_cell_cycle": 0.15,
        "generic_inflammation": 0.30,
        "stress_response": 0.40,
        "batch_ribosomal": 0.30,
        "immune_infiltration": 0.25,
    },
    "senescence": {
        "proliferation_cell_cycle": -0.25,
        "generic_inflammation": 0.45,
        "stress_response": 0.50,
        "batch_ribosomal": 0.35,
        "immune_infiltration": 0.30,
    },
    "mixed": {
        "proliferation_cell_cycle": 0.35,
        "generic_inflammation": 0.45,
        "stress_response": 0.35,
        "batch_ribosomal": 0.35,
        "immune_infiltration": 0.40,
    },
}


def load_signatures(sig_path: Path) -> dict[str, list[tuple[str, str, float]]]:
    """Load all_signatures.tsv -> {sig_id: [(gene, direction, weight), ...]}."""
    df = pd.read_csv(sig_path, sep="\t")
    sigs: dict[str, list[tuple[str, str, float]]] = {}
    for _, row in df.iterrows():
        sid = row["signature_id"]
        sigs.setdefault(sid, []).append(
            (row["gene_symbol"], row["direction"], float(row["weight"]))
        )
    return sigs


def load_confounders(conf_path: Path) -> dict[str, list[str]]:
    """Load confounder_panel.yaml -> {set_name: [gene, ...]}."""
    with open(conf_path) as f:
        data = yaml.safe_load(f)
    result = {}
    for set_name, info in data["sets"].items():
        result[set_name] = [g.upper() for g in info["genes_up"]]
    return result


def load_cohort_manifest(manifest_path: Path) -> pd.DataFrame:
    """Load cohort_manifest.tsv."""
    return pd.read_csv(manifest_path, sep="\t")


def build_gene_universe(
    signatures: dict[str, list[tuple[str, str, float]]],
    confounders: dict[str, list[str]],
) -> list[str]:
    """Build a gene universe: all signature + confounder + background genes."""
    genes = set()

    # Add all signature genes
    for gene_list in signatures.values():
        for gene, _, _ in gene_list:
            genes.add(gene.upper())

    # Add all confounder genes
    for gene_list in confounders.values():
        for gene in gene_list:
            genes.add(gene.upper())

    # Add background genes (skip any that overlap)
    bg_set = set()
    for g in BACKGROUND_GENES:
        g_upper = g.upper()
        if g_upper not in genes:
            bg_set.add(g_upper)

    # We want ~10000 total genes
    target = 10000
    current = len(genes)
    bg_needed = max(0, target - current)
    bg_sorted = sorted(bg_set)
    genes.update(bg_sorted[:bg_needed])

    # If still short (unlikely with our list), add numbered placeholder genes
    if len(genes) < target:
        idx = 1
        while len(genes) < target:
            placeholder = f"BKGD{idx:04d}"
            if placeholder not in genes:
                genes.add(placeholder)
            idx += 1

    return sorted(genes)


def get_signature_effect_size(
    sig_id: str, cohort_program: str, rng: np.random.Generator
) -> float:
    """Determine the Cohen's d for a signature in a given cohort.

    Returns the effect size (may be 0 for no effect).
    """
    sig_program = SIGNATURE_PROGRAM_MAP.get(sig_id)
    if sig_program is None:
        return 0.0

    # ---- Low-coverage sentinels: no expression effect at all ----
    if sig_program == "lowcov":
        return 0.0

    # ---- Brittle signatures: at most 1 cohort shows random weak effect ----
    if sig_program == "brittle":
        # Use a deterministic hash to pick at most 1 cohort for each brittle sig
        # Hash the sig_id to pick a "lucky" cohort index
        sig_hash = hash(sig_id) % 14  # 14 cohorts
        cohort_programs = [
            "inflammation", "interferon", "hypoxia", "proliferation",
            "emt", "senescence", "mixed",
            "inflammation", "interferon", "hypoxia", "proliferation",
            "emt", "senescence", "mixed",
        ]
        if cohort_programs[sig_hash] == cohort_program:
            return rng.uniform(0.15, 0.30)  # weak signal in just 1 cohort
        return rng.uniform(-0.05, 0.05)  # noise

    # ---- Confounded signatures ----
    if sig_program.startswith("confounded_"):
        # These are driven by confounder genes, not their own gene sets
        # The signature genes themselves show the CONFOUNDER effects
        confounder_name = {
            "confounded_prolif": "proliferation_cell_cycle",
            "confounded_ribo": "batch_ribosomal",
            "confounded_immune": "immune_infiltration",
        }.get(sig_program)
        if confounder_name and cohort_program in CONFOUNDER_EFFECT_MAP:
            return CONFOUNDER_EFFECT_MAP[cohort_program].get(confounder_name, 0.0)
        return 0.0

    # ---- Stealth confounded: real signal + confounder ----
    if sig_program.startswith("stealth_"):
        # Real biological signal component
        bio_program = {
            "stealth_hypoxia": "hypoxia",
            "stealth_inflam": "inflammation",
            "stealth_emt": "emt",
        }.get(sig_program, "")
        if cohort_program in EFFECT_SIZE_MAP:
            return EFFECT_SIZE_MAP[cohort_program].get(bio_program, 0.0)
        return 0.0

    # ---- Mixed signatures: partial signal ----
    if sig_program.startswith("mixed_"):
        bio_program = {
            "mixed_emt": "emt",
            "mixed_hypoxia": "hypoxia",
            "mixed_senescence": "senescence",
        }.get(sig_program, "")
        if cohort_program in EFFECT_SIZE_MAP:
            base_d = EFFECT_SIZE_MAP[cohort_program].get(bio_program, 0.0)
            # Mixed sigs have incoherent genes, so effective d is reduced
            return base_d * 0.6
        return 0.0

    # ---- Durable signatures: straightforward lookup ----
    if cohort_program in EFFECT_SIZE_MAP:
        return EFFECT_SIZE_MAP[cohort_program].get(sig_program, 0.0)

    return 0.0


def generate_cohort(
    cohort_row: pd.Series,
    gene_universe: list[str],
    signatures: dict[str, list[tuple[str, str, float]]],
    confounders: dict[str, list[str]],
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Generate expression matrix and phenotype table for one cohort.

    Returns (expression_df, phenotype_df).
    """
    cohort_id = cohort_row["cohort_id"]
    platform = cohort_row["platform"]
    cohort_program = cohort_row["biological_program"]
    n_samples = int(cohort_row["sample_count"])
    case_label = cohort_row["case_label"]
    control_label = cohort_row["control_label"]

    n_genes = len(gene_universe)

    # Split samples roughly 50/50 case/control
    n_case = n_samples // 2
    n_control = n_samples - n_case

    # Platform-specific baseline
    if platform == "affymetrix":
        baseline_mu = 5.5
    else:
        baseline_mu = 6.5
    baseline_sd = 1.5

    # Per-gene variance (log2-scale)
    gene_variance = rng.uniform(0.5, 2.0, size=n_genes)
    gene_sd = np.sqrt(gene_variance)

    # Per-gene baseline (slight gene-to-gene variation)
    gene_baseline = rng.normal(baseline_mu, baseline_sd, size=n_genes)
    # Clip to realistic range
    gene_baseline = np.clip(gene_baseline, 2.0, 12.0)

    # Initialize expression matrix: genes x samples
    # Control samples first, then case samples
    expr = np.zeros((n_genes, n_samples))
    gene_to_idx = {g: i for i, g in enumerate(gene_universe)}

    for j in range(n_samples):
        expr[:, j] = rng.normal(gene_baseline, gene_sd)

    # Build set of genes that belong to each signature, for effect application
    sig_genes_sets: dict[str, set[str]] = {}
    for sig_id, gene_list in signatures.items():
        sig_genes_sets[sig_id] = {g.upper() for g, _, _ in gene_list}

    # Collect ALL genes that appear in ANY brittle signature.
    # These genes must NOT receive differential effects from ANY source,
    # ensuring that brittle signatures measure as noise across all cohorts.
    # Durable signatures lose a few shared genes but have enough unique
    # genes to maintain their signal.
    brittle_sig_ids = {
        sid for sid, prog in SIGNATURE_PROGRAM_MAP.items()
        if prog == "brittle"
    }
    all_brittle_genes: set[str] = set()
    for sid in brittle_sig_ids:
        if sid in sig_genes_sets:
            all_brittle_genes.update(sig_genes_sets[sid])

    # Apply signature effects to CASE samples only
    for sig_id, gene_list in signatures.items():
        d = get_signature_effect_size(sig_id, cohort_program, rng)
        if abs(d) < 0.01:
            continue

        # Skip brittle signatures entirely (their effect is ~0 by design)
        sig_prog = SIGNATURE_PROGRAM_MAP.get(sig_id, "")
        if sig_prog == "brittle":
            continue

        for gene, direction, weight in gene_list:
            gene_upper = gene.upper()
            idx = gene_to_idx.get(gene_upper)
            if idx is None:
                continue

            # Skip genes that appear in ANY brittle signature to prevent
            # brittle sigs from piggybacking on durable/confounded effects
            if gene_upper in all_brittle_genes:
                continue

            # Effect magnitude = d * weight * gene_sd[idx]
            effect = d * weight * gene_sd[idx]
            if direction == "down":
                effect = -effect

            # Apply to case samples (indices n_control:)
            expr[idx, n_control:] += effect

    # Apply confounder effects to case samples
    if cohort_program in CONFOUNDER_EFFECT_MAP:
        for conf_set, conf_genes in confounders.items():
            d_conf = CONFOUNDER_EFFECT_MAP[cohort_program].get(conf_set, 0.0)
            if abs(d_conf) < 0.01:
                continue
            for gene in conf_genes:
                gene_upper = gene.upper()
                idx = gene_to_idx.get(gene_upper)
                if idx is None:
                    continue
                # Skip brittle genes
                if gene_upper in all_brittle_genes:
                    continue
                # Only apply if this gene is NOT already in a durable signature
                # (to avoid double-counting; durable sigs already got their effect)
                is_in_durable = False
                for sid, gene_set in sig_genes_sets.items():
                    sig_prog_inner = SIGNATURE_PROGRAM_MAP.get(sid, "")
                    if sig_prog_inner in (
                        "inflammation", "interferon", "hypoxia",
                        "proliferation", "emt", "senescence",
                    ) and gene_upper in gene_set:
                        is_in_durable = True
                        break
                if not is_in_durable:
                    effect = d_conf * gene_sd[idx]
                    expr[idx, n_control:] += effect

    # Clip to realistic log2 range
    expr = np.clip(expr, 1.5, 15.0)

    # Build sample IDs
    sample_ids = []
    phenotypes = []
    for j in range(n_control):
        sample_ids.append(f"{cohort_id}_ctrl_{j + 1:03d}")
        phenotypes.append(control_label)
    for j in range(n_case):
        sample_ids.append(f"{cohort_id}_case_{j + 1:03d}")
        phenotypes.append(case_label)

    # Create DataFrames
    expr_df = pd.DataFrame(expr, index=gene_universe, columns=sample_ids)
    expr_df.index.name = "gene_symbol"

    pheno_df = pd.DataFrame(
        {"sample": sample_ids, "phenotype": phenotypes}
    )

    return expr_df, pheno_df


def compute_effect_sizes_audit(
    expr_df: pd.DataFrame,
    pheno_df: pd.DataFrame,
    signatures: dict[str, list[tuple[str, str, float]]],
    case_label: str,
    control_label: str,
) -> dict[str, float]:
    """Compute observed Cohen's d for each signature in a cohort."""
    case_mask = pheno_df["phenotype"] == case_label
    ctrl_mask = pheno_df["phenotype"] == control_label
    case_samples = pheno_df.loc[case_mask, "sample"].tolist()
    ctrl_samples = pheno_df.loc[ctrl_mask, "sample"].tolist()

    result = {}
    for sig_id, gene_list in signatures.items():
        genes = [g.upper() for g, _, _ in gene_list]
        genes_present = [g for g in genes if g in expr_df.index]
        if len(genes_present) < 3:
            result[sig_id] = float("nan")
            continue

        case_vals = expr_df.loc[genes_present, case_samples].values.mean(axis=0)
        ctrl_vals = expr_df.loc[genes_present, ctrl_samples].values.mean(axis=0)

        pooled_sd = np.sqrt(
            (np.var(case_vals, ddof=1) + np.var(ctrl_vals, ddof=1)) / 2
        )
        if pooled_sd < 1e-8:
            result[sig_id] = 0.0
        else:
            result[sig_id] = (case_vals.mean() - ctrl_vals.mean()) / pooled_sd

    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate synthetic cohort data")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    print(f"Loading signatures from {CURATION_DIR / 'all_signatures.tsv'}")
    signatures = load_signatures(CURATION_DIR / "all_signatures.tsv")
    print(f"  Loaded {len(signatures)} signatures")

    print(f"Loading confounders from {CONFIG_DIR / 'confounder_panel.yaml'}")
    confounders = load_confounders(CONFIG_DIR / "confounder_panel.yaml")
    print(f"  Loaded {len(confounders)} confounder sets")

    print(f"Loading cohort manifest from {CONFIG_DIR / 'cohort_manifest.tsv'}")
    manifest = load_cohort_manifest(CONFIG_DIR / "cohort_manifest.tsv")
    print(f"  Loaded {len(manifest)} cohorts")

    print("Building gene universe...")
    gene_universe = build_gene_universe(signatures, confounders)
    print(f"  Gene universe: {len(gene_universe)} genes")

    # Verify all signature and confounder genes are in universe
    missing_sig = set()
    for gene_list in signatures.values():
        for gene, _, _ in gene_list:
            if gene.upper() not in set(gene_universe):
                missing_sig.add(gene.upper())
    missing_conf = set()
    for gene_list in confounders.values():
        for gene in gene_list:
            if gene.upper() not in set(gene_universe):
                missing_conf.add(gene.upper())
    if missing_sig:
        print(f"  WARNING: {len(missing_sig)} signature genes missing from universe: {missing_sig}")
    if missing_conf:
        print(f"  WARNING: {len(missing_conf)} confounder genes missing from universe: {missing_conf}")

    # Create output directories
    MATRIX_DIR.mkdir(parents=True, exist_ok=True)
    PHENO_DIR.mkdir(parents=True, exist_ok=True)

    total_samples = 0
    all_effect_sizes = {}

    for _, cohort_row in manifest.iterrows():
        cohort_id = cohort_row["cohort_id"]
        print(f"\nGenerating cohort: {cohort_id}")
        print(f"  Platform: {cohort_row['platform']}, Program: {cohort_row['biological_program']}, "
              f"Samples: {cohort_row['sample_count']}")

        expr_df, pheno_df = generate_cohort(
            cohort_row, gene_universe, signatures, confounders, rng
        )

        # Save
        matrix_path = MATRIX_DIR / f"{cohort_id}.tsv"
        pheno_path = PHENO_DIR / f"{cohort_id}.tsv"
        expr_df.to_csv(matrix_path, sep="\t", float_format="%.4f")
        pheno_df.to_csv(pheno_path, sep="\t", index=False)

        n = int(cohort_row["sample_count"])
        total_samples += n
        print(f"  Saved: {matrix_path.name} ({len(gene_universe)} genes x {n} samples)")

        # Audit: compute effect sizes for key signatures
        effects = compute_effect_sizes_audit(
            expr_df, pheno_df, signatures,
            cohort_row["case_label"], cohort_row["control_label"]
        )
        all_effect_sizes[cohort_id] = effects

        # Print a few key effect sizes
        key_sigs = [
            "hallmark_inflammatory_response", "hallmark_ifng_response",
            "hallmark_hypoxia", "hallmark_e2f_targets", "hallmark_emt",
            "curated_senescence", "brittle_overfit_noise",
            "confounded_proliferation",
        ]
        for s in key_sigs:
            d = effects.get(s, float("nan"))
            if not np.isnan(d):
                print(f"    {s}: d={d:.3f}")

    # Write freeze_audit.json
    audit = {
        "cohort_count": len(manifest),
        "total_samples": total_samples,
        "gene_universe_size": len(gene_universe),
        "generation_seed": args.seed,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "cohort_ids": manifest["cohort_id"].tolist(),
        "method": "synthetic_cohorts_modeling_geo_characteristics",
        "description": (
            "Biologically realistic synthetic cohort expression matrices. "
            "Durable signatures show cross-cohort consistency (Cohen's d ~0.5-1.0). "
            "Brittle signatures show effect in at most 1 cohort. "
            "Confounded signatures are driven by confounder gene overlap. "
            "Stealth confounded signatures have real signal + confounder contamination."
        ),
    }
    audit_path = FREEZE_DIR / "freeze_audit.json"
    with open(audit_path, "w") as f:
        json.dump(audit, f, indent=2)
    print(f"\nFreeze audit written to {audit_path}")

    print(f"\n{'=' * 60}")
    print(f"DONE: {len(manifest)} cohorts, {total_samples} total samples, "
          f"{len(gene_universe)} genes")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
