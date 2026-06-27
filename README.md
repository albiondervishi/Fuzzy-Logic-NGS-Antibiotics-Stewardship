# Fuzzy Antibiotic Stewardship Engine (FASE)

A clinician-readable, rule-based decision support framework for antibiotic stewardship in critically ill patients integrating biomarker kinetics, microbiology, and patient severity into actionable therapeutic recommendations.

## Overview

The **Fuzzy Antibiotic Stewardship Engine (FASE)** is an explainable clinical decision support system designed to assist antimicrobial stewardship in the ICU. The framework combines longitudinal changes in inflammatory biomarkers, severity scores, and microbiological findings using a fuzzy logic approach to generate four stewardship recommendations:

* 🔴 **Escalate therapy**
* 🟢 **De-escalate therapy**
* 🟠 **Stop therapy**
* 🔵 **No change / Continue therapy**

Unlike black-box machine learning approaches, FASE is fully interpretable and produces clinician-readable IF–THEN rules that can be directly audited and adapted to local stewardship protocols.

---

## Clinical Motivation

Antibiotic stewardship decisions in critically ill patients are often challenging because clinicians must simultaneously interpret:

* Dynamic biomarker trajectories (PCT, CRP, leukocytes, temperature)
* Severity of illness (SAPS II)
* Blood culture and metagenomic next-generation sequencing (mNGS) results
* Clinical response over time

Current stewardship recommendations are frequently qualitative and difficult to operationalize in a reproducible manner. FASE addresses this gap by translating complex clinical information into transparent decision rules.

---

## Input Variables

### Biomarkers

* Procalcitonin (PCT)
* C-reactive protein (CRP)
* Leukocyte count
* Body temperature

### Clinical Severity

* SAPS II score
* Infection complexity
* Age

### Microbiology

* Blood culture results
* Metagenomic next-generation sequencing (mNGS)

---

## Core Concepts

The algorithm computes three interpretable fuzzy domains:

### Danger Score

Reflects evidence of ongoing infection, deterioration, or alarm signals.

Examples:

* Rising PCT or CRP
* Worsening SAPS II
* Positive microbiology
* Persistent fever

---

### Response Score

Reflects evidence of clinical improvement.

Examples:

* Falling PCT or CRP
* Defervescence
* Normalizing leukocyte count
* Improving organ dysfunction

---

### Normality Score

Represents the degree of physiological normalization.

Examples:

* Biomarkers returning to normal ranges
* Afebrile state
* Stable clinical condition

---

## Decision Logic

Priority order:

```text
Escalate
    ↓
Stop
    ↓
De-escalate
    ↓
No change
```

### Escalate

```text
IF hard alarm present
OR danger score ≥ threshold
THEN escalate therapy
```

---

### Stop

```text
IF adequate response
AND low danger score
AND negative microbiology
THEN stop therapy
```

---

### De-escalate

```text
IF response score is high
AND danger score remains low
THEN de-escalate therapy
```

---

### No Change

```text
IF signals are mixed or uncertain
THEN continue therapy and reassess
```

---

## Validation

The project includes:

* Multinomial classification models
* Leave-One-Out Cross-Validation (LOOCV)
* One-vs-rest ROC analysis
* Confusion matrix evaluation
* Clinician-readable IF–THEN rules

---

## Key Features

✅ Explainable AI (XAI)

✅ Human-readable decision rules

✅ Integration of mNGS and blood cultures

✅ Longitudinal biomarker interpretation

✅ Suitable for antimicrobial stewardship programs

✅ Easily adaptable to local clinical protocols

---

## Repository Structure

```text
├── data/
├── scripts/
├── figures/
├── results/
├── manuscript/
└── README.md
```

---

## Potential Applications

* Intensive Care Units (ICU)
* Antimicrobial Stewardship Programs (AMS)
* Sepsis management
* Decision support systems
* Prospective clinical validation studies

---

---

## Disclaimer

This software is intended for research purposes only and should not replace clinical judgment. Clinical implementation requires prospective validation and local governance approval.
