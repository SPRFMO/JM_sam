# JM_sam Context

This repository is the canonical home for State-Space Assessment Model analyses supporting the SPRFMO jack mackerel assessment and MSE work.

## Language

**JM_sam**:
The repository that owns the SAM workflow, Quarto source reports, rendered diagnostics, and supporting SAM data products.
_Avoid_: temporary SAM folder, copied SAM folder

**JM_SCW_prep/sam**:
The former SAM working folder in the SCW preparation repository.
_Avoid_: canonical SAM source

**Provenance Artifact**:
A frozen rendered output retained only to document what existed in the former SAM working folder at migration time.
_Avoid_: active report, source of truth

**SAM Diagnostics Report**:
The Quarto-authored report that presents the active SAM diagnostic analysis.
_Avoid_: index page

## Relationships

- **JM_sam** owns the active **SAM Diagnostics Report**.
- **JM_SCW_prep/sam** may be referenced by **Provenance Artifacts**, but it is not the source of truth after migration.
- A **Provenance Artifact** can be removed after the migration history is no longer needed for comparison.

## Example dialogue

> **Dev:** "Should I update the SAM report in `JM_SCW_prep/sam` too?"
> **Domain expert:** "No. Update the **SAM Diagnostics Report** in **JM_sam**; `JM_SCW_prep/sam` is only the former working location."

## Flagged ambiguities

- "reconcile repositories" was resolved to mean fully migrate SAM ownership to **JM_sam**, with any retained `JM_SCW_prep/sam` output treated only as a **Provenance Artifact**.
