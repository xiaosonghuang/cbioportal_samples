cBioPortal Sample Query Script
==============================

cbioportal_samples.py queries the public cBioPortal REST API
(https://www.cbioportal.org/api/swagger-ui/index.html) directly using Python's
requests library (without bravado or cbio_py) to list samples for a given
cancer type.

Output columns:
  - Cancer Type
  - Patient ID
  - Sample ID
  - Study ID


Requirements
------------
- Python 3.9+
- requests

Install dependencies:
  pip install -r requirements.txt


Usage
-----
Write tab-delimited (TSV) output to stdout (log messages go to stderr):
  python3 cbioportal_samples.py acc > acc_samples.tsv

Write TSV or CSV output to a file:
  python3 cbioportal_samples.py "Breast Cancer" -o breast_samples.tsv
  python3 cbioportal_samples.py "Adrenocortical Carcinoma" --format csv -o acc_samples.csv

Match only a single OncoTree cancer type without child subtypes:
  python3 cbioportal_samples.py breast --exact-type -o breast_exact.tsv

Show all command-line options:
  python3 cbioportal_samples.py --help


Example Output (python3 cbioportal_samples.py acc)
--------------------------------------------------
Cancer Type	Patient ID	Sample ID	Study ID
Adrenocortical Carcinoma	TCGA-OR-A5J1	TCGA-OR-A5J1-01	acc_tcga
Adrenocortical Carcinoma	TCGA-OR-A5J2	TCGA-OR-A5J2-01	acc_tcga
Adrenocortical Carcinoma	TCGA-OR-A5J3	TCGA-OR-A5J3-01	acc_tcga


Notes on cBioPortal Data
------------------------
- Cancer types in /api/cancer-types are organized in the OncoTree hierarchy via
  parent links. By default, querying an organ or mid-level category (such as
  breast, nsclc, coadread, or Melanoma) also includes child subtypes so studies
  annotated under more specific histological subtypes (such as brca, luad, or
  skcm) are included. Pass --exact-type to disable subtype expansion.
- Organ entries in OncoTree use names like "Breast", "Lung", and "Prostate"
  without the word "Cancer", so inputs like "Breast Cancer" or "Prostate Cancer"
  are matched to the corresponding OncoTree entry. If an input is not found,
  close matches from /api/cancer-types are listed in the error message.
- Studies are filtered by exact study["cancerTypeId"] membership from
  /api/studies rather than /api/studies?keyword=..., since keyword performs
  substring matching across study names and IDs.
- Pan-cancer studies assigned cancerTypeId "mixed" do not belong to a single
  cancer type at the study level and are excluded unless "mixed" is queried
  directly.
- Samples are fetched per study from /api/studies/{studyId}/samples using
  pageNumber and pageSize pagination, and each study's retrieved sample count is
  checked against study["allSampleCount"]. Rows are kept per (Study ID,
  Sample ID) because the same TCGA patient/sample can appear across multiple
  study releases (e.g., acc_tcga, acc_tcga_gdc, acc_tcga_pan_can_atlas_2018).
