# Gene2Phenotype (G2P) scripts

The Gene2Phenotype scripts repo provides a set of scripts to add, delete or update G2P data.

It is part of the Gene2Phenotype project, established by David FitzPatrick in 2012 with the aim of accelerating the diagnosis of children with developmental disorders. In 2014, the database moved from the University of Edinburgh to EMBL-EBI and a dedicated website was launched to improve data accessibility. G2P has been generalised to cover other disease areas and extended to capture additional information to provide a more detailed understanding of disease mechanism. To read more about the project please visit the [G2P website](https://www.ebi.ac.uk/gene2phenotype/about/project).

### Requirements

- Python 3.10

## `load_data/load_records.py`

Validates and imports G2P records from a CSV or XLSX file. The script:

- reads the input file
- preloads reference data from the G2P MySQL database
- validates publications through the internal G2P API
- validates phenotypes through the JAX HPO API
- creates draft curations and publishes them through the internal G2P API

### Input file format

Supported formats:

- `.csv`
- `.xlsx`

Mandatory columns:

- `gene symbol`
- `disease name`
- `allelic requirement`
- `molecular mechanism`
- `confidence`
- `publication`
- `panel`
- `inferred variant consequence`
- `evidence based variant consequence`

At least one of these two variant consequence columns must be populated for each row.

Optional columns:

- `publication molecular mechanism evidence`
- `molecular mechanism categorisation`
- `publication comment`
- `publication families`
- `publication affected individuals`
- `publication consanguinity`
- `publication families ancestry`
- `disease mim`
- `disease MONDO`
- `cross cutting modifier`
- `publication phenotypes`
- `publication variant descriptions`
- `publication variant types`
- `public comment`
- `private comment`

### Column notes

- `publication` must contain a PMID. Numeric spreadsheet values such as `12345` or `12345.0` are accepted.
- `disease name` is expected to contain `<gene symbol>-related`.
- `inferred variant consequence` and `evidence based variant consequence` accept one or more G2P ontology terms separated by `;`. Input values may use `_`; the script converts `_` to spaces before validation.
- `publication molecular mechanism evidence`, `cross cutting modifier`, `publication phenotypes`, `publication variant descriptions`, and `publication variant types` accept multiple values separated by `;`.
- `publication variant types` may include inheritance qualifiers in parentheses, for example `missense_variant (de_novo)` or `inframe_insertion (unknown_inheritance)`.
- `publication phenotypes` expects HPO accessions such as `HP:0001897`.
- `disease mim` is validated as an OMIM identifier and `disease MONDO` as a MONDO identifier when those columns are present.

## `create_drafts/create_draft_records.py`

Creates automatic curation drafts in G2P from source JSON files such as ClinGen or PanelApp exports. The script:

- reads a source JSON file
- fetches publication details through the internal G2P API
- fetches disease cross-reference details through the internal G2P API
- checks existing automatic drafts to avoid duplicate session names
- inserts new automatic drafts through the internal G2P API

### CLI inputs

Required arguments:

- `--source`
- `--input_file`
- `--config`
- `--panel`

Optional arguments:

- `--output_file`

Config requirements:

- `[api] api_url`
- `[api] api_username`
- `[api] api_password`

### Input file format

Supported format:

- `.json`

The script expects a list of source records. Common input fields used by the script are:

- `gene_symbol`
- `pmids`
- `allelic_requirement`
- `disease`
- `disease_id`
- `mondo_id`
- `mechanism`
- `evidence`
- `phenotypes`
- `confidence`
- `url`
- `comments`
- `status`
- `type`

### Input behavior

- A record is considered for draft creation when `status` is `approved` or `needs changes`, or when `type` is present and not equal to `in_g2p`.
- Records are skipped as already represented when `comments.g2p_data` contains text such as `matches the existing g2p` or `draft matches`.
- `pmids` is expected to be a list of PMIDs. Each PMID is resolved through the internal G2P publication endpoint.
- `allelic_requirement` is mapped from source labels such as `autosomal recessive`, `autosomal dominant`, and `X-linked` to G2P values.
- `mechanism` is checked against a limited supported list: `loss of function`, `gain of function`, `dominant negative`, `undetermined non loss of function`, and `undetermined`.
- `disease_id` may contain comma-separated OMIM identifiers. `mondo_id` is resolved as a MONDO identifier.

### Draft content

The script creates automatic drafts with:

- `status` set to `automatic`
- `panels` set from `--panel`
- `source_data.name` set from `--source`
- publication details populated from the internal G2P API
- disease cross-reference details populated from the internal G2P API

Important limitations:

- The script creates drafts, not published records.
- It does not build a fully publishable record. Fields such as `confidence`, `disease.disease_name`, `molecular_mechanism.name`, `molecular_mechanism.support`, and `variant_consequences` may remain empty and require curator completion later.
- Duplicate checking before insertion is based on automatic draft `session_name`, not a full JSON equality comparison.
