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
