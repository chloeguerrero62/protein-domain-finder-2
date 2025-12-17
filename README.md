# Protein Domain Finder 2

A toolkit for identifying protein domains using machine learning and sequence analysis.

## Features
- Jupyter Notebooks and Python scripts for end-to-end domain prediction workflow
- Customizable with your own sequence data
- Integrated evaluation and optimization documentation

## Repository Structure

```
.
├── src/                # Core Python modules and functions
├── notebooks/          # Jupyter notebooks for analysis, prototyping, results
├── scripts/            # Shell and utility scripts for workflow automation
├── data/               # (Not tracked) Raw or processed data; add README.md for download instructions
├── docs/               # Project documentation (.md files)
├── environment.yml     # Conda environment description
├── requirements.txt    # Python dependencies
├── .gitignore
├── README.md
```

## Quickstart

1. **Install dependencies**
    ```bash
    conda env create -f environment.yml
    conda activate protein-domain-finder-2
    # or
    pip install -r requirements.txt
    ```
2. **Run the main workflow**
    Open `notebooks/protein_domain_workflow.ipynb` in Jupyter and follow the steps, OR run `src/main.py` if available.

3. **Data**
    - Place data files in the `data/` directory. See `data/README.md` for info.
    - Large data files are not stored in this repository.

## Documentation

- See `docs/` for evaluation, organization, and optimization plans.

## Scripts

- `scripts/complete_mmseqs2_workflow.sh` for end-to-end analysis pipeline.

## Contributing

Pull requests and suggestions welcome! Check issues for ideas.

## License

[Specify your license here]
