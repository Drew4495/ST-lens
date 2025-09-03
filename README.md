# fPCA for Spatial Transcriptomics

Functional Principal Component Analysis (fPCA) is a fast and accurate method for spatial domain detection in spatial transcriptomics. By modeling gene expression patterns as smooth functions over a mesh representation of tissue, fPCA better captures complex geometries (including holes and non-convexities) than methods relying on Euclidean distance or graph structures. Compared to leading algorithms, it offers greater interpretability through a single smoothing parameter, reduced runtime, and improved accuracy across both simulated and real datasets (e.g., HER2 breast cancer, human cortex, mouse cerebellum). These features make fPCA a scalable and robust tool for exploring tissue organization and gene expression architecture

![fPCA Overview](README_images/ST_Method_Overview_Vertical.png)

A non-convexity (concavity) ....
![Convexity and Non-Convexity Example](README_images/convex.png)

Photo adapted from https://undergroundmathematics.org/glossary/convex-shape




## Tutorials

Tutorials created by Marco Galliani (acknowledged in manuscript)
- [fPCA 2D Tutorial (code on GitHub)](https://github.com/Drew4495/ST-lens/blob/main/Tutorials/fPCA_2D.html)  
- [fPCA 2D Tutorial (rendered)](https://htmlpreview.github.io/?https://raw.githubusercontent.com/Drew4495/ST-lens/main/Tutorials/fPCA_2D.html)




## Manuscript scripts

- Talk about manuscript (in review). Breifly talk about 3 different tissues
- For full scripts that recreated the manuscript, see additional readme.md in ???? and the scripts in 











TEMPLATE
# Foobar

Foobar is a Python library for dealing with word pluralization.

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install foobar.

```bash
pip install foobar
```

## Usage

```python
import foobar

# returns 'words'
foobar.pluralize('word')

# returns 'geese'
foobar.pluralize('goose')

# returns 'phenomenon'
foobar.singularize('phenomena')
```

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

[MIT](https://choosealicense.com/licenses/mit/)
