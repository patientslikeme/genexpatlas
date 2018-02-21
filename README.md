# genexpatlas
A package to programmatically access the [EMBL-EBI's Gene Expression Atlas](https://www.ebi.ac.uk/gxa/home) via Python

## Overview
This project began as a translation of the [ExpressionAtlas](http://www.bioconductor.org/packages/release/bioc/html/ExpressionAtlas.html)
R package into Python 2.x at [PatientsLikeMe](https://blog.patientslikeme.com). We are open-sourcing this tool to make it 
available to researchers, engineers, and individuals who wish to use Python to programmatically access the Gene 
Expression Atlas. 

A note on data formats: the original ExpressionAtlas package downloads an `.Rdata` summary file for research. This isn't
yet supported here, but what we do provide is a native-ish Python data structure for describing experiments in the GXA. 
Each experiment is represented by a dictionary (and collectively are in a list of dictionaries), which contains metadata 
describing the experiment itself, as well as data describing comparisons made between conditions, which are loaded into 
a Pandas dataframe.

**NOTE**: So far this has been mostly tested on comparison studies (e.g. [E-GEOD-10315](https://www.ebi.ac.uk/gxa/experiments/E-GEOD-10315/Results)) with `*analytics.tsv` files, 
support for baseline studies is, as of the current release, untested and the lack of an `*analytics.tsv` file will 
cause an error.

## Installation
Install with `pip` via Github:  
`pip install git+git://github.com/patientslikeme/genexpatlas.git`

## Contributing
We welcome any contributions (including but not limited to bug reports, bug fixes, documentation and enhancements) to 
ensure the robustness of this package and to ensure that it covers as many usecases as possible while still being 
concise and easy to use. To contribute, please make a fork of this project, develop there, and when ready and tested, make a pull request to this repo. 

Please check the repo's [issues](https://github.com/patientslikeme/genexpatlas/issues) for ideas on where to start, as 
well as for conversation on the package.

We strongly recommend setting up a virtual environment for your testing, and using `pip install -e` to install a
development executable.

## Usage
This package is a collection of methods that can be used on their own without any object instantiation.

```python
import genexpatlas as gea

# Search for data on humans containing the phrases 'term1' OR 'term2'
experiments = gea.search_atlas_experiments(search=['term1', 'term2'], species='homo sapiens')

# Pull experiment data into a list of dictionaries
loaded_data = gea.get_atlas_experiments(experiments)
```

## License
Like the ExpressionAtlas R package that this is based on, this software is licensed under [GPL v3](https://www.gnu.org/licenses/gpl-3.0.en.html).

## Future Directions
* Flexibility to process non-comparison data elegantly
* Support binary download of `.Rdata` summary files
* Support translation of `.Rdata` summary files
* Verify Python 3.x compatibility
* Add to PyPI