# Visium stitcher

The repository houses a rudimentary implementation of the Visium stitcher, allowing for a simple joining of adjacent tissue slices into a joint AnnData object for analysis. This is a two-step process, with the first step requiring aligning the images in Fiji. The resulting `.xml` file is subsequently used to transform the images and coordinates in Python.

Please refer to the [example notebook](notebooks/demo.ipynb) for exact Fiji instructions and a subsequent Python demonstration. The demo data can be found at `/lustre/scratch126/cellgen/team205/kp9/jupyter/ny1/23-03-06-stitcher/visium_stitcher`

## Installation

```bash
pip install git+https://github.com/Teichlab/visium_stitcher.git
```
