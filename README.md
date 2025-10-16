# Visium stitcher

The repository houses a Visium stitcher, allowing for a simple joining of adjacent tissue slices into a joint AnnData object for analysis. This is a two-step process, with the first step requiring aligning the images and extracting the transformation matrices in Fiji ([TrakEM2](https://imagej.net/plugins/trakem2/)). The resulting `.xml` file is subsequently used to transform the images and coordinates in Python.

## Installation

The package is available on PyPi
```bash
pip install visium_stitcher
```

## Usage

Please refer to the [example notebook](https://nbviewer.org/github/Teichlab/visium_stitcher/blob/main/notebooks/demo.ipynb) for exact Fiji instructions and a subsequent Python demonstration. The demo data is available at `ftp://ftp.sanger.ac.uk/pub/users/kp9/demo_data.tar.gz`.

The function docstrings are provided on [ReadTheDocs](https://visium-stitcher.readthedocs.io/en/latest/).

## Citation

Please cite our [paper](https://www.nature.com/articles/s41586-023-06806-x)
```
@article{zhang2024human,
  title={A human embryonic limb cell atlas resolved in space and time},
  author={Zhang, Bao and He, Peng and Lawrence, John EG and Wang, Shuaiyu and Tuck, Elizabeth and Williams, Brian A and Roberts, Kenny and Kleshchevnikov, Vitalii and Mamanova, Lira and Bolt, Liam and others},
  journal={Nature},
  volume={635},
  number={8039},
  pages={668--678},
  year={2024},
  publisher={Nature Publishing Group UK London}
}
```
