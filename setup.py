import setuptools
from lec import __version__

with open("readme.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="lec",
    version=__version__,
    description="Local Three-Dimensional Magnetostatic Equilibrium Code",
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering",
    ],
    # url="",
    author="Ke Liu",
    author_email="lk2020@mail.ustc.edu.cn",
    license="GNU 3.0",
    packages=setuptools.find_packages(),
)

