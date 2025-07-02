from setuptools import dist, setup, find_packages, Extension


def readme():
    with open("README.md") as f:
        return f.read()


setup(
    name="vc",
    version="0.1",
    author="Caoqi Fan",
    author_email="fcq1116@gmail.com",
    description="Virtual Cells",
    long_description=readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/Ephraim-usc/vc.git",
    packages=find_packages(),
    python_requires=">=3",
    install_requires=[
        "numpy>=1.14.5",
        "pandas",
    ],
    ext_modules=[],
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
