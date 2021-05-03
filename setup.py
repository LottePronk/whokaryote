import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="contig_classifier",  # Replace with your own username
    version="0.0.1",
    author="Lotte Pronk",
    author_email="lotte.pronk@wur.nl",
    description="A fast and easy tool to predict if a metagenomic contig is eukaryotic or prokaryotic",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://gitlab.com/",
    project_urls={
        "Bug Tracker": "https://gitlab.com",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Linux or MacOS",
    ],
    #  package_dir={"": "src"},
    #  packages=setuptools.find_packages(where="src"),
    packages=['contig_classifier'],
    include_package_data=True,
    python_requires=">=3.8",
    scripts=['bin/commandline.py'],
    install_requires=['pandas', 'pathlib', 'biopython', 'numpy', 'joblib', 'sklearn']
)
