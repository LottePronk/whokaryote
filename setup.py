import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Whokaryote",
    version="0.0.1",
    author="Lotte Pronk",
    author_email="lotte.pronk@wur.nl",
    description="Whokaryote: A fast and easy tool to predict if a metagenomic contig is eukaryotic or prokaryotic",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://git.wageningenur.nl/lotte.pronk/whokaryote",
    project_urls={
        "Bug Tracker": "https://git.wageningenur.nl/lotte.pronk/whokaryote",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: GNU Affero General Public License v3.0",
        "Operating System :: Linux or MacOS",
    ],
    packages=['whokaryote_scripts'],
    include_package_data=True,
    package_data={'whokaryote_scripts': ['data/*.joblib']},
    python_requires=">=3.8",
    scripts=['bin/whokaryote.py'],
    install_requires=['pandas', 'pathlib', 'biopython', 'numpy', 'joblib', 'scikit-learn==1.0.2']
)
